"""Own aperture photometry from a TESScut cutout (written from scratch)."""
import numpy as np
from astropy.io import fits
from core import robust_sigma

def cut_phot(path, half=1, strict=True, verbose=False):
    with fits.open(path) as h:
        d = h[1].data
        t = np.array(d["TIME"], float); q = np.array(d["QUALITY"]).astype(np.int64)
        fl = np.array(d["FLUX"], float)
    med = np.nanmedian(fl, axis=0)
    # target = brightest pixel within the central 5x5
    sub = med[3:8, 3:8]; iy, ix = np.unravel_index(np.nanargmax(sub), sub.shape); iy += 3; ix += 3
    ap = np.zeros(med.shape, bool); ap[iy - half:iy + half + 1, ix - half:ix + half + 1] = True
    # background pixels: faintest 20% of the cutout in the median image, outside the aperture
    thr = np.nanpercentile(med[~ap], 20)
    bg = (med <= thr) & ~ap
    bkg = np.nanmedian(fl[:, bg], axis=1)
    flux = np.nansum(fl[:, ap], axis=1) - ap.sum() * bkg
    # flux-weighted centroid in the aperture
    yy, xx = np.mgrid[0:med.shape[0], 0:med.shape[1]]
    net = fl - bkg[:, None, None]
    cx = np.nansum(net[:, ap] * xx[ap], axis=1) / np.nansum(net[:, ap], axis=1)
    cy = np.nansum(net[:, ap] * yy[ap], axis=1) / np.nansum(net[:, ap], axis=1)
    ok = np.isfinite(t) & np.isfinite(flux) & (flux > 0)
    if strict: ok &= (q == 0)
    else: ok &= (q & (1 | 2 | 4 | 8 | 16 | 32 | 128 | 16384)) == 0
    if verbose:
        print(f"  target pixel (row,col)=({iy},{ix}); aperture {2*half+1}x{2*half+1}; n_bg_pix={bg.sum()}; median net flux={np.nanmedian(flux[ok]):.0f} e/s; bkg median={np.nanmedian(bkg[ok]):.1f} e/s/px (range {np.nanmin(bkg[ok]):.0f}-{np.nanmax(bkg[ok]):.0f})")
    return dict(t=t[ok], f=flux[ok], bkg=bkg[ok], cx=cx[ok], cy=cy[ok], q=q[ok], med=med, target=(iy, ix))
