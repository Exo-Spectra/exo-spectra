"""Independent re-check library (written from scratch, no project code used).

Loading of SPOC / TESS-SPOC / QLP light-curve FITS files with astropy only,
quality masking, binning, detrending and the box-depth statistic.
"""
import glob
import os
import numpy as np
from astropy.io import fits

ROOT = "E:/projects/exo-spectra/data/tess_hunt/qlp/mastDownload"
LK_DEFAULT = 17087  # lightkurve 2.6.0 TessQualityFlags.DEFAULT_BITMASK

TARGETS = {
    253929649: dict(P=4.814602, T0=1329.651, T14=1.1 / 24.0, depth=432e-6),
    142825877: dict(P=7.521089, T0=1684.36146, T14=2.0 / 24.0, depth=407e-6),
}


def robust_sigma(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 3:
        return np.nan
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def list_files(tic):
    out = []
    for f in sorted(glob.glob(f"{ROOT}/TESS/*{tic}*/*_lc.fits")):
        out.append(("SPOC", f))
    for f in sorted(glob.glob(f"{ROOT}/HLSP/hlsp_tess-spoc*{tic}*/*_lc.fits")):
        out.append(("TESS-SPOC", f))
    for f in sorted(glob.glob(f"{ROOT}/HLSP/hlsp_qlp*{tic}*/*_llc.fits")):
        out.append(("QLP", f))
    return out


def load_raw(tic):
    """Return list of dicts with raw columns (no masking)."""
    res = []
    for author, f in list_files(tic):
        with fits.open(f) as h:
            p = h[0].header
            d = h[1].data
            names = [c.name for c in h[1].columns]
            t = np.array(d["TIME"], float)
            q = np.array(d["QUALITY"]).astype(np.int64)
            sap = np.array(d["SAP_FLUX"], float)
            bkg = np.array(d["SAP_BKG"], float)
            if author == "QLP":
                col = "KSPSAP_FLUX" if "KSPSAP_FLUX" in names else "DET_FLUX"
                det = np.array(d[col], float)
                orbit = np.array(d["ORBITID"]).astype(int)
            else:
                col = "PDCSAP_FLUX"
                det = np.array(d[col], float)
                orbit = np.zeros(t.size, int)
            res.append(dict(author=author, sector=int(p["SECTOR"]), file=os.path.basename(f),
                            t=t, q=q, sap=sap, det=det, bkg=bkg, detcol=col, orbit=orbit,
                            cam=p.get("CAMERA"), ccd=p.get("CCD")))
    return res


def quality_keep(q, mode):
    if mode == "lkdef":
        return (q & LK_DEFAULT) == 0
    if mode == "strict":
        return q == 0
    if mode == "none":
        return np.ones(q.size, bool)
    raise ValueError(mode)


def bin_series(t, f, width=600.0 / 86400.0, native_dt=None):
    """Bin to `width` on a fixed grid; keep bins with >= half nominal occupancy."""
    if native_dt is None:
        native_dt = np.nanmedian(np.diff(t))
    if native_dt >= 0.9 * width:
        return t.copy(), f.copy()
    nominal = width / native_dt
    idx = np.floor((t - t[0]) / width).astype(np.int64)
    n = np.bincount(idx)
    st = np.bincount(idx, weights=t)
    sf = np.bincount(idx, weights=f)
    ok = n >= max(1.0, 0.5 * nominal)
    return st[ok] / n[ok], sf[ok] / n[ok]


def running_median(t, f, window):
    """Time-windowed running median, never bridging gaps > 0.5 d."""
    out = np.empty_like(f)
    seg = np.concatenate([[0], np.where(np.diff(t) > 0.5)[0] + 1, [t.size]])
    for a, b in zip(seg[:-1], seg[1:]):
        tt, ff = t[a:b], f[a:b]
        lo = np.searchsorted(tt, tt - window / 2.0, side="left")
        hi = np.searchsorted(tt, tt + window / 2.0, side="right")
        for i in range(tt.size):
            out[a + i] = np.median(ff[lo[i]:hi[i]])
    return out


def spline_trend(t, f, knot=0.75):
    """Least-squares cubic spline with knot spacing `knot` d per segment (gaps > 0.5 d)."""
    from scipy.interpolate import LSQUnivariateSpline
    out = np.empty_like(f)
    seg = np.concatenate([[0], np.where(np.diff(t) > 0.5)[0] + 1, [t.size]])
    for a, b in zip(seg[:-1], seg[1:]):
        tt, ff = t[a:b], f[a:b]
        span = tt[-1] - tt[0]
        if tt.size < 12 or span < 2 * knot:
            out[a:b] = np.polyval(np.polyfit(tt - tt[0], ff, 1 if tt.size > 3 else 0), tt - tt[0])
            continue
        nk = int(span // knot)
        knots = np.linspace(tt[0], tt[-1], nk + 1)[1:-1]
        try:
            out[a:b] = LSQUnivariateSpline(tt, ff, knots, k=3)(tt)
        except Exception:
            out[a:b] = running_median(tt, ff, knot)
    return out


def prepare(s, mask="strict", flux="auto", detrend="runmed", window=0.75, clip=5.0, bin_width=600.0 / 86400.0):
    """Mask -> bin (to 10 min) -> normalise -> detrend -> symmetric sigma clip.

    flux: 'auto' = PDCSAP for SPOC/TESS-SPOC and SAP for QLP (own detrending of raw photometry);
          'det'  = pipeline flux (PDCSAP; QLP KSPSAP/DET);
    detrend: 'runmed' | 'spline' | 'none' (only median normalisation; used by the local-baseline statistic)
    Returns dict(t, r, sigma) with r = residual relative flux (zero median).
    """
    if flux == "auto":
        f = s["sap"] if s["author"] == "QLP" else s["det"]
    elif flux == "det":
        f = s["det"]
    elif flux == "sap":
        f = s["sap"]
    else:
        raise ValueError(flux)
    keep = quality_keep(s["q"], mask) & np.isfinite(f) & np.isfinite(s["t"])
    t, f = s["t"][keep], f[keep]
    if t.size < 50:
        return None
    o = np.argsort(t)
    t, f = t[o], f[o]
    native = np.nanmedian(np.diff(s["t"][np.isfinite(s["t"])]))
    t, f = bin_series(t, f, bin_width, native)
    f = f / np.nanmedian(f)
    good = np.ones(t.size, bool)
    r = f - 1.0
    for _ in range(2):
        if detrend == "runmed":
            tr = running_median(t[good], f[good], window)
        elif detrend == "spline":
            tr = spline_trend(t[good], f[good], window)
        elif detrend == "none":
            tr = np.full(good.sum(), np.median(f[good]))
        else:
            raise ValueError(detrend)
        r_good = f[good] / tr - 1.0
        sig = robust_sigma(r_good)
        r = np.full(t.size, np.nan)
        r[good] = r_good
        if clip is None:
            break
        newgood = good & (np.abs(np.nan_to_num(r, nan=1e9)) < clip * sig)
        if newgood.sum() == good.sum():
            break
        good = newgood
    ok = np.isfinite(r) & good
    t, r = t[ok], r[ok]
    r = r - np.median(r)
    return dict(t=t, r=r, sigma=robust_sigma(r), f=f[ok])


def in_transit(t, P, T0, T14):
    ph = (t - T0 + 0.5 * P) % P - 0.5 * P
    return np.abs(ph) < 0.5 * T14


def box_depth(t, r, P, T0, T14, sigma=None):
    """Depth = mean(out) - mean(in) on residuals. Returns depth, formal err, n_in."""
    m = in_transit(t, P, T0, T14)
    n = int(m.sum())
    if n < 2:
        return np.nan, np.nan, n
    if sigma is None:
        sigma = robust_sigma(r)
    d = np.mean(r[~m]) - np.mean(r[m])
    return d, sigma / np.sqrt(n), n


def local_depth(t, f, P, T0, T14, half=3.0, sigma=None):
    """Alternative 'detrending': per-event straight-line baseline fitted to out-of-transit
    points within +-half*T14 of each predicted mid-time; depth = -mean in-transit residual.
    Events need >=2 in-transit points and >=3 out-of-transit points on each side.
    Returns weighted depth, formal err, n_in(total), n_events."""
    if sigma is None:
        sigma = robust_sigma(np.diff(f)) / np.sqrt(2.0)
    n0 = int(np.floor((t[0] - T0) / P)) - 1
    n1 = int(np.ceil((t[-1] - T0) / P)) + 1
    num = den = 0.0
    nin = nev = 0
    for n in range(n0, n1 + 1):
        tc = T0 + n * P
        a = np.searchsorted(t, tc - half * T14)
        b = np.searchsorted(t, tc + half * T14)
        if b - a < 8:
            continue
        tt, ff = t[a:b] - tc, f[a:b]
        m = np.abs(tt) < 0.5 * T14
        left = (~m) & (tt < 0)
        right = (~m) & (tt > 0)
        if m.sum() < 2 or left.sum() < 3 or right.sum() < 3:
            continue
        c = np.polyfit(tt[~m], ff[~m], 1)
        d = -np.mean(ff[m] - np.polyval(c, tt[m]))
        w = m.sum()
        num += w * d
        den += w
        nin += int(m.sum())
        nev += 1
    if den == 0:
        return np.nan, np.nan, 0, 0
    return num / den, sigma / np.sqrt(nin), nin, nev


def label(s):
    return f"{s['author']}-s{s['sector']}"
