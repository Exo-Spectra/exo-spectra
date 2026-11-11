"""s106 alias test replicated with the anchor the authors' own tables imply: T0 = 1329.63239 (BLS t0,
'2018-anchored epoch') for all five periods. (The reference of timing_per_series.csv is P=4.814600,
T0=1329.63239 - derived from T_obs - dt of all six rows.)"""
import numpy as np, warnings; warnings.filterwarnings("ignore")
from core import *
from tesscut_phot import cut_phot
from astropy.io import fits
CUT = "E:/projects/exo-spectra/data/tess_hunt/tesscut/tesscut/tess-s0106-1-3_338.369633_-35.683456_11x11_astrocut.fits"
T14 = 1.1 / 24; T0A = 1329.63239
ALIASES = [4.878611, 4.852472, 4.814564, 4.869591, 4.789319]
AUTH = {4.878611: (-180.9, 146.0), 4.852472: (-116.7, 156.1), 4.814564: (368.0, 144.2), 4.869591: (25.8, 146.0), 4.789319: (140.6, 146.0)}
def prep_tf(t, f, window=0.75, clip=5.0):
    o = np.argsort(t); t, f = t[o], f[o]
    t, f = bin_series(t, f, 600 / 86400, np.median(np.diff(t)))
    f = f / np.median(f); good = np.ones(t.size, bool)
    for _ in range(2):
        tr = running_median(t[good], f[good], window)
        r = np.full(t.size, np.nan); r[good] = f[good] / tr - 1
        sig = robust_sigma(r[good]); good = good & (np.abs(np.nan_to_num(r, nan=9)) < clip * sig)
    r = r[good] - np.median(r[good])
    return dict(t=t[good], r=r, sigma=robust_sigma(r))
def centered_phot():
    # variant: 3x3 aperture on the geometric cutout centre (5,5) instead of the brightest pixel (6,5)
    with fits.open(CUT) as h:
        d = h[1].data; t = np.array(d["TIME"], float); q = np.array(d["QUALITY"]).astype(np.int64); fl = np.array(d["FLUX"], float)
    med = np.nanmedian(fl, axis=0); ap = np.zeros(med.shape, bool); ap[4:7, 4:7] = True
    bg = (med <= np.nanpercentile(med[~ap], 20)) & ~ap
    flux = np.nansum(fl[:, ap], axis=1) - ap.sum() * np.nanmedian(fl[:, bg], axis=1)
    ok = np.isfinite(t) & np.isfinite(flux) & (flux > 0) & (q == 0)
    return t[ok], flux[ok]
variants = {}
for apn, half, strict in (("3x3 on star, q==0", 1, True), ("3x3 on star, keep flag 2048", 1, False), ("5x5 on star, q==0", 2, True)):
    p = cut_phot(CUT, half, strict)
    for w in (0.75, 1.5):
        variants[f"{apn}, runmed {w} d"] = prep_tf(p["t"], p["f"], w)
t, f = centered_phot(); variants["3x3 on cutout centre (5,5), q==0, runmed 0.75 d"] = prep_tf(t, f, 0.75)
print("authors (alias_test_s106.csv): " + "; ".join(f"P={P}: {a[0]:+.0f}+-{a[1]:.0f}" for P, a in AUTH.items()))
allp = []
for vn, s in variants.items():
    print(f"\n[{vn}] n={s['t'].size} sigma(10 min)={s['sigma']*1e6:.0f} ppm")
    zs = []
    for P in ALIASES:
        d, err, n = box_depth(s["t"], s["r"], P, T0A, T14, s["sigma"])
        sh = np.arange(T14 / 4, P, T14 / 4); sh = sh[(sh > 2 * T14) & (sh < P - 2 * T14)]
        null = np.array([box_depth(s["t"], s["r"], P, T0A + x, T14, s["sigma"])[0] for x in sh]); null = null[np.isfinite(null)]
        p1 = (1 + np.sum(null >= d)) / (1 + null.size); zs.append((d / null.std(), p1, P, null))
        print(f"   P={P:.6f}: depth@predicted={d*1e6:6.0f} +- {err*1e6:3.0f} formal (n_in={n}); phase-null SD={null.std()*1e6:3.0f} ppm; z_emp={d/null.std():5.2f}; single-trial p={p1:.3f}")
    # best-of-five null: draw 5 independent phases (one per period), take the max z; repeat
    rng = np.random.default_rng(1); sims = []
    for _ in range(20000):
        sims.append(max(rng.choice(z[3]) / z[3].std() for z in zs))
    zobs = max(z[0] for z in zs); pbest = np.mean(np.array(sims) >= zobs)
    z4814 = [z for z in zs if z[2] == 4.814564][0]
    print(f"   best of five: P={max(zs)[2]} z_emp={zobs:.2f}; chance that the best of five random-phase periods reaches this = {pbest:.3f}"
          f"   | P=4.814564 alone: z_emp={z4814[0]:.2f}, single p={z4814[1]:.3f}, 'any of five' p={1-(1-z4814[1])**5:.3f}")
    allp.append((vn, z4814[0], z4814[1], 1 - (1 - z4814[1]) ** 5))
print("\nGaussian reading of the authors' own number: z=2.55 -> one-sided p=0.0054; any of five: %.3f" % (1 - (1 - 0.0054) ** 5))
