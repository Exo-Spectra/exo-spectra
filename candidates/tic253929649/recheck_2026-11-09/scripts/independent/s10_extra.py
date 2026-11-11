"""Extra checks: local-linear-baseline statistic (alternative detrending) with nulls, per-event depths,
leave-one-sector-out, odd/even, comparison with the authors' depth_by_pipeline.csv at THEIR BLS ephemeris,
band-limited box search rank."""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from astropy.timeseries import BoxLeastSquares
from core import *
rng = np.random.default_rng(11)
PRIORITY = {"SPOC": 0, "TESS-SPOC": 1, "QLP": 2}
def prep_all(raw, cfg):
    out = []
    for s in raw:
        pr = prepare(s, **cfg)
        if pr is None: continue
        pr.update(author=s["author"], sector=s["sector"], name=label(s)); out.append(pr)
    return out
def ops(series):
    return [min([x for x in series if x["sector"] == sec], key=lambda x: PRIORITY[x["author"]]) for sec in sorted({x["sector"] for x in series})]
def pooled(series, P, T0, T14):
    num = den = 0.0
    for s in series:
        m = in_transit(s["t"], P, T0, T14); w = 1 / s["sigma"] ** 2
        num += w * s["r"][m].sum(); den += w * m.sum()
    return (-num / den if den > 0 else np.nan)
def pooled_local(series, P, T0, T14):
    num = den = 0.0
    for s in series:
        d, err, nin, nev = local_depth(s["t"], s["f"], P, T0, T14, sigma=s["sigma"])
        if nin == 0: continue
        w = nin / s["sigma"] ** 2; num += w * d; den += w
    return num / den if den > 0 else np.nan
def pv(o, null):
    null = np.asarray(null); null = null[np.isfinite(null)]
    return (1 + np.sum(null >= o)) / (1 + null.size), null.std()

AUTH = {253929649: dict(P=4.878611, T0=1329.63239, T14=0.96 / 24, csv="E:/projects/exo-spectra/reports/tess_hunt/deep_253929649/depth_by_pipeline.csv"),
        142825877: dict(P=7.521091, T0=1684.37428, T14=1.92 / 24, csv="E:/projects/exo-spectra/reports/tess_hunt/deep_142825877/depth_by_pipeline.csv")}
for tic, e in TARGETS.items():
    P, T0, T14 = e["P"], e["T0"], e["T14"]
    raw = load_raw(tic)
    print(f"\n{'='*20} TIC {tic} {'='*20}")
    # --- 1. local-linear statistic on un-detrended flux (strict mask; PDCSAP / QLP pipeline-detrended flux)
    loc = prep_all(raw, dict(mask="strict", flux="det", detrend="none"))
    rand = [(p, T0 + u * p) for p, u in zip(rng.uniform(0.7 * P, 1.4 * P, 1500), rng.uniform(0, 1, 1500))]
    sh = np.arange(T14 / 4, P, T14 / 4); sh = sh[(sh > 2 * T14) & (sh < P - 2 * T14)]
    for gname, g in (("all series", loc), ("one per sector (SPOC>TESS-SPOC>QLP)", ops(loc)), ("QLP only", [x for x in loc if x["author"] == "QLP"])):
        d = pooled_local(g, P, T0, T14)
        p1, s1 = pv(d, [pooled_local(g, P, T0 + x, T14) for x in sh]); p2, s2 = pv(d, [pooled_local(g, p, t, T14) for p, t in rand])
        print(f"  [local-linear baseline, strict] {gname:38s}: depth={d*1e6:5.0f} ppm | phase null SD={s1*1e6:4.0f} z={d/s1:5.2f} p={p1:.4f} | random-ephemeris null SD={s2*1e6:4.0f} z={d/s2:5.2f} p={p2:.4f}")
    # --- 2. leave-one-sector-out and odd/even (strict, running median, one per sector)
    B = prep_all(raw, dict(mask="strict", flux="auto", detrend="runmed", window=0.75)); g = ops(B)
    d = pooled(g, P, T0, T14); p2, s2 = pv(d, [pooled(g, p, t, T14) for p, t in rand])
    print(f"  [running median, strict] one per sector: depth={d*1e6:.0f} ppm z={d/s2:.2f} p={p2:.4f}")
    for x in g:
        gg = [y for y in g if y is not x]
        d = pooled(gg, P, T0, T14); p2, s2 = pv(d, [pooled(gg, p, t, T14) for p, t in rand])
        print(f"     without sector {x['sector']:3d}: depth={d*1e6:5.0f} ppm z={d/s2:5.2f} p={p2:.4f}")
    for par, nm in ((0, "even"), (1, "odd")):
        # odd/even = ephemeris with 2P anchored at T0 (+P)
        d = pooled(g, 2 * P, T0 + par * P, T14); p2, s2 = pv(d, [pooled(g, 2 * p, t + par * p, T14) for p, t in rand])
        print(f"     {nm} events only: depth={d*1e6:5.0f} ppm (null SD {s2*1e6:.0f}, z={d/s2:.2f})")
    # --- 3. per-event depths (one per sector)
    print("  per-event depths (one per sector, strict, running median; local-linear in brackets):")
    lg = ops(loc)
    for x, xl in zip(g, lg):
        n0 = int(np.ceil((x["t"][0] - T0) / P)); n1 = int(np.floor((x["t"][-1] - T0) / P))
        for n in range(n0, n1 + 1):
            tc = T0 + n * P; m = np.abs(x["t"] - tc) < T14 / 2
            if m.sum() < 2: 
                print(f"     {x['name']:14s} cycle {n:4d} tc={tc:9.3f}: n_in={int(m.sum())} (gap)"); continue
            dep = -x["r"][m].mean(); err = x["sigma"] / np.sqrt(m.sum())
            ml = np.abs(xl["t"] - tc) < 3 * T14
            dl = local_depth(xl["t"][ml], xl["f"][ml], P, tc, T14, sigma=xl["sigma"])[0] if ml.sum() > 8 else np.nan
            print(f"     {x['name']:14s} cycle {n:4d} tc={tc:9.3f}: n_in={int(m.sum()):2d} depth={dep*1e6:6.0f} +- {err*1e6:4.0f} ppm  [{dl*1e6:6.0f}]")
    # --- 4. authors' depth_by_pipeline.csv at their BLS ephemeris
    a = AUTH[tic]; df = pd.read_csv(a["csv"])
    F = {x["name"]: x for x in prep_all(raw, dict(mask="lkdef", flux="det", detrend="runmed", window=0.75))}
    Bd = {x["name"]: x for x in B}
    print(f"  comparison with authors' depth_by_pipeline.csv at their BLS ephemeris P={a['P']}, t0={a['T0']}, dur={a['T14']*24:.2f} h:")
    print("     series          authors(depth+-err, n_in) | mine lkdef+pipeline flux | mine strict+own detrend")
    for _, r in df.iterrows():
        nm = f"{r['author']}-s{int(r['sector'])}"
        o = []
        for D in (F, Bd):
            if nm in D:
                d, err, n = box_depth(D[nm]["t"], D[nm]["r"], a["P"], a["T0"], a["T14"], D[nm]["sigma"]); o.append(f"{d*1e6:6.0f}+-{err*1e6:3.0f} (n={n:2d})")
            else: o.append("   n/a")
        print(f"     {nm:15s} {r['depth_ppm']:6.1f}+-{r['err_ppm']:5.1f} (n={int(r['n_in']):2d}) | {o[0]} | {o[1]}")
    # --- 5. band-limited box search (0.75P..1.33P), one per sector strict
    t = np.concatenate([x["t"] for x in g]); r = np.concatenate([x["r"] for x in g]); dy = np.concatenate([np.full(x["t"].size, x["sigma"]) for x in g])
    o = np.argsort(t); t, r, dy = t[o], r[o], dy[o]
    c = T14 / (4.0 * (t[-1] - t[0])); periods = np.exp(np.arange(np.log(0.75 * P), np.log(1.33 * P), c))
    out = {}
    for sign in (1, -1):
        pw = BoxLeastSquares(t, sign * r, dy).power(periods, [T14], objective="snr", method="fast", oversample=8)
        out[sign] = np.asarray(pw.depth_snr)
    near = np.abs(periods - P) < 0.002; cl = out[1][near].max()
    # count distinct peaks (local maxima separated by >0.5% in period) above the claimed one
    def npeaks(snr, thr):
        idx = np.where(snr > thr)[0]; 
        if idx.size == 0: return 0
        groups = 1 + int(np.sum(np.diff(periods[idx]) / periods[idx][:-1] > 0.002)); return groups
    print(f"  band-limited box search {0.75*P:.2f}-{1.33*P:.2f} d, fixed T14 (one per sector, strict): claimed-period SNR={cl:.2f};"
          f" max dip SNR in band={out[1].max():.2f} at P={periods[np.argmax(out[1])]:.5f}; max bump SNR in band={out[-1].max():.2f};"
          f" distinct dip peaks above claimed: {npeaks(out[1], cl)}; distinct bump peaks above claimed: {npeaks(out[-1], cl)}; n_periods={periods.size}")
