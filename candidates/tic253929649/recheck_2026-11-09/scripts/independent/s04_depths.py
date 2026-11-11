"""Per-series and pooled box depth at the claimed linear ephemeris, with empirical nulls.
Null 1: same P, phase shifted over a grid (step T14/4), excluding |shift| < 2*T14.
Null 2: random ephemerides, P' ~ U(0.7P, 1.4P), random phase (coherent across all series).
"""
import sys, json
import numpy as np, pandas as pd
from core import *

rng = np.random.default_rng(20261109)
CONFIGS = {
    "A_lkdef_runmed":   dict(mask="lkdef",  flux="auto", detrend="runmed", window=0.75),
    "B_strict_runmed":  dict(mask="strict", flux="auto", detrend="runmed", window=0.75),
    "C_strict_pipedet": dict(mask="strict", flux="det",  detrend="runmed", window=0.75),
    "D_strict_spline":  dict(mask="strict", flux="auto", detrend="spline", window=0.75),
    "F_lkdef_pipedet":  dict(mask="lkdef",  flux="det",  detrend="runmed", window=0.75),
    "G_strict_runmed1.5": dict(mask="strict", flux="auto", detrend="runmed", window=1.5),
    "H_strict_runmed0.4": dict(mask="strict", flux="auto", detrend="runmed", window=0.4),
}
PRIORITY = {"SPOC": 0, "TESS-SPOC": 1, "QLP": 2}


def sums(t, r, P, T0, T14):
    m = in_transit(t, P, T0, T14)
    return r[m].sum(), int(m.sum())


def pooled(series, P, T0, T14):
    num = den = 0.0
    for s in series:
        sm, n = sums(s["t"], s["r"], P, T0, T14)
        w = 1.0 / s["sigma"] ** 2
        num += w * sm
        den += w * n
    return (-num / den if den > 0 else np.nan), den


def null_draws(P, T0, T14, nrand=4000):
    k = np.arange(1, int(P / (T14 / 4)))
    shifts = k * T14 / 4
    shifts = shifts[(shifts > 2 * T14) & (shifts < P - 2 * T14)]
    grid = [(P, T0 + sh) for sh in shifts]
    rand = [(p, T0 + u * p) for p, u in zip(rng.uniform(0.7 * P, 1.4 * P, nrand), rng.uniform(0, 1, nrand))]
    return grid, rand


def pval(obs, null):
    null = np.asarray(null); null = null[np.isfinite(null)]
    return (1 + np.sum(null >= obs)) / (1 + null.size), np.nanstd(null), null.size


rows = []
for tic, e in TARGETS.items():
    P, T0, T14 = e["P"], e["T0"], e["T14"]
    raw = load_raw(tic)
    grid, rand = null_draws(P, T0, T14)
    for cname, cfg in CONFIGS.items():
        series = []
        for s in raw:
            pr = prepare(s, **cfg)
            if pr is None: continue
            pr.update(author=s["author"], sector=s["sector"], name=label(s))
            series.append(pr)
        print(f"\n##### TIC {tic}  config {cname}")
        print(f"{'series':15s} {'n':>5s} {'sig10m':>7s} {'n_in':>4s} {'depth':>7s} {'formal':>6s} {'nullSD_ph':>9s} {'p_phase':>8s} {'nullSD_rnd':>10s} {'p_rand':>8s}")
        for s in series:
            d, err, n = box_depth(s["t"], s["r"], P, T0, T14, s["sigma"])
            ng = [box_depth(s["t"], s["r"], p, t0, T14, s["sigma"])[0] for p, t0 in grid]
            nr = [box_depth(s["t"], s["r"], p, t0, T14, s["sigma"])[0] for p, t0 in rand]
            pg, sg, _ = pval(d, ng); pr_, sr, _ = pval(d, nr)
            print(f"{s['name']:15s} {s['t'].size:5d} {s['sigma']*1e6:7.0f} {n:4d} {d*1e6:7.0f} {err*1e6:6.0f} {sg*1e6:9.0f} {pg:8.4f} {sr*1e6:10.0f} {pr_:8.4f}")
            rows.append(dict(tic=tic, config=cname, series=s["name"], n=s["t"].size, sigma_ppm=s["sigma"]*1e6, n_in=n,
                             depth_ppm=d*1e6, formal_ppm=err*1e6, nullsd_phase_ppm=sg*1e6, p_phase=pg, nullsd_rand_ppm=sr*1e6, p_rand=pr_))
        groups = {
            "POOLED all series": series,
            "POOLED one per sector (SPOC>TESS-SPOC>QLP)": [min([x for x in series if x["sector"] == sec], key=lambda x: PRIORITY[x["author"]]) for sec in sorted({x["sector"] for x in series})],
            "POOLED one per sector (QLP>TESS-SPOC>SPOC)": [max([x for x in series if x["sector"] == sec], key=lambda x: PRIORITY[x["author"]]) for sec in sorted({x["sector"] for x in series})],
            "POOLED SPOC 2-min only": [x for x in series if x["author"] == "SPOC"],
            "POOLED TESS-SPOC only": [x for x in series if x["author"] == "TESS-SPOC"],
            "POOLED QLP only": [x for x in series if x["author"] == "QLP"],
        }
        for gname, g in groups.items():
            if not g: continue
            d, den = pooled(g, P, T0, T14)
            formal = 1 / np.sqrt(den)
            ng = [pooled(g, p, t0, T14)[0] for p, t0 in grid]
            nr = [pooled(g, p, t0, T14)[0] for p, t0 in rand]
            pg, sg, n1 = pval(d, ng); pr_, sr, n2 = pval(d, nr)
            print(f"{gname:45s} [{'+'.join(x['name'] for x in g)}]")
            print(f"     depth={d*1e6:6.0f} ppm formal={formal*1e6:5.0f} | phase-null SD={sg*1e6:5.0f} (N={n1}) z={d/sg:5.2f} p={pg:.4f} | random-ephemeris null SD={sr*1e6:5.0f} (N={n2}) z={d/sr:5.2f} p={pr_:.5f}")
            rows.append(dict(tic=tic, config=cname, series=gname, n=sum(x["t"].size for x in g), sigma_ppm=np.nan, n_in=np.nan,
                             depth_ppm=d*1e6, formal_ppm=formal*1e6, nullsd_phase_ppm=sg*1e6, p_phase=pg, nullsd_rand_ppm=sr*1e6, p_rand=pr_))
pd.DataFrame(rows).to_csv("out_depths.csv", index=False)
