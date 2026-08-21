"""Limb-asymmetry study — targeted emcee escalation for >2 sigma rows.

Stage 2 ran without --mcmc, so the promised escalation for |dd| > 2 sigma
candidates never fired. This script re-fits those rows from the persisted
light curves (no raw-archive re-extraction) and adds dd_mcmc_ppm /
dd_mcmc_err_ppm columns to the catalog in place.

The catalog does not store the ramp coefficients, so the refit recovers
them from the `ramp` column persisted alongside the light curve (noiseless
-> exact); starting the ramp from zero made the first-pass refits drift
into other local minima (quad_exp settling exponential is degenerate with
the ingress depth), which showed up even on dd~0 control rows.

Usage: python limb_asym_mcmc.py [--thr SIGMA] [--planet NAME] [--all]
                                [--refit-only]
Run after: limb_asym_run.py
"""
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import limb_asym_model as lam                                # noqa: E402
from survey_analyze import ROOT                              # noqa: E402

OUT = ROOT / "reports" / "limb_asymmetry"
LC_DIR = OUT / "lightcurves"


def recover_ramp_coeffs(lc, cfg):
    """The persisted `ramp` column is the stage-2 ramp model evaluated on
    the (clipped) time grid — noiseless, so a least-squares fit recovers
    the coefficients exactly up to the x-conditioning of THIS cfg."""
    from scipy.optimize import least_squares
    t = lc.t_bmjd.to_numpy()
    y = lc.ramp.to_numpy()
    x = (t - cfg.t_ref) / cfg.span
    xe = t - cfg.t_min
    npar = lam.RAMP_NPAR[cfg.ramp]

    def solve(c0):
        return least_squares(
            lambda c: lam.ramp_flux(cfg.ramp, c, x, xe) - y, c0,
            method="lm" if cfg.ramp != "quad_exp" else "trf")

    if cfg.ramp != "quad_exp":
        sol = solve(np.zeros(npar))
    else:
        # the settling timescale creates local minima — multi-start over tau
        best = None
        for tau in np.geomspace(1e-3, 0.4, 12):
            c0 = np.zeros(npar)
            c0[3] = tau
            s = solve(c0)
            if best is None or s.cost < best.cost:
                best = s
        sol = best
    resid = float(np.max(np.abs(sol.fun))) * 1e6
    return sol.x, resid


def run_one(row, per, refit_only=False):
    lc_path = LC_DIR / (f"{row.target}__{row.visit}_{row.det}__"
                        f"{row.planet.replace(' ', '_')}.csv")
    if not lc_path.exists():
        return None
    lc = pd.read_csv(lc_path)
    t = lc.t_bmjd.to_numpy()
    f = lc.flux.to_numpy()

    cfg = lam.TransitConfig(per=per, e=row.e, w_deg=row.w_deg,
                            ramp=row.ramp_model,
                            t_ref=float(np.median(t)),
                            span=float(t.max() - t.min()),
                            t_min=float(t.min()))
    dur_d = row.t14_h / 24.0
    init = {"t0": row.t0_bmjd, "k": row.k, "aRs": row.aRs, "b": row.b,
            "q1": row.q1, "q2": row.q2}
    coeffs, ramp_resid = recover_ramp_coeffs(lc, cfg)
    init.update({f"c{i}": c for i, c in enumerate(coeffs)})
    if ramp_resid > 5.0:
        print(f"  WARN: ramp recovery residual {ramp_resid:.1f} ppm",
              flush=True)
    priors = {"aRs": (row.aRs, 0.1 * row.aRs), "b": (row.b, 0.1)}

    sym, ts, fs = lam.fit_transit(t, f, init, cfg, priors, dur_d,
                                  n_restarts=3)
    asym, ta, fa = lam.fit_asymmetric(ts, fs, sym, priors, dur_d)
    dd_refit = lam.delta_depth_ppm(asym)
    if abs(dd_refit - row.dd_ppm) > 0.5 * abs(row.dd_ppm) + 20.0:
        print(f"  WARN: refit dd {dd_refit:+.0f} ppm vs catalog "
              f"{row.dd_ppm:+.0f} ppm — LSQ refit drifted", flush=True)
    if refit_only:
        return dd_refit, np.nan, np.nan

    dd_mc, dd_mc_err, _ = lam.mcmc_asym(ta, fa, asym, priors, dur_d)
    return dd_refit, dd_mc, dd_mc_err


def main():
    args = sys.argv[1:]
    thr = 2.0
    if "--thr" in args:
        thr = float(args[args.index("--thr") + 1])
    planet = None
    if "--planet" in args:
        planet = args[args.index("--planet") + 1]
    refit_only = "--refit-only" in args
    take_all = "--all" in args

    cat_path = OUT / "limb_asym_catalog.csv"
    cat = pd.read_csv(cat_path)
    eph = pd.read_csv(ROOT / "data" / "processed" / "survey_ephemerides.csv")
    sel = (cat.partial_flag == 0) & cat.dd_ppm.notna()
    if not take_all:
        sel &= cat.dd_sigma.abs() > thr
    if planet:
        sel &= cat.planet == planet
    print(f"{int(sel.sum())} rows to escalate "
          f"(thr={'all' if take_all else thr}, planet={planet or 'any'}, "
          f"{'refit only' if refit_only else 'refit + emcee'})")

    for i in cat.index[sel]:
        row = cat.loc[i]
        per = float(eph[eph.pl_name == row.planet].pl_orbper.iloc[0])
        print(f"{row.visit} {row.det} {row.planet}: catalog dd "
              f"{row.dd_ppm:+.0f} +- {row.dd_err_ppm:.0f} ppm "
              f"({row.dd_sigma:+.1f} sigma) -> "
              f"{'refit' if refit_only else 'emcee'}...", flush=True)
        res = run_one(row, per, refit_only)
        if res is None:
            print("  no light curve — skipped", flush=True)
            continue
        dd_refit, dd_mc, dd_mc_err = res
        if refit_only:
            print(f"  refit: dd {dd_refit:+.0f} ppm "
                  f"(catalog {row.dd_ppm:+.0f})", flush=True)
            continue
        cat.loc[i, "dd_mcmc_ppm"] = dd_mc
        cat.loc[i, "dd_mcmc_err_ppm"] = dd_mc_err
        z = dd_mc / dd_mc_err if dd_mc_err > 0 else np.nan
        print(f"  refit dd {dd_refit:+.0f}; emcee: dd {dd_mc:+.0f} +- "
              f"{dd_mc_err:.0f} ppm ({z:+.1f} sigma)", flush=True)

    if not refit_only:
        cat.to_csv(cat_path, index=False)
        print(f"catalog updated -> {cat_path}")


if __name__ == "__main__":
    main()
