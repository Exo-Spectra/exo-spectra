"""Limb-asymmetry study — injection-recovery sensitivity calibration.

For every gated visit-detector from limb_asym_run.py: take the saved
light curve, replace the fitted transit with an asymmetric one at a known
Delta-depth, re-add the (cyclically shifted) real residuals — preserving
the visit's red noise — and refit with the cheap pipeline. Products per
visit-detector:
  * sigma_cal — empirical scatter of recovered dd at zero injection
    (THE calibrated error; decides whether the bootstrap or the
    beta-inflated error from the catalog is the honest one)
  * bias slope of recovered vs injected dd
  * 95% upper limit on |Delta depth|
-> reports/limb_asymmetry/upper_limits.csv + one PNG per visit-detector.

Usage: python limb_asym_inject.py [targets...] [--n-amp N] [--n-real N]
Run after: limb_asym_run.py
"""
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import limb_asym_model as lam                                # noqa: E402
from survey_analyze import ROOT                              # noqa: E402

OUT = ROOT / "reports" / "limb_asymmetry"
LC_DIR = OUT / "lightcurves"


def split_k(k: float, dd_ppm: float):
    """Mean-depth-preserving ingress/egress radii for a given Delta depth."""
    d = k * k
    half = 0.5 * dd_ppm * 1e-6
    if d - half <= 0 or d + half <= 0:
        return np.nan, np.nan
    return np.sqrt(d - half), np.sqrt(d + half)


def run_one(row, n_amp: int, n_real: int, seed: int = 3):
    lc_path = LC_DIR / (f"{row.target}__{row.visit}_{row.det}__"
                        f"{row.planet.replace(' ', '_')}.csv")
    if not lc_path.exists():
        return None
    lc = pd.read_csv(lc_path)
    t = lc.t_bmjd.to_numpy()
    ramp = lc.ramp.to_numpy()
    res_asym = (lc.flux - lc.model_asym).to_numpy()

    cfg = lam.TransitConfig(per=row_per(row), e=row.e, w_deg=row.w_deg,
                            ramp=row.ramp_model, asym=True,
                            t_ref=float(np.median(t)),
                            span=float(t.max() - t.min()),
                            t_min=float(t.min()))
    dur_d = row.t14_h / 24.0
    init = {"t0": row.t0_bmjd, "k_in": row.k, "k_eg": row.k,
            "aRs": row.aRs, "b": row.b, "q1": row.q1, "q2": row.q2}
    priors = {"aRs": (row.aRs, 0.1 * row.aRs), "b": (row.b, 0.1)}
    inc = lam.impact_to_inc(row.b, row.aRs, row.e, row.w_deg)
    u1, u2 = lam.ld_coeffs(row.q1, row.q2)
    z, front = lam.proj_sep(t, row.t0_bmjd, row_per(row), row.aRs, inc,
                            row.e, row.w_deg)

    # injection amplitudes scaled to the bootstrap error, up to +-4x
    step = max(row.dd_err_boot_ppm, 10.0)
    amps = np.concatenate([[0.0], np.linspace(-4, 4, n_amp) * step])
    amps = np.unique(np.round(amps, 3))

    rng = np.random.default_rng(seed + hash(row.visit + row.det) % 10000)
    recs = []
    for amp in amps:
        k_in, k_eg = split_k(row.k, amp)
        if not np.isfinite(k_in):
            continue
        tr = np.ones_like(t)
        for kk, side in ((k_in, t <= row.t0_bmjd), (k_eg, t > row.t0_bmjd)):
            m = front & side & (z < 1.0 + kk)
            tr[m] = 1.0 - lam.occult(z[m], kk, u1, u2)
        base = ramp * tr
        for r in range(n_real):
            shift = 0 if (amp == 0.0 and r == 0) else \
                rng.integers(1, len(t) - 1)
            f_inj = base + np.roll(res_asym, shift)
            try:
                fit, _, _ = lam.fit_transit(t, f_inj, init, cfg, priors,
                                            dur_d, n_restarts=1)
                recs.append({"amp": amp,
                             "dd_hat": lam.delta_depth_ppm(fit)})
            except Exception:
                continue
    if len(recs) < 10:
        return None
    rec = pd.DataFrame(recs)

    at0 = rec[rec.amp == 0.0].dd_hat
    sigma_cal = float(at0.std(ddof=1)) if len(at0) > 3 else np.nan
    coef = np.polyfit(rec.amp, rec.dd_hat, 1)
    # pooled variant: scatter around the bias line over ALL amplitudes —
    # same quantity as sigma_cal under a linear response, but ~n_amp x
    # more samples (the zero-amp std alone has ~24% error at n_real=10)
    sigma_resid = float((rec.dd_hat - np.polyval(coef, rec.amp)).std(ddof=2))
    ul95 = abs(row.dd_ppm) + 1.645 * sigma_cal
    out = {"target": row.target, "planet": row.planet, "visit": row.visit,
           "det": row.det, "n_fits": len(rec), "sigma_cal_ppm": sigma_cal,
           "sigma_cal_resid_ppm": sigma_resid,
           "bias_slope": float(coef[0]), "bias_intercept_ppm": float(coef[1]),
           "dd_obs_ppm": row.dd_ppm,
           "cal_vs_boot": sigma_cal / row.dd_err_boot_ppm,
           "cal_vs_beta": sigma_cal / row.dd_err_ppm,
           "ul95_abs_dd_ppm": ul95}

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(rec.amp, rec.dd_hat, ".", ms=5, alpha=0.6, color="C0")
    lim = rec.amp.abs().max() * 1.1
    ax.plot([-lim, lim], [-lim, lim], "-", color="0.7", lw=0.8,
            label="1:1")
    ax.plot([-lim, lim], np.polyval(coef, [-lim, lim]), "--", color="C3",
            lw=1.0, label=f"slope {coef[0]:.3f}")
    ax.set_xlabel("injected dd [ppm]")
    ax.set_ylabel("recovered dd [ppm]")
    ax.set_title(f"{row.planet} {row.visit} {row.det} — "
                 f"sigma_cal {sigma_cal:.0f} ppm")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / (f"inject__{row.target}__{row.visit}_{row.det}__"
                       f"{row.planet.replace(' ', '_')}.png"), dpi=150)
    plt.close(fig)
    return out


_EPH = None


def row_per(row):
    """Orbital period from the ephemerides CSV (not in the catalog row)."""
    global _EPH
    if _EPH is None:
        _EPH = pd.read_csv(ROOT / "data" / "processed" /
                           "survey_ephemerides.csv")
    return float(_EPH[_EPH.pl_name == row.planet].pl_orbper.iloc[0])


def main():
    args = sys.argv[1:]
    n_amp, n_real = 8, 10
    if "--n-amp" in args:
        n_amp = int(args[args.index("--n-amp") + 1])
        del args[args.index("--n-amp"):args.index("--n-amp") + 2]
    if "--n-real" in args:
        n_real = int(args[args.index("--n-real") + 1])
        del args[args.index("--n-real"):args.index("--n-real") + 2]
    only = {a for a in args if not a.startswith("--")}

    cat = pd.read_csv(OUT / "limb_asym_catalog.csv")
    cat = cat[(cat.partial_flag == 0) & cat.dd_ppm.notna()]
    if only:
        cat = cat[cat.target.isin(only)]
    print(f"{len(cat)} visit-detector rows to calibrate "
          f"({n_amp}+1 amps x {n_real} realizations)")

    outs = []
    for _, row in cat.iterrows():
        res = run_one(row, n_amp, n_real)
        if res is None:
            print(f"{row.visit} {row.det} {row.planet}: no light curve / "
                  f"too few fits — skipped", flush=True)
            continue
        outs.append(res)
        print(f"{row.visit} {row.det} {row.planet}: sigma_cal "
              f"{res['sigma_cal_ppm']:.0f} ppm (boot x{res['cal_vs_boot']:.2f}, "
              f"beta-infl x{res['cal_vs_beta']:.2f}), slope "
              f"{res['bias_slope']:.3f}, UL95 {res['ul95_abs_dd_ppm']:.0f} ppm",
              flush=True)

    if outs:
        df = pd.DataFrame(outs)
        df.to_csv(OUT / "upper_limits.csv", index=False)
        print(f"\n{len(df)} rows -> {OUT / 'upper_limits.csv'}")
        print(f"median calibration: sigma_cal/boot = "
              f"{df.cal_vs_boot.median():.2f}, sigma_cal/beta-inflated = "
              f"{df.cal_vs_beta.median():.2f}")


if __name__ == "__main__":
    main()
