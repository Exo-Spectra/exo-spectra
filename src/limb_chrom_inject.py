"""Chromatic limb-asymmetry study (gap #8), stage 3a: injection-recovery
calibration of the per-bin dd errors.

Mirrors limb_asym_inject.py, adapted to the divide-white per-bin pipeline:
for every calibratable bin of limb_chrom_catalog.csv (dd_err_boot below
ERR_SKIP -- the red-end garbage bins are excluded up front), rebuild the
common-mode-corrected curve, refit sym+asym once to recover the fitted
model, ramp and residuals (the catalog does not store ramp coefficients --
the study-05 lesson), then inject known asymmetries at the bin's own k,
re-add cyclically shifted residuals, refit, and measure what comes back.

Products per bin: sigma_cal (scatter of recovered dd at zero injection),
pooled sigma_resid, bias slope, UL95 = |dd_obs| + 1.645 sigma_cal.
-> reports/limb_chrom/limb_chrom_calibration.csv (append/resume safe)

Usage: python limb_chrom_inject.py [targets...]
(amplitudes/realizations are the N_AMP/N_REAL/N_REAL_ZERO constants below)
Run after: limb_chrom_fit.py
"""
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import numpy as np    # noqa: E402
import pandas as pd   # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import limb_asym_model as lam                                  # noqa: E402
from survey_analyze import ROOT                                # noqa: E402
from limb_asym_inject import split_k, row_per                  # noqa: E402
from limb_chrom_fit import PIN, CAT, WHITE_LC, CHROM           # noqa: E402

ERR_SKIP = 500.0    # ppm; bins with a worse bootstrap error are hopeless
N_AMP = 4           # non-zero amplitudes
N_REAL = 7          # realizations per amplitude
N_REAL_ZERO = 12    # zero-amplitude realizations (sets sigma_cal)


def corrected_curve(wrow, crow):
    """(t, f) of the divide-white corrected bin curve (as in stage 2)."""
    b = pd.read_csv(CHROM / "lightcurves" /
                    f"{crow.target}__{crow.visit}_{crow.det}__bins.csv")
    w = pd.read_csv(WHITE_LC / (f"{crow.target}__{crow.visit}_{crow.det}__"
                                f"{crow.planet.replace(' ', '_')}.csv"))
    m = pd.merge(b.assign(_t=b.t_bmjd.round(8)),
                 w.assign(_t=w.t_bmjd.round(8)), on="_t", suffixes=("", "_w"))
    cm = m.flux * m.ramp / m.model_sym
    f = (m[f"bin{int(crow.bin):02d}"] / cm).to_numpy(float)
    return m.t_bmjd.to_numpy(float), f


def run_one(wrow, crow, per, rng):
    t, f = corrected_curve(wrow, crow)
    dur_d = wrow.t14_h / 24.0
    init = {"t0": wrow.t0_bmjd, "k": wrow.k, "aRs": wrow.aRs, "b": wrow.b,
            "q1": wrow.q1, "q2": wrow.q2}
    priors = {n: (init[n], sd) for n, sd in PIN.items()}
    cfg = lam.TransitConfig(per=per, e=wrow.e, w_deg=wrow.w_deg,
                            ramp=crow.ramp_model,
                            t_ref=float(np.median(t[np.isfinite(f)])),
                            span=float(t.max() - t.min()),
                            t_min=float(t.min()))
    sym, tt, ff = lam.fit_transit(t, f, init, cfg, priors, dur_d,
                                  n_restarts=1)
    asym, tt, ff = lam.fit_asymmetric(tt, ff, sym, priors, dur_d)
    res = asym.residuals
    # noiseless ramp of the fitted asym model: model / transit part
    k_bin = float(np.sqrt(0.5 * (asym["k_in"] ** 2 + asym["k_eg"] ** 2)))
    inc = lam.impact_to_inc(asym["b"], asym["aRs"], cfg.e, cfg.w_deg)
    u1, u2 = lam.ld_coeffs(asym["q1"], asym["q2"])
    z, front = lam.proj_sep(tt, asym["t0"], per, asym["aRs"], inc,
                            cfg.e, cfg.w_deg)
    tr_fit = np.ones_like(tt)
    for kk, side in ((asym["k_in"], tt <= asym["t0"]),
                     (asym["k_eg"], tt > asym["t0"])):
        m = front & side & (z < 1.0 + kk)
        tr_fit[m] = 1.0 - lam.occult(z[m], kk, u1, u2)
    ramp_fit = asym.model / tr_fit

    step = max(crow.dd_err_boot_ppm, 10.0)
    amps = np.concatenate([[0.0], np.linspace(-3, 3, N_AMP) * step])
    amps = np.unique(np.round(amps, 3))
    init_a = {n: asym[n] for n in asym.names}

    recs = []
    for amp in amps:
        k_in, k_eg = split_k(k_bin, amp)
        if not np.isfinite(k_in):
            continue
        tr = np.ones_like(tt)
        for kk, side in ((k_in, tt <= asym["t0"]), (k_eg, tt > asym["t0"])):
            m = front & side & (z < 1.0 + kk)
            tr[m] = 1.0 - lam.occult(z[m], kk, u1, u2)
        base = ramp_fit * tr
        n_real = N_REAL_ZERO if amp == 0.0 else N_REAL
        for r in range(n_real):
            shift = rng.integers(1, len(tt) - 1)
            f_inj = base + np.roll(res, shift)
            try:
                fit, _, _ = lam.fit_transit(tt, f_inj, init_a, asym.cfg,
                                            priors, dur_d, n_restarts=1)
                recs.append({"amp": amp, "dd_hat": lam.delta_depth_ppm(fit)})
            except Exception:
                continue
    if len(recs) < 10:
        return None
    rec = pd.DataFrame(recs)
    at0 = rec[rec.amp == 0.0].dd_hat
    sigma_cal = float(at0.std(ddof=1)) if len(at0) > 3 else np.nan
    coef = np.polyfit(rec.amp, rec.dd_hat, 1)
    sigma_resid = float((rec.dd_hat - np.polyval(coef, rec.amp)).std(ddof=2))
    return {"target": crow.target, "planet": crow.planet,
            "visit": crow.visit, "det": crow.det, "bin": int(crow.bin),
            "wave_lo_um": crow.wave_lo_um, "wave_hi_um": crow.wave_hi_um,
            "n_fits": len(rec), "sigma_cal_ppm": sigma_cal,
            "sigma_resid_ppm": sigma_resid, "bias_slope": float(coef[0]),
            "dd_obs_ppm": crow.dd_ppm,
            "cal_vs_boot": sigma_cal / crow.dd_err_boot_ppm,
            "ul95_abs_dd_ppm": abs(crow.dd_ppm) + 1.645 * sigma_cal}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    only = set(args)
    wcat = pd.read_csv(CAT)
    wcat = wcat[wcat.partial_flag == 0]
    ccat = pd.read_csv(CHROM / "limb_chrom_catalog.csv")
    out_path = CHROM / "limb_chrom_calibration.csv"
    out_rows, done = [], set()
    if out_path.exists():
        prev = pd.read_csv(out_path)
        out_rows = prev.to_dict("records")
        done = {(r["visit"], r["det"], r["planet"], r["bin"])
                for r in out_rows}

    todo = ccat[ccat.dd_err_boot_ppm < ERR_SKIP]
    print(f"{len(todo)} bins to calibrate "
          f"(skipping {len(ccat) - len(todo)} with err >= {ERR_SKIP} ppm)",
          flush=True)
    for _, crow in todo.iterrows():
        if only and crow.target not in only:
            continue
        key = (crow.visit, crow.det, crow.planet, int(crow.bin))
        if key in done:
            continue
        wsel = wcat[(wcat.visit == crow.visit) & (wcat.det == crow.det)
                    & (wcat.planet == crow.planet)]
        if not len(wsel):
            continue
        wrow = wsel.iloc[0]
        per = row_per(crow)
        rng = np.random.default_rng(11 + hash(crow.visit + crow.det
                                              + str(crow.bin)) % 10000)
        try:
            res = run_one(wrow, crow, per, rng)
        except Exception as exc:
            print(f"{crow.visit} {crow.det} bin{int(crow.bin)}: "
                  f"CAL FAILED ({exc})", flush=True)
            continue
        if res is None:
            continue
        out_rows.append(res)
        pd.DataFrame(out_rows).to_csv(out_path, index=False)
        print(f"{crow.visit} {crow.det} {crow.planet} bin{int(crow.bin)}: "
              f"sigma_cal {res['sigma_cal_ppm']:.0f} ppm "
              f"(boot x{res['cal_vs_boot']:.2f}, slope "
              f"{res['bias_slope']:.3f})", flush=True)
    print(f"calibration: {len(out_rows)} rows -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
