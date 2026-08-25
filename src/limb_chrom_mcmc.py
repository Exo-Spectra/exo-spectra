"""Chromatic limb-asymmetry study (gap #8), stage 3c: MCMC escalation.

Study-05 discipline: every per-bin fit with |z_cal| >= Z_MIN gets a full
emcee posterior of the asymmetric model (free white-noise jitter). The
corrected curve and the sym+asym LSQ solution are rebuilt on the fly with
the SAME full ramp menu as stage 2 (the catalog does not store ramp
coefficients; rebuilding avoids the study-05 refit-drift gotcha because the
posterior starts from the fresh MAP).

-> reports/limb_chrom/limb_chrom_mcmc.csv (append/resume safe)
Usage: python limb_chrom_mcmc.py [--z-min 2.0]
Run after: limb_chrom_report.py (needs limb_chrom_fdr.csv)
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
from limb_asym_inject import row_per                           # noqa: E402
from limb_chrom_fit import PIN, CAT, CHROM                     # noqa: E402
from limb_chrom_inject import corrected_curve                  # noqa: E402

Z_MIN = 2.0


def main():
    z_min = Z_MIN
    if "--z-min" in sys.argv:
        z_min = float(sys.argv[sys.argv.index("--z-min") + 1])
    fdr = pd.read_csv(CHROM / "limb_chrom_fdr.csv")
    wcat = pd.read_csv(CAT)
    wcat = wcat[wcat.partial_flag == 0]

    out_path = CHROM / "limb_chrom_mcmc.csv"
    out_rows, done = [], set()
    if out_path.exists():
        prev = pd.read_csv(out_path)
        out_rows = prev.to_dict("records")
        done = {(r["visit"], r["det"], r["planet"], r["bin"])
                for r in out_rows}

    sel = fdr[np.abs(fdr.z_cal) >= z_min]
    print(f"{len(sel)} bins with |z_cal| >= {z_min}", flush=True)
    for _, crow in sel.iterrows():
        key = (crow.visit, crow.det, crow.planet, int(crow.bin))
        if key in done:
            continue
        wrow = wcat[(wcat.visit == crow.visit) & (wcat.det == crow.det)
                    & (wcat.planet == crow.planet)].iloc[0]
        per = row_per(crow)
        t, f = corrected_curve(wrow, crow)
        dur_d = wrow.t14_h / 24.0
        init = {"t0": wrow.t0_bmjd, "k": wrow.k, "aRs": wrow.aRs,
                "b": wrow.b, "q1": wrow.q1, "q2": wrow.q2}
        priors = {n: (init[n], sd) for n, sd in PIN.items()}
        try:
            (sym, tt, ff), _ = lam.fit_with_ramp_selection(
                t, f, init, priors, dur_d, per, wrow.e, wrow.w_deg,
                ramps=("linear", "quad", "quad_exp"))
            asym, tt, ff = lam.fit_asymmetric(tt, ff, sym, priors, dur_d)
            dd_map = lam.delta_depth_ppm(asym)
            dd_mc, dd_mc_err, _ = lam.mcmc_asym(tt, ff, asym, priors, dur_d)
        except Exception as exc:
            print(f"{crow.visit} {crow.det} bin{int(crow.bin)}: "
                  f"MCMC FAILED ({exc})", flush=True)
            continue
        out_rows.append({
            "target": crow.target, "planet": crow.planet,
            "visit": crow.visit, "det": crow.det, "bin": int(crow.bin),
            "band": crow.band, "ramp": asym.cfg.ramp,
            "dd_lsq_ppm": crow.dd_ppm, "sigma_cal_ppm": crow.sigma_cal_ppm,
            "z_cal": crow.z_cal, "dd_map_ppm": dd_map,
            "dd_mcmc_ppm": dd_mc, "dd_mcmc_err_ppm": dd_mc_err,
            "z_mcmc": dd_mc / dd_mc_err})
        pd.DataFrame(out_rows).to_csv(out_path, index=False)
        print(f"{crow.planet} {crow.visit} {crow.det} bin{int(crow.bin)} "
              f"({crow.band}): z_cal {crow.z_cal:+.2f} -> mcmc "
              f"{dd_mc:+.0f}+-{dd_mc_err:.0f} ppm (z {dd_mc/dd_mc_err:+.2f})",
              flush=True)
    print(f"mcmc: {len(out_rows)} rows -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
