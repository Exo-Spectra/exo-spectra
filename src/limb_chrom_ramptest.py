"""Chromatic limb-asymmetry study (gap #8) — ramp-model sensitivity test.

HISTORICAL TEST: the original stage 2 restricted the per-bin baseline to
linear/quad (divide-white should remove the visit ramp); this test showed
the restriction biased dd, so the released limb_chrom_fit.py now uses the
full menu — re-running this against it yields near-zero shifts by design. But the study-05 lesson says a settling exponential
is partially degenerate with the ingress depth — if a bin retains a
RESIDUAL exponential relative to white, it could leak into dd. For every
band bin of the FDR-passing planets: refit with the full ramp menu
(linear/quad/quad_exp, BIC selection) and compare dd.

-> reports/limb_chrom/ramp_sensitivity.csv (+ stdout table)
Usage: python limb_chrom_ramptest.py [planets...; default: FDR-passing]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import limb_asym_model as lam                                  # noqa: E402
from limb_asym_inject import row_per                           # noqa: E402
from limb_chrom_fit import PIN, CAT, CHROM                     # noqa: E402
from limb_chrom_inject import corrected_curve                  # noqa: E402


def refit_full_ramps(wrow, crow, per):
    t, f = corrected_curve(wrow, crow)
    dur_d = wrow.t14_h / 24.0
    init = {"t0": wrow.t0_bmjd, "k": wrow.k, "aRs": wrow.aRs, "b": wrow.b,
            "q1": wrow.q1, "q2": wrow.q2}
    priors = {n: (init[n], sd) for n, sd in PIN.items()}
    (sym, tt, ff), bics = lam.fit_with_ramp_selection(
        t, f, init, priors, dur_d, per, wrow.e, wrow.w_deg,
        ramps=("linear", "quad", "quad_exp"))
    asym, tt, ff = lam.fit_asymmetric(tt, ff, sym, priors, dur_d)
    return lam.delta_depth_ppm(asym), sym.cfg.ramp


def main():
    fdr = pd.read_csv(CHROM / "limb_chrom_fdr.csv")
    contrasts = pd.read_csv(CHROM / "band_contrast.csv")
    wcat = pd.read_csv(CAT)
    wcat = wcat[wcat.partial_flag == 0]

    planets = set(sys.argv[1:]) or set(contrasts[contrasts.fdr_pass].planet)
    sel = fdr[fdr.planet.isin(planets) & (fdr.band != "cont")]
    print(f"ramp test: {len(sel)} band bins of {sorted(planets)}",
          flush=True)

    rows = []
    for _, crow in sel.iterrows():
        wrow = wcat[(wcat.visit == crow.visit) & (wcat.det == crow.det)
                    & (wcat.planet == crow.planet)].iloc[0]
        try:
            dd_full, ramp = refit_full_ramps(wrow, crow, row_per(crow))
        except Exception as exc:
            print(f"{crow.visit} {crow.det} bin{int(crow.bin)}: "
                  f"FAILED ({exc})", flush=True)
            continue
        rows.append({"planet": crow.planet, "visit": crow.visit,
                     "det": crow.det, "bin": int(crow.bin),
                     "band": crow.band, "dd_stage2_ppm": crow.dd_ppm,
                     "ramp_stage2": crow.ramp_model,
                     "dd_fullramp_ppm": dd_full, "ramp_full": ramp,
                     "shift_sigma": (dd_full - crow.dd_ppm)
                     / crow.sigma_cal_ppm})
        print(f"{crow.planet} {crow.visit} {crow.det} bin{int(crow.bin)} "
              f"({crow.band}): dd {crow.dd_ppm:+.0f} ({crow.ramp_model}) -> "
              f"{dd_full:+.0f} ppm ({ramp}), przesuniecie "
              f"{rows[-1]['shift_sigma']:+.2f} sigma_cal", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(CHROM / "ramp_sensitivity.csv", index=False)
    if len(df):
        print(f"\nmedian |shift| = {df.shift_sigma.abs().median():.2f} "
              f"sigma_cal; max = {df.shift_sigma.abs().max():.2f}; "
              f"quad_exp wybrany w {int((df.ramp_full == 'quad_exp').sum())}"
              f"/{len(df)} binach")


if __name__ == "__main__":
    main()
