"""Chromatic limb-asymmetry study (gap #8) — limb-darkening sensitivity test.

The stage-2 fits pin q1/q2 to the WHITE-light values in every wavelength
bin, but limb darkening is chromatic: a wrong per-bin LD could in principle
masquerade as a wavelength-dependent dd. For every bin of the planets whose
band contrasts pass FDR, refit with q1/q2 FREE (only t0/aRs/b stay pinned)
and compare dd. If dd moves by less than ~1 sigma_cal, the pinned-LD
systematic cannot explain the band contrast.

-> reports/limb_chrom/ld_sensitivity.csv (+ stdout table)
Usage: python limb_chrom_ldtest.py [planet substrings...; default: the
       FDR-passing contrast planets]
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

PIN_NOLD = {n: sd for n, sd in PIN.items() if n not in ("q1", "q2")}


def refit(wrow, crow, per, pin):
    t, f = corrected_curve(wrow, crow)
    dur_d = wrow.t14_h / 24.0
    init = {"t0": wrow.t0_bmjd, "k": wrow.k, "aRs": wrow.aRs, "b": wrow.b,
            "q1": wrow.q1, "q2": wrow.q2}
    priors = {n: (init[n], sd) for n, sd in pin.items()}
    cfg = lam.TransitConfig(per=per, e=wrow.e, w_deg=wrow.w_deg,
                            ramp=crow.ramp_model,
                            t_ref=float(np.median(t[np.isfinite(f)])),
                            span=float(t.max() - t.min()),
                            t_min=float(t.min()))
    sym, tt, ff = lam.fit_transit(t, f, init, cfg, priors, dur_d,
                                  n_restarts=1)
    asym, tt, ff = lam.fit_asymmetric(tt, ff, sym, priors, dur_d)
    return lam.delta_depth_ppm(asym), asym["q1"], asym["q2"]


def main():
    fdr = pd.read_csv(CHROM / "limb_chrom_fdr.csv")
    contrasts = pd.read_csv(CHROM / "band_contrast.csv")
    wcat = pd.read_csv(CAT)
    wcat = wcat[wcat.partial_flag == 0]

    planets = set(sys.argv[1:]) or \
        set(contrasts[contrasts.fdr_pass].planet)
    sel = fdr[fdr.planet.isin(planets) & (fdr.band != "cont")]
    print(f"LD test: {len(sel)} band bins of {sorted(planets)}", flush=True)

    rows = []
    for _, crow in sel.iterrows():
        wrow = wcat[(wcat.visit == crow.visit) & (wcat.det == crow.det)
                    & (wcat.planet == crow.planet)].iloc[0]
        per = row_per(crow)
        try:
            dd_free, q1f, q2f = refit(wrow, crow, per, PIN_NOLD)
        except Exception as exc:
            print(f"{crow.visit} {crow.det} bin{int(crow.bin)}: "
                  f"FAILED ({exc})", flush=True)
            continue
        rows.append({"planet": crow.planet, "visit": crow.visit,
                     "det": crow.det, "bin": int(crow.bin),
                     "band": crow.band, "dd_pinned_ppm": crow.dd_ppm,
                     "dd_freeLD_ppm": dd_free,
                     "shift_sigma": (dd_free - crow.dd_ppm)
                     / crow.sigma_cal_ppm,
                     "q1_white": wrow.q1, "q1_free": q1f,
                     "q2_white": wrow.q2, "q2_free": q2f})
        print(f"{crow.planet} {crow.visit} {crow.det} bin{int(crow.bin)} "
              f"({crow.band}): dd {crow.dd_ppm:+.0f} -> {dd_free:+.0f} ppm "
              f"(przesuniecie {rows[-1]['shift_sigma']:+.2f} sigma_cal)",
              flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(CHROM / "ld_sensitivity.csv", index=False)
    if len(df):
        print(f"\nmedian |shift| = {df.shift_sigma.abs().median():.2f} "
              f"sigma_cal; max = {df.shift_sigma.abs().max():.2f}")


if __name__ == "__main__":
    main()
