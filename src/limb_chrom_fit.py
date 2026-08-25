"""Chromatic limb-asymmetry study (gap #8), stage 2: per-bin asymmetric fits.

For every full-transit (visit, det, planet) of study-05 and every wavelength
bin of stage 1: divide-white common-mode correction, then the study-05
asymmetric transit fit with the geometry pinned to the white solution.

Divide-white: cm = white_flux x white_ramp / white_model_sym (the visit's
common-mode systematics INCLUDING the ramp and the shared noise); the bin
curve divided by cm keeps its own transit and loses the common systematics.
Consequence worth stating: any asymmetry present in the WHITE curve lives in
cm, so the per-bin dd measured here is the CHROMATIC CONTRAST relative to
the band-averaged (white) asymmetry -- which is the physically interesting
quantity, and the white study already showed the band average is null.

Per (row, bin): symmetric fit (ramp by BIC over the full
linear/quad/quad_exp menu; t0/aRs/b/q1/q2
pinned to the white catalog values by tight Gaussian priors, k free) ->
asymmetric refit (k_in/k_eg) -> dd = (k_eg^2 - k_in^2) x 1e6 [ppm] with a
cyclic-shift residual bootstrap error. Calibration of these errors by
injection-recovery is stage 3 -- as in study-05, the raw errors are NOT to
be trusted for significance claims.

Output: reports/limb_chrom/limb_chrom_catalog.csv
Usage: python limb_chrom_fit.py [target ...] [--visit jwN...] [--n-boot 200]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from survey_analyze import ROOT                                  # noqa: E402
import limb_asym_model as lam                                    # noqa: E402

import numpy as np    # noqa: E402
import pandas as pd   # noqa: E402

CAT = ROOT / "reports" / "limb_asymmetry" / "limb_asym_catalog.csv"
WHITE_LC = ROOT / "reports" / "limb_asymmetry" / "lightcurves"
CHROM = ROOT / "reports" / "limb_chrom"
N_BOOT = 200
PIN = {"t0": 1e-5, "aRs": 1e-3, "b": 1e-3, "q1": 1e-3, "q2": 1e-3}


def fit_bin(t, f, row, per, n_boot):
    """Symmetric+asymmetric fit of one common-mode-corrected bin curve."""
    dur_d = row.t14_h / 24.0
    init = {"t0": row.t0_bmjd, "k": row.k, "aRs": row.aRs, "b": row.b,
            "q1": row.q1, "q2": row.q2}
    priors = {n: (init[n], sd) for n, sd in PIN.items()}
    # full ramp menu incl. quad_exp: the 2026-08-22 sensitivity test showed
    # the settling exponential shifts dd by up to 1.4 sigma_cal in bins where
    # BIC prefers it (LHS 1140 b CH4 band) -- restricting to linear/quad
    # biases exactly the interesting bins (the study-05 degeneracy lesson)
    (sym, tt, ff), bics = lam.fit_with_ramp_selection(
        t, f, init, priors, dur_d, per, row.e, row.w_deg,
        ramps=("linear", "quad", "quad_exp"))
    asym, tt, ff = lam.fit_asymmetric(tt, ff, sym, priors, dur_d)
    dd = lam.delta_depth_ppm(asym)
    dd_err, _ = lam.bootstrap_delta_depth(tt, ff, asym, priors, dur_d,
                                          n_boot=n_boot)
    return {
        "dd_ppm": dd, "dd_err_boot_ppm": dd_err,
        "k_in": asym["k_in"], "k_eg": asym["k_eg"],
        "depth_ppm": sym["k"] ** 2 * 1e6,
        "ramp_model": sym.cfg.ramp, "sigma_ppm": asym.sigma * 1e6,
        "chi2_red": asym.chi2 / asym.dof, "n_pts": asym.dof + len(asym.names),
    }


def main() -> None:
    argv = sys.argv[1:]
    n_boot, only_visit = N_BOOT, None
    if "--n-boot" in argv:
        i = argv.index("--n-boot")
        n_boot = int(argv[i + 1])
        argv = argv[:i] + argv[i + 2:]
    if "--visit" in argv:
        i = argv.index("--visit")
        only_visit = argv[i + 1]
        argv = argv[:i] + argv[i + 2:]
    only = set(argv)

    cat = pd.read_csv(CAT)
    cat = cat[cat.partial_flag == 0]
    eph = pd.read_csv(ROOT / "data" / "processed" / "survey_ephemerides.csv")
    idx = pd.read_csv(CHROM / "bins_index.csv")

    out_rows = []
    out_path = CHROM / "limb_chrom_catalog.csv"
    done = set()
    if out_path.exists():
        prev = pd.read_csv(out_path)
        out_rows = prev.to_dict("records")
        done = {(r["visit"], r["det"], r["planet"], r["bin"])
                for r in out_rows}

    for _, row in cat.iterrows():
        if only and row.target not in only:
            continue
        if only_visit and row.visit != only_visit:
            continue
        bins_f = CHROM / "lightcurves" / \
            f"{row.target}__{row.visit}_{row.det}__bins.csv"
        white_f = WHITE_LC / (f"{row.target}__{row.visit}_{row.det}__"
                              f"{row.planet.replace(' ', '_')}.csv")
        if not bins_f.exists() or not white_f.exists():
            print(f"{row.visit} {row.det} {row.planet}: missing input "
                  f"-- skipped", flush=True)
            continue
        b = pd.read_csv(bins_f)
        w = pd.read_csv(white_f)
        m = pd.merge(b.assign(_t=b.t_bmjd.round(8)),
                     w.assign(_t=w.t_bmjd.round(8)), on="_t",
                     suffixes=("", "_w"))
        if len(m) < 0.9 * len(w):
            print(f"{row.visit} {row.det}: bin/white time-grid mismatch "
                  f"({len(m)}/{len(w)}) -- skipped", flush=True)
            continue
        cm = m.flux * m.ramp / m.model_sym
        per = float(eph[eph.pl_name == row.planet].pl_orbper.iloc[0])
        binfo = idx[(idx.visit == row.visit) & (idx.det == row.det)]

        for k in sorted(int(x) for x in binfo["bin"].unique()):
            key = (row.visit, row.det, row.planet, k)
            if key in done:
                continue
            f = (m[f"bin{k:02d}"] / cm).to_numpy(float)
            t = m.t_bmjd.to_numpy(float)
            bi = binfo[binfo["bin"] == k].iloc[0]
            try:
                res = fit_bin(t, f, row, per, n_boot)
            except Exception as exc:   # keep the sweep alive, log the cell
                print(f"{row.visit} {row.det} {row.planet} bin{k}: "
                      f"FIT FAILED ({exc})", flush=True)
                continue
            res.update({"target": row.target, "planet": row.planet,
                        "visit": row.visit, "det": row.det, "bin": k,
                        "wave_lo_um": bi.wave_lo_um,
                        "wave_hi_um": bi.wave_hi_um,
                        "spot_flag": row.spot_flag,
                        "dd_white_ppm": row.dd_ppm})
            out_rows.append(res)
            print(f"{row.visit} {row.det} {row.planet} bin{k}: "
                  f"dd={res['dd_ppm']:+.0f}+-{res['dd_err_boot_ppm']:.0f} ppm"
                  f" ({res['ramp_model']}, chi2r {res['chi2_red']:.2f})",
                  flush=True)
            pd.DataFrame(out_rows).to_csv(out_path, index=False)

    print(f"catalog: {len(out_rows)} rows -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
