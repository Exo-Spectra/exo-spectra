"""Study 8 v1.1 release export: absolute evidences + exact retrieval inputs.

Why: the v1.0 release published only lnB vs flat for 108 non-flat cells.
A machine-learning group that wants to benchmark against this matrix needs
the absolute lnZ of all 144 runs and the spectra exactly as fitted (cut to
the pair's common wavelength range, ppm). This script writes them from the
same build_cells() the retrievals used, so nothing is re-derived by hand.

Hard number gate (exit code != 0 blocks the release): the flat model's lnZ
is a 1-D integral, so it is recomputed here by quadrature from the exported
spectra and compared with every dynesty flat lnZ (36 cells: 18 spectrum x
error-variant baselines x 2 codes). Agreement within 4 x lnZ_err proves the
exported spectra/errors are the fitted inputs.

Runs only inside the working repo (reads reports/retrieval_robustness/ and
data/processed/); a public clone can re-check the flat lnZ from the exported
CSVs with flat_lnz() below.

Usage: python src/gap7_export_release.py --out <public>/studies/08-retrieval-robustness/results
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gap7_retrieval as g  # noqa: E402

COLS = ["spectrum", "pair", "planet", "axis", "side", "code", "model",
        "err_variant", "n_points", "lnZ", "lnZ_err", "max_lnL", "ncall",
        "nlive_used", "capped", "best_params", "post_mean", "post_std"]


def flat_lnz(depth, err):
    """Exact lnZ of the flat model: uniform prior [min-5max(err), max+5max(err)]."""
    lo, hi = depth.min() - 5 * err.max(), depth.max() + 5 * err.max()
    # lnL(c) is a Gaussian in c: -0.5*S0*(c-mu)^2 + const -> closed-form integral
    w = 1.0 / err ** 2
    s0, s1, s2 = w.sum(), (w * depth).sum(), (w * depth ** 2).sum()
    mu = s1 / s0
    const = -0.5 * (s2 - s1 ** 2 / s0) - 0.5 * np.log(2 * np.pi * err ** 2).sum()
    k = np.sqrt(s0)
    mass = ndtr((hi - mu) * k) - ndtr((lo - mu) * k)      # truncation by the prior
    return const + 0.5 * np.log(2 * np.pi / s0) + np.log(mass) - np.log(hi - lo)


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--out", required=True)
    out = Path(a.parse_args().out)
    (out / "input_spectra").mkdir(parents=True, exist_ok=True)

    ev = pd.concat([pd.read_csv(g.OUT / f"evidences_{c}.csv") for c in ("taurex", "platon")],
                   ignore_index=True)
    if len(ev) != 144 or ev.duplicated(["spectrum", "code", "model", "err_variant"]).any():
        sys.exit(f"GATE FAIL: expected 144 unique runs, got {len(ev)}")

    cells = g.build_cells()
    setup, fails = [], []
    for cell in cells:
        sp = pd.DataFrame(dict(wave_um=cell["wave"], dwave_um=cell["dwave"],
                               depth_ppm=cell["depth"], err_ppm=cell["err"]))
        sp.to_csv(out / "input_spectra" / f"{cell['spectrum']}.csv", index=False,
                  float_format="%.10g")
        for variant in cell["variants"]:
            err = cell["err"] * (np.sqrt(g.INFLATE[cell["pair"]]) if variant == "infl" else 1.0)
            rows = ev[(ev.spectrum == cell["spectrum"]) & (ev.err_variant == variant)]
            if (rows.n_points != len(sp)).any():
                fails.append(f"{cell['spectrum']}/{variant}: n_points {rows.n_points.unique()} vs {len(sp)}")
            for _, r in rows[rows.model == "flat"].iterrows():
                z = flat_lnz(cell["depth"], err)
                d = z - r.lnZ
                ok = abs(d) <= 4 * r.lnZ_err
                print(f"flat {cell['spectrum']:<20} {variant:<4} {r.code:<6} "
                      f"dynesty {r.lnZ:12.4f} +/- {r.lnZ_err:.3f}  quad {z:12.4f}  "
                      f"diff {d:+.3f}  {'OK' if ok else 'FAIL'}")
                if not ok:
                    fails.append(f"flat lnZ {cell['spectrum']}/{variant}/{r.code}: diff {d:+.3f}")
            t_lo = max(100.0, g.T_FACT[0] * cell["teq"]); t_hi = min(3400.0, g.T_FACT[1] * cell["teq"])
            setup.append(dict(
                spectrum=cell["spectrum"], pair=cell["pair"], planet=cell["planet"],
                side=cell["side"], err_variant=variant,
                err_scale=float(np.sqrt(g.INFLATE[cell["pair"]])) if variant == "infl" else 1.0,
                source_file=Path(cell["file"]).name, n_points=len(sp),
                wave_min_um=float(cell["wave"].min()), wave_max_um=float(cell["wave"].max()),
                model_lam_lo_um=cell["lam_lo"], model_lam_hi_um=cell["lam_hi"],
                mp_me=cell["mp_me"], rp_re=cell["rp_re"], rs_sun=cell["rs_sun"],
                teff_k=cell["teff"], teq_k=cell["teq"],
                models=",".join(cell["models"]),
                prior_flat=json.dumps([round(float(cell["depth"].min() - 5 * err.max()), 6),
                                       round(float(cell["depth"].max() + 5 * err.max()), 6)]),
                prior_rp_re=json.dumps([g.PRIOR_RP_FACT[0] * cell["rp_re"], g.PRIOR_RP_FACT[1] * cell["rp_re"]]),
                prior_T_k=json.dumps([t_lo, t_hi]),
                prior_logx=json.dumps(list(g.PRIOR_LOGX)),
                prior_logpc_pa=json.dumps(list(g.PRIOR_LOGPC)),
            ))
    keys = {c["spectrum"] for c in cells}
    if set(ev.spectrum) != keys:
        fails.append(f"spectrum keys differ: {sorted(set(ev.spectrum) ^ keys)}")

    # 144 cells come from fewer distinct nested runs: the flat baseline does not
    # depend on the RT code (one run serves both codes in most cells), and
    # TOI-5205 b side a is the same spectrum in two pairs. Those were computed
    # separately but with identical input and seed 42, so they are identical
    # runs (same lnZ and ncall); run_id makes that explicit for benchmark users.
    ev = ev[COLS].sort_values(["pair", "side", "err_variant", "code", "model"])
    for c in ("n_points", "ncall", "nlive_used", "capped"):
        ev[c] = ev[c].astype("Int64")          # integers, empty = not recorded
    sig = ev.lnZ.round(8).astype(str) + "|" + ev.ncall.astype(str)
    ev.insert(0, "run_id", pd.factorize(sig)[0] + 1)
    print(f"distinct nested-sampling runs: {ev.run_id.nunique()} (rows: {len(ev)}; "
          f"capped at the call budget: {int((ev.drop_duplicates('run_id').capped == 1).sum())})")
    ev.to_csv(out / "evidences_all.csv", index=False)
    pd.DataFrame(setup).to_csv(out / "retrieval_setup.csv", index=False)
    print(f"\nwritten: evidences_all.csv ({len(ev)} rows), retrieval_setup.csv "
          f"({len(setup)} rows), input_spectra/ ({len(cells)} files)")
    if fails:
        print("GATE FAIL:\n  " + "\n  ".join(fails)); sys.exit(1)
    print("GATE OK")


if __name__ == "__main__":
    main()
