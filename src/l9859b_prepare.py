"""Gap #1 verification, step 1: unify every L 98-59 b spectrum variant into
one format (wave_um, dwave_um, depth_ppm, err_ppm) for the retrieval matrix.

Variants written to data/processed/l9859b/:
  pub_firefly.csv   spec_id 921, Bello-Arufe 2025, FIREFLy reduction (80 pts)
  pub_eureka.csv    spec_id 922, Bello-Arufe 2025, Eureka! reduction (218 pts)
  own_v1..own_v4.csv  our uniform re-extraction per visit (survey stage 1)
  own_avg.csv       inverse-variance weighted average of the 4 own visits
  own_avg_infl.csv  same, errors inflated per-bin by sqrt(max(1, chi2_red))
                    -- carries the visit-to-visit inconsistency we measured
                    into the averaged error bars (published analyses do not)

Caveat (by design, documented in the study README, "Caveats"): our per-visit
errors are white out-of-transit scatter only, no correlated-noise term.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from spectra_io import load_spectrum  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed" / "l9859b"
SURVEY = ROOT / "reports" / "survey" / "spectra"

PUB = {  # spec_id -> (file, out name)
    921: ("80-53-84-86__L_98_59_b_3.12372_5578_1.tbl", "pub_firefly"),
    922: ("44-59-69-92__L_98_59_b_3.12372_5578_2.tbl", "pub_eureka"),
}
OWN_VISITS = ["jw03942001001", "jw03942002001", "jw03942003001", "jw03942004001"]


def load_pub(fname):
    df = load_spectrum(ROOT / "data" / "spectra" / fname, "Transmission")
    # value is PL_TRANDEP in %, symmetric errors after symmetrisation
    err = 0.5 * (df["err_hi"].abs() + df["err_lo"].abs())
    out = pd.DataFrame({
        "wave_um": df["wave"],
        "dwave_um": df["dwave"],
        "depth_ppm": df["value"] * 1e4,
        "err_ppm": err * 1e4,
    })
    out = out.dropna(subset=["wave_um", "depth_ppm", "err_ppm"]).reset_index(drop=True)
    # Eureka spectrum has null BANDWIDTH -> fill with the median grid spacing
    if out["dwave_um"].isna().any():
        step = float(np.median(np.diff(out["wave_um"])))
        out["dwave_um"] = out["dwave_um"].fillna(step)
    return out


def load_own(visit):
    df = pd.read_csv(SURVEY / f"L-98-59__{visit}__L_98-59_b.csv")
    return pd.DataFrame({
        "wave_um": df["wave"],
        "dwave_um": df["dwave"],
        "depth_ppm": df["depth_pct"] * 1e4,
        "err_ppm": df["err_pct"] * 1e4,
        "det": df["det"],
    })


def average_visits(visits):
    """Weighted average on visit 1's bin grid (grids differ by <1e-3 um)."""
    ref = visits[0]
    rows = []
    for i in range(len(ref)):
        w0, dw0, det0 = ref["wave_um"][i], ref["dwave_um"][i], ref["det"][i]
        d, e = [], []
        for v in visits:
            sel = v[(v["det"] == det0)
                    & (np.abs(v["wave_um"] - w0) < 0.5 * dw0)]
            if len(sel) == 1:
                d.append(float(sel["depth_ppm"].iloc[0]))
                e.append(float(sel["err_ppm"].iloc[0]))
        if len(d) < 2:
            continue
        d, e = np.array(d), np.array(e)
        w = 1.0 / e**2
        mean = np.sum(w * d) / np.sum(w)
        err = np.sqrt(1.0 / np.sum(w))
        chi2_red = np.sum(w * (d - mean) ** 2) / (len(d) - 1)
        rows.append((w0, dw0, mean, err, err * np.sqrt(max(1.0, chi2_red)),
                     len(d), chi2_red))
    return pd.DataFrame(rows, columns=[
        "wave_um", "dwave_um", "depth_ppm", "err_ppm", "err_infl_ppm",
        "n_visits", "chi2_red"])


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    for _, (fname, name) in PUB.items():
        df = load_pub(fname)
        df.to_csv(OUT / f"{name}.csv", index=False)
        print(f"{name}: {len(df)} pts, {df['wave_um'].min():.3f}-"
              f"{df['wave_um'].max():.3f} um, median err "
              f"{df['err_ppm'].median():.0f} ppm")

    visits = [load_own(v) for v in OWN_VISITS]
    for i, (v, df) in enumerate(zip(OWN_VISITS, visits), 1):
        df.drop(columns="det").to_csv(OUT / f"own_v{i}.csv", index=False)
        print(f"own_v{i} ({v}): {len(df)} pts, median err "
              f"{df['err_ppm'].median():.0f} ppm")

    avg = average_visits(visits)
    avg[["wave_um", "dwave_um", "depth_ppm", "err_ppm"]].to_csv(
        OUT / "own_avg.csv", index=False)
    infl = avg[["wave_um", "dwave_um", "depth_ppm", "err_infl_ppm"]].rename(
        columns={"err_infl_ppm": "err_ppm"})
    infl.to_csv(OUT / "own_avg_infl.csv", index=False)
    print(f"own_avg: {len(avg)} bins, median err {avg['err_ppm'].median():.0f}"
          f" ppm (inflated {avg['err_infl_ppm'].median():.0f} ppm), "
          f"median per-bin chi2_red {avg['chi2_red'].median():.2f}")

    # injection (negative) test: flat truth + white noise at own_avg errors;
    # the framework must NOT prefer an atmosphere on this by >2 sigma
    rng = np.random.default_rng(20260818)
    synth = avg[["wave_um", "dwave_um", "err_ppm"]].copy()
    synth["depth_ppm"] = (avg["depth_ppm"].median()
                          + rng.normal(0, avg["err_ppm"]))
    synth[["wave_um", "dwave_um", "depth_ppm", "err_ppm"]].to_csv(
        OUT / "synth_flat.csv", index=False)
    print(f"synth_flat: {len(synth)} bins (truth = flat "
          f"{avg['depth_ppm'].median():.0f} ppm)")


if __name__ == "__main__":
    main()
