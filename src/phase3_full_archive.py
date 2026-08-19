"""Phase 3: run the pair statistic over the whole archive.

Steps:
  1. compute chi2/p for every candidate pair (data/processed/pairs.csv)
  2. classify pairs:
       epoch_class (Eclipse only, from per-point OBS_DATE ranges):
         independent  date ranges disjoint -> genuinely different epochs
         same_obs     date ranges overlap  -> same observation(s)
         unknown      dates missing
       shared_suspect: statistically "too consistent" (p > 0.999 with n >= 8)
         -> likely re-reductions of the same underlying data
  3. Benjamini-Hochberg FDR over the clean-test subset
  4. write full results + summary + overlay plots for top anomalies

Outputs: data/processed/pair_results.csv, reports/phase3_summary.md,
         reports/anomaly_<rank>_<planet>.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from compare import compare_pair
from spectra_io import load_spectrum, usable_points

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed"
REPORTS = ROOT / "reports"

FDR_LEVEL = 0.01
SHARED_P = 0.999   # p above this with n >= SHARED_N -> suspected shared data
SHARED_N = 8
TOP_PLOTS = 10


def bh_fdr(p: np.ndarray, alpha: float) -> np.ndarray:
    """Benjamini-Hochberg: boolean mask of rejected (significant) tests."""
    n = len(p)
    order = np.argsort(p)
    thresh = alpha * (np.arange(1, n + 1) / n)
    passed = p[order] <= thresh
    k = np.max(np.nonzero(passed)[0]) + 1 if passed.any() else 0
    out = np.zeros(n, bool)
    out[order[:k]] = True
    return out


def epoch_class(row, summary) -> str:
    """Eclipse pairs: compare per-spectrum OBS_DATE ranges."""
    a = summary.loc[row.spec_id_a]
    b = summary.loc[row.spec_id_b]
    if pd.isna(a.obs_date_min) or pd.isna(b.obs_date_min):
        return "unknown"
    # disjoint date ranges (0.5 d tolerance) -> different epochs
    if a.obs_date_max < b.obs_date_min - 0.5 or b.obs_date_max < a.obs_date_min - 0.5:
        return "independent"
    return "same_obs"


def main() -> None:
    summary = pd.read_csv(OUT / "spectra_summary.csv").set_index("spec_id")
    pairs = pd.read_csv(OUT / "pairs.csv")
    REPORTS.mkdir(exist_ok=True)

    # load every needed spectrum once
    needed = sorted(set(pairs.spec_id_a) | set(pairs.spec_id_b))
    specs = {}
    for sid in needed:
        r = summary.loc[sid]
        specs[sid] = usable_points(load_spectrum(ROOT / "data" / "spectra" / r.file, r.spec_type))
    print(f"loaded {len(specs)} spectra")

    rows = []
    for p in pairs.itertuples():
        res = compare_pair(specs[p.spec_id_a], specs[p.spec_id_b])
        row = {
            "pl_name": p.pl_name, "spec_type": p.spec_type,
            "spec_id_a": p.spec_id_a, "spec_id_b": p.spec_id_b,
            "authors_a": p.authors_a, "authors_b": p.authors_b,
            "instrument_a": p.instrument_a, "instrument_b": p.instrument_b,
            "same_instrument": p.same_instrument,
            "epoch_class": epoch_class(p, summary) if p.spec_type == "Eclipse" else "unknown",
            "tested": res is not None,
        }
        if res is not None:
            row.update(n=res.n, n_interp=res.n_interp, offset=res.offset,
                       offset_err=res.offset_err, chi2=res.chi2, dof=res.dof,
                       p_value=res.p_value, chi2_red=res.chi2 / res.dof,
                       slope=res.slope, slope_err=res.slope_err,
                       p_slope=res.p_slope,
                       chi2_red_slope=(res.chi2_slope / res.dof_slope
                                       if res.dof_slope else None))
        rows.append(row)
    df = pd.DataFrame(rows)

    tested = df[df.tested].copy()
    tested["shared_suspect"] = (tested.p_value > SHARED_P) & (tested.n >= SHARED_N)
    # same_obs eclipse pairs are also shared data by construction
    tested.loc[tested.epoch_class == "same_obs", "shared_suspect"] = True

    # clean weather test: same instrument, not suspected shared data,
    # exclude Direct Imaging (FLAM units differ between papers; offset model inadequate)
    clean = tested[
        tested.same_instrument & ~tested.shared_suspect & (tested.spec_type != "Direct Imaging")
    ].copy()
    clean["significant"] = bh_fdr(clean.p_value.to_numpy(), FDR_LEVEL)

    # cross-instrument view (weaker: offset may not absorb all systematics)
    cross = tested[
        ~tested.same_instrument & ~tested.shared_suspect & (tested.spec_type != "Direct Imaging")
    ].copy()
    cross["significant"] = bh_fdr(cross.p_value.to_numpy(), FDR_LEVEL)

    full = pd.concat([clean, cross, tested[tested.shared_suspect | (tested.spec_type == "Direct Imaging")]])
    full = full.sort_values("p_value").reset_index(drop=True)
    full.to_csv(OUT / "pair_results.csv", index=False)

    sig_clean = clean[clean.significant].sort_values("p_value")
    sig_cross = cross[cross.significant].sort_values("p_value")

    lines = [
        "# Phase 3 — full-archive pair comparison\n",
        f"- candidate pairs: {len(df)}; tested (>=3 matched points): {len(tested)}",
        f"- flagged as shared data / re-reductions: {int(tested.shared_suspect.sum())}"
        f" (incl. {int((tested.epoch_class == 'same_obs').sum())} same-obs eclipse pairs by OBS_DATE)",
        f"- eclipse pairs with confirmed independent epochs: {int((tested.epoch_class == 'independent').sum())}",
        f"\n## Clean weather test (same instrument, shared-data excluded, no DI)",
        f"- pairs: {len(clean)}, significant at FDR {FDR_LEVEL}: {len(sig_clean)}",
        f"- planets with >=1 significant pair: {sig_clean.pl_name.nunique()}",
        f"\n## Cross-instrument (systematics-limited, interpret with care)",
        f"- pairs: {len(cross)}, significant at FDR {FDR_LEVEL}: {len(sig_cross)}",
        f"\n## Top clean anomalies",
    ]
    cols = ["pl_name", "spec_type", "authors_a", "authors_b", "epoch_class", "n", "offset", "chi2_red", "p_value"]
    lines.append(sig_clean[cols].head(20).to_string(index=False))
    (REPORTS / "phase3_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines[:12]))

    # overlay plots for top clean anomalies
    for rank, r in enumerate(sig_clean.head(TOP_PLOTS).itertuples(), 1):
        sa, sb = specs[r.spec_id_a], specs[r.spec_id_b]
        fig, ax = plt.subplots(figsize=(9, 5))
        for s, meta, color in [(sa, (r.authors_a, r.instrument_a), "C0"),
                               (sb, (r.authors_b, r.instrument_b), "C1")]:
            ax.errorbar(s["wave"].to_numpy(), s["value"].to_numpy(),
                        yerr=[s["err_lo"].to_numpy(), s["err_hi"].to_numpy()],
                        fmt="o", ms=3, lw=0.8, alpha=0.8, color=color,
                        label=f"{meta[0]} ({meta[1].split('(')[0].strip()})")
        ax.set_xlabel("wavelength [microns]")
        ax.set_ylabel("depth [%]")
        ax.set_title(f"#{rank} {r.pl_name} [{r.spec_type}] p={r.p_value:.2e}, offset={r.offset:+.4f}")
        ax.legend(fontsize=8)
        fig.tight_layout()
        slug = r.pl_name.replace(" ", "_")
        fig.savefig(REPORTS / f"anomaly_{rank:02d}_{slug}.png", dpi=150)
        plt.close(fig)
    print(f"\nplots for top {min(TOP_PLOTS, len(sig_clean))} anomalies -> reports/")


if __name__ == "__main__":
    main()
