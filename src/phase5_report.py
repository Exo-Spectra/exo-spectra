"""Phase 5: assemble the summary report from Tier A/B/C outputs.

Also cross-checks Tier B oddball ranks against Phase 3/4 pair discrepancies:
if a spectrum is shape-anomalous within its cohort, is it also the one that
disagrees with other epochs of the same planet? (Spearman rank correlation
between the per-spectrum oddball percentile and the max chi2_red over the
pairs it participates in.)

Outputs: reports/phase5_summary.md, reports/phase5_vs_phase3.png
"""
import html
import sys
from pathlib import Path

sys.stdout.reconfigure(errors="replace")  # cp1250 console vs unescaped names (Cañas)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from phase3_full_archive import FDR_LEVEL

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed"
REPORTS = ROOT / "reports"


def unescape(df: pd.DataFrame, col: str = "authors") -> pd.DataFrame:
    df = df.copy()
    df[col] = df[col].map(lambda s: html.unescape(str(s)))
    return df


def main() -> None:
    feats = unescape(pd.read_csv(OUT / "phase5_features.csv"))
    anoms = unescape(pd.read_csv(OUT / "phase5_point_anomalies.csv"))
    hs = pd.read_csv(OUT / "phase5_instrument_hotspots.csv")
    scores = unescape(pd.read_csv(OUT / "phase5_cohort_scores.csv"))
    pairs = pd.read_csv(OUT / "pair_results.csv")

    # ---- cross-check vs phase 3: oddball percentile vs max pair chi2_red ----
    scores["oddball_pct"] = scores.groupby("cohort")["oddball_rank"].transform(
        lambda r: 1.0 - (r - 1) / (len(r) - 1) if len(r) > 1 else 0.5)
    tested = pairs[pairs.tested & pairs.chi2_red.notna()]
    long = pd.concat([
        tested[["spec_id_a", "chi2_red"]].rename(columns={"spec_id_a": "spec_id"}),
        tested[["spec_id_b", "chi2_red"]].rename(columns={"spec_id_b": "spec_id"}),
    ])
    max_chi2 = long.groupby("spec_id").chi2_red.max().rename("max_pair_chi2_red")
    xc = scores.merge(max_chi2, on="spec_id", how="inner")
    rho, pval = stats.spearmanr(xc.oddball_pct, np.log10(xc.max_pair_chi2_red))

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(xc.oddball_pct, xc.max_pair_chi2_red, s=14, alpha=0.5)
    ax.set_yscale("log")
    ax.set_xlabel("Tier B oddball percentile within cohort (1 = most anomalous)")
    ax.set_ylabel("max chi2_red over phase-3 pairs")
    ax.set_title(f"Shape anomaly vs epoch discrepancy: Spearman rho={rho:.2f} (p={pval:.1e})")
    fig.tight_layout()
    fig.savefig(REPORTS / "phase5_vs_phase3.png", dpi=150)
    plt.close(fig)

    # ---- summary markdown ------------------------------------------------
    lines = ["# Phase 5 — model-independent anomaly search over the archive\n"]

    lines.append("## Tier A — per-spectrum structure statistics")
    for st in ("Transmission", "Eclipse"):
        f = feats[feats.spec_type == st]
        lines.append(f"- {st}: {len(f)} spectra (>=5 usable points), "
                     f"structured @FDR{FDR_LEVEL}: {int(f.structured.sum())} "
                     f"({f.structured.mean():.0%}), median chi2_red_flat "
                     f"{f.chi2_red_flat.median():.2f}")
    top_struct = feats.nlargest(15, "chi2_red_flat")
    cols = ["pl_name", "spec_type", "instrument", "authors", "n_used",
            "chi2_red_flat", "chi2_red_slope", "acf_lag1_snr"]
    lines.append("\nTop 15 by structure amplitude (chi2_red vs flat):\n")
    lines.append(top_struct[cols].to_string(index=False))

    lines.append(f"\n## Tier C — point anomalies (|z_local| > 4)")
    lines.append(f"- {len(anoms)} anomalous points in {anoms.spec_id.nunique()} spectra; "
                 f"{int((anoms.n_confirming > 0).sum())} confirmed by >=1 other spectrum "
                 f"of the same planet, {int((anoms.n_contradicting > 0).sum())} contradicted")
    lines.append(f"- instrument hotspots (>=3 anomalies from >=2 planets in a 1% "
                 f"wavelength bin): {len(hs)} — recurring bins are suspected "
                 f"instrument/reduction systematics, not astrophysics\n")
    acols = ["pl_name", "spec_type", "instrument", "authors", "wave", "z_local",
             "n_other_specs", "n_confirming", "n_contradicting"]
    lines.append("Top 15 point anomalies by |z|:\n")
    lines.append(anoms.head(15)[acols].to_string(index=False))
    lines.append("\nInstrument hotspots:\n")
    lines.append(hs.to_string(index=False))

    lines.append(f"\n## Tier B — cohort shape oddballs ({scores.cohort.nunique()} cohorts, "
                 f"{len(scores)} spectra)")
    for c, g in scores.groupby("cohort"):
        top = g.nsmallest(3, "oddball_rank")
        entries = "; ".join(f"{r.pl_name} ({r.authors}, amp_snr {r.amp_snr:.1f})"
                            for r in top.itertuples())
        lines.append(f"- **{c}** (n={len(g)}): {entries}")

    lines.append(f"\n## Cross-check vs phase 3 pair test")
    lines.append(f"- spectra in both analyses: {len(xc)}; Spearman(oddball percentile, "
                 f"log max pair chi2_red) = {rho:.2f} (p = {pval:.1e})")
    lines.append("- plot: phase5_vs_phase3.png")

    (REPORTS / "phase5_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines[:14]))
    print(f"\nwrote reports/phase5_summary.md; cross-check rho={rho:.2f} p={pval:.1e}")


if __name__ == "__main__":
    main()
