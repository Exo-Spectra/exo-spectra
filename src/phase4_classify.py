"""Phase 4: classify tested pairs using the index `note` field.

JWST-era papers publish per-visit spectra and per-pipeline reductions
separately, labeled in `note` (e.g. "Eureka! reduction; visit 2").
This lets us split pairs into:

    epoch_same_pipeline   same paper, different visit, same pipeline
                          -> GOLD: pure epoch-to-epoch (weather/activity) test
    same_data_diff_pipe   same paper, same visit, different pipeline
                          -> measures reduction-pipeline systematics
    epoch_diff_pipeline   same paper, different visit AND different pipeline
    derived               one side is a coadd/joint/average (shares data)
    cross_paper           different papers (provenance unknown without
                          reading them; may be new data or re-reduction)
    within_paper_other    same paper, notes not parseable

Outputs: data/processed/pair_results_classified.csv + stats to stdout.
"""
import html
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from phase3_full_archive import bh_fdr, FDR_LEVEL

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed"

DERIVED_KEYS = ["coadd", "joint", "average", "weighted mean", "combined"]
PIPELINES = [("eureka", "Eureka!"), ("exotic", "ExoTiC-JEDI"), ("jedi", "ExoTiC-JEDI"),
             ("tiberius", "Tiberius"), ("firefly", "FIREFLy")]


def parse_note(note):
    if not isinstance(note, str):
        return None, None, False
    n = note.lower()
    derived = any(k in n for k in DERIVED_KEYS)
    m = re.search(r"visit\s*(\d+)", n)
    visit = m.group(1) if m else None
    pipe = next((label for k, label in PIPELINES if k in n), None)
    return visit, pipe, derived


def classify(row, notes):
    na, nb = notes.get(row.spec_id_a), notes.get(row.spec_id_b)
    bca, bcb = notes_bib.get(row.spec_id_a), notes_bib.get(row.spec_id_b)
    if bca != bcb:
        return "cross_paper"
    va, pa, da = na
    vb, pb, db = nb
    if da or db:
        return "derived"
    if va and vb:
        if va != vb and pa == pb and pa:
            return "epoch_same_pipeline"
        if va == vb and pa != pb and pa and pb:
            return "same_data_diff_pipe"
        if va != vb:
            return "epoch_diff_pipeline"
        return "same_visit_same_pipe"
    return "within_paper_other"


index = pd.read_csv(ROOT / "data" / "raw" / "spectra_index.csv")
index["pl_name"] = index["pl_name"]  # spec_id = row number (as in build_summary)
notes = {i: parse_note(r) for i, r in index["note"].items()}
notes_bib = index["bibcode"].to_dict()

res = pd.read_csv(OUT / "pair_results.csv")
res["pair_class"] = [classify(r, notes) for r in res.itertuples()]
res["authors_a"] = res["authors_a"].map(lambda s: html.unescape(s) if isinstance(s, str) else s)
res["authors_b"] = res["authors_b"].map(lambda s: html.unescape(s) if isinstance(s, str) else s)

tested = res[res.tested & (res.spec_type != "Direct Imaging")].copy()

print("pair classes (tested, no DI):")
print(tested.pair_class.value_counts().to_string())

for cls, label in [("epoch_same_pipeline", "GOLD epoch pairs (weather test)"),
                   ("same_data_diff_pipe", "same-data pipeline-vs-pipeline")]:
    sub = tested[tested.pair_class == cls].copy()
    if not len(sub):
        continue
    sub["significant"] = bh_fdr(sub.p_value.to_numpy(), FDR_LEVEL)
    nsig = int(sub.significant.sum())
    print(f"\n== {label}: {len(sub)} pairs, significant @FDR1%: {nsig} "
          f"({100*nsig/len(sub):.0f}%), planets: {sub.pl_name.nunique()}, "
          f"median chi2_red {sub.chi2_red.median():.2f}")
    cols = ["pl_name", "spec_type", "n", "offset", "chi2_red", "p_value"]
    print(sub.sort_values("p_value")[cols].head(12).to_string(index=False))
    sub_out = sub.sort_values("p_value")
    sub_out.to_csv(OUT / f"pairs_{cls}.csv", index=False)

res.to_csv(OUT / "pair_results_classified.csv", index=False)
print("\nsaved pair_results_classified.csv")
