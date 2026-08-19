"""Phase 2: prototype the pair statistic on well-studied planets.

Usage: python src/phase2_case_study.py "GJ 1214 b" [spec_type]

Produces:
    reports/<planet>_overlay.png      all usable spectra overlaid
    reports/<planet>_pairs.csv        per-pair statistics
    stdout summary
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from compare import compare_pair
from spectra_io import load_spectrum, usable_points

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "reports"


def main(planet: str, spec_type: str = "Transmission") -> None:
    summary = pd.read_csv(ROOT / "data" / "processed" / "spectra_summary.csv")
    pairs = pd.read_csv(ROOT / "data" / "processed" / "pairs.csv")
    sel = summary[(summary.pl_name == planet) & (summary.spec_type == spec_type) & (summary.n_usable > 0)]
    print(f"{planet} [{spec_type}]: {len(sel)} usable spectra")
    REPORTS.mkdir(exist_ok=True)

    # load all spectra once
    specs = {}
    for r in sel.itertuples():
        specs[r.spec_id] = usable_points(load_spectrum(ROOT / "data" / "spectra" / r.file, spec_type))

    # overlay plot
    fig, ax = plt.subplots(figsize=(11, 6))
    for r in sel.itertuples():
        s = specs[r.spec_id]
        label = f"{r.authors} ({r.instrument.split('(')[0].strip()}, n={len(s)})"
        ax.errorbar(s["wave"].to_numpy(), s["value"].to_numpy(),
                    yerr=[s["err_lo"].to_numpy(), s["err_hi"].to_numpy()],
                    fmt="o", ms=3, lw=0.8, capsize=0, alpha=0.7,
                    label=label if len(sel) <= 14 else None)
    ax.set_xscale("log")
    ax.set_xlabel("wavelength [microns]")
    ax.set_ylabel("transit depth [%]" if spec_type == "Transmission" else "eclipse depth [%]")
    ax.set_title(f"{planet} — all published {spec_type.lower()} spectra")
    if len(sel) <= 14:
        ax.legend(fontsize=7, loc="best")
    slug = planet.replace(" ", "_")
    fig.tight_layout()
    fig.savefig(REPORTS / f"{slug}_overlay.png", dpi=150)
    print(f"overlay -> reports/{slug}_overlay.png")

    # pair statistics
    psel = pairs[(pairs.pl_name == planet) & (pairs.spec_type == spec_type)]
    rows = []
    for p in psel.itertuples():
        if p.spec_id_a not in specs or p.spec_id_b not in specs:
            continue
        res = compare_pair(specs[p.spec_id_a], specs[p.spec_id_b])
        if res is None:
            continue
        rows.append({
            "authors_a": p.authors_a, "authors_b": p.authors_b,
            "instrument_a": p.instrument_a, "instrument_b": p.instrument_b,
            "same_instrument": p.same_instrument,
            "n": res.n, "n_interp": res.n_interp,
            "offset": res.offset, "offset_err": res.offset_err,
            "chi2": res.chi2, "dof": res.dof, "p_value": res.p_value,
            "chi2_red": res.chi2 / res.dof,
        })
    out = pd.DataFrame(rows).sort_values("p_value")
    out.to_csv(REPORTS / f"{slug}_pairs.csv", index=False)
    print(f"\npairs tested: {len(out)} (of {len(psel)} candidates; rest had <3 matched points)")
    if len(out):
        print(f"consistent (p>=0.01): {(out.p_value >= 0.01).sum()}, discrepant (p<0.01): {(out.p_value < 0.01).sum()}")
        print(f"  of discrepant, same-instrument: {((out.p_value < 0.01) & out.same_instrument).sum()}")
        cols = ["authors_a", "authors_b", "same_instrument", "n", "offset", "chi2_red", "p_value"]
        print("\nmost discrepant pairs:")
        print(out[cols].head(8).to_string(index=False))
        print("\nmost consistent pairs:")
        print(out[cols].tail(3).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "GJ 1214 b",
         sys.argv[2] if len(sys.argv) > 2 else "Transmission")
