"""Phase 1: parse all downloaded spectra and build the comparison base.

Outputs (data/processed/):
    spectra_summary.csv  one row per spectrum: id, planet, type, instrument,
                         usable points, wavelength range, obs_date range
    pairs.csv            all comparable pairs: same planet + same spec_type
                         + overlapping wavelength range, with overlap stats
"""
import csv
import itertools
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from spectra_io import load_spectrum, local_name, usable_points

ROOT = Path(__file__).resolve().parents[1]
SPECTRA = ROOT / "data" / "spectra"
OUT = ROOT / "data" / "processed"


def main() -> None:
    index = pd.read_csv(ROOT / "data" / "raw" / "spectra_index.csv")
    OUT.mkdir(parents=True, exist_ok=True)

    rows, errors = [], []
    for i, r in index.iterrows():
        fname = local_name(r["spec_path"])
        try:
            spec = load_spectrum(SPECTRA / fname, r["spec_type"])
        except Exception as e:
            errors.append((fname, str(e)))
            continue
        use = usable_points(spec)
        rows.append({
            "spec_id": i,  # row number in the index = stable id
            "file": fname,
            "pl_name": r["pl_name"],
            "spec_type": r["spec_type"],
            "authors": r["authors"],
            "instrument": r["instrument"],
            "facility": r["facility"],
            "bibcode": r["bibcode"],
            "n_points": len(spec),
            "n_usable": len(use),
            "wave_min": use["wave"].min() if len(use) else None,
            "wave_max": use["wave"].max() if len(use) else None,
            "obs_date_min": use["obs_date"].min() if len(use) else None,
            "obs_date_max": use["obs_date"].max() if len(use) else None,
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "spectra_summary.csv", index=False)
    print(f"parsed {len(summary)}/{len(index)} spectra, {len(errors)} errors")
    if errors:
        with open(OUT / "parse_errors.log", "w", encoding="utf-8") as f:
            for fname, err in errors:
                f.write(f"{fname}\t{err}\n")
        print(f"errors logged to {OUT / 'parse_errors.log'}")

    # comparable pairs: same planet, same type, overlapping wavelength range
    ok = summary[summary["n_usable"] >= 1].copy()
    pair_rows = []
    for (pl, st), grp in ok.groupby(["pl_name", "spec_type"]):
        for a, b in itertools.combinations(grp.itertuples(), 2):
            # inclusive interval overlap (also matches single-point spectra)
            if not (a.wave_min <= b.wave_max and b.wave_min <= a.wave_max):
                continue
            lo = max(a.wave_min, b.wave_min)
            hi = min(a.wave_max, b.wave_max)
            pair_rows.append({
                "pl_name": pl,
                "spec_type": st,
                "spec_id_a": a.spec_id, "spec_id_b": b.spec_id,
                "file_a": a.file, "file_b": b.file,
                "authors_a": a.authors, "authors_b": b.authors,
                "instrument_a": a.instrument, "instrument_b": b.instrument,
                "same_instrument": a.instrument == b.instrument,
                "overlap_lo": lo, "overlap_hi": hi,
                "n_usable_a": a.n_usable, "n_usable_b": b.n_usable,
            })
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(OUT / "pairs.csv", index=False)

    print(f"\ncomparable pairs: {len(pairs)}")
    if len(pairs):
        print(pairs.groupby("spec_type").size().to_string())
        print(f"\nsame-instrument pairs (clean weather test): {pairs['same_instrument'].sum()}")
        print(f"planets with >=1 pair: {pairs['pl_name'].nunique()}")
        top = pairs.groupby("pl_name").size().sort_values(ascending=False).head(10)
        print("\ntop 10 planets by pair count:\n" + top.to_string())


if __name__ == "__main__":
    main()
