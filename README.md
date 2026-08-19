# exo-spectra

**Homogeneous statistics on every published exoplanet spectrum.**

The NASA Exoplanet Archive collects the published spectra of exoplanet
atmospheres: 1826 spectra of 289 planets (transmission, eclipse, and direct
imaging). Almost all of them were analyzed one paper at a time, with different
tools and different assumptions. This project treats the archive as a single
dataset and runs uniform, model-free statistics across all of it — on a home
PC, using only public data.

## Studies

| # | study | one-line result | release |
|---|---|---|---|
| 1 | [Epoch variability](studies/01-epoch-variability/) | In pairs of spectra that differ only by observation epoch, 17% disagree beyond their stated errors — all around M-dwarf stars. | `study-01-v1.0` |
| 2 | [Blind anomaly search](studies/02-blind-anomaly-search/) | A model-free scan of the real archive: it blindly rediscovers known detections (KELT-9 b Balmer lines), known disputes (GJ 1132 b), and known instrument artifacts — plus a ranked catalog of what remains. | `study-02-v1.0` |

Each study page has four layers: a TL;DR for everyone, a summary for the
technical reader, an expert section, and the result catalogs with column
documentation.

## Repository layout

- `src/` — the shared pipeline: archive downloader, IPAC table parser,
  pair statistics, and the per-study analysis scripts
- `studies/<study>/results/` — result catalogs (CSV) and figures
- `data/spectra_summary.csv` — the base catalog: one row per published spectrum
- `docs/METHODS.md` — the statistical methods in one place
- `docs/REPRODUCE.md` — how to re-run everything from scratch

## Quick start

```bash
pip install -r requirements.txt
python src/download_spectra.py      # fetch all 1826 spectra (~112 MB) from the archive
python src/build_summary.py         # parse and index them
python src/phase3_full_archive.py   # study 1: pair statistics
python src/phase5_features.py       # study 2: structure stats + point anomalies
python src/phase5_cohort.py         # study 2: cohort shape oddballs
```

## License and citation

Code: MIT. Result catalogs and figures: CC BY 4.0. To cite, see
[CITATION.cff](CITATION.cff).

## Acknowledgments

This research has made use of the NASA Exoplanet Archive, which is operated by
the California Institute of Technology, under contract with the National
Aeronautics and Space Administration under the Exoplanet Exploration Program.
The spectra themselves are the work of their original authors — every catalog
row carries the source bibcode.
