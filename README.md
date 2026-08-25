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
| 3 | [JWST mini-survey](studies/03-jwst-mini-survey/) | 22 JWST/NIRSpec transit visits of 7 planets re-reduced from detector frames with one shared code: 14 of 28 visit pairs (50%) disagree beyond errors, clustered on active M dwarfs; a quiet-star control shows none. | `study-03-v1.0` |
| 4 | [L 98-59 b SO₂ verification](studies/04-l9859b-so2-verification/) | The reported ~3σ SO₂ detection re-tested with 54 retrievals (9 spectrum variants × 3 models × 2 codes): it does not survive a change of reduction pipeline, and in our own re-extraction it is largely driven by one anomalous visit. The two retrieval codes agree — the data, not the software, decide. | `study-04-v1.0` |
| 5 | [Limb asymmetry](studies/05-limb-asymmetry/) | A separate ingress/egress depth fitted to 33 JWST white-light transit light curves of 7 M-dwarf planets: zero significant morning/evening asymmetries; injection-calibrated errors give the first uniform upper-limit table for this sample (median 95% limit 246 ppm). | `study-05-v1.0` |
| 6 | [Transit timing & planet hunt](studies/06-transit-timing/) | A by-product sweep of the 22 re-reduced JWST visits: 22 fresh mid-transit times for 7 planets (K2-18 b runs ~51 min late on the archive ephemeris, TOI-776 b/c 12–16 min) and a blind SNR ≥ 5 dip search that finds zero credible new planets. | `study-06-v1.0` |
| 7 | [Chromatic limb asymmetry](studies/07-chromatic-limb-asymmetry/) | The Study 5 transits re-fitted in 8 wavelength bins per detector: no molecular-band-vs-continuum asymmetry contrast passes FDR 1%; two ~3.3σ follow-up candidates (GJ 1132 b CH₄, TOI-776 c CO) survive every systematic test; first per-band upper-limit catalog for small M-dwarf planets (median 95% limit 170 ppm). | `study-07-v1.0` |

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
python src/phase4_classify.py       # study 1: pair provenance classes (gold sample)
python src/phase5_features.py       # study 2: structure stats + point anomalies
python src/phase5_cohort.py         # study 2: cohort shape oddballs
python src/phase5_report.py         # study 2: summary + cross-check vs study 1
```

Study 3 works on JWST detector data (~155 GB of downloads) and has its own
pipeline: see [docs/REPRODUCE.md](docs/REPRODUCE.md), section 6. Study 4 runs
Bayesian retrievals (TauREx + PLATON) on the Study 3 spectra of L 98-59 b:
see section 7 (~1–2 days of CPU time). Study 5 fits asymmetric transits to
the white-light curves of Studies 1 and 3: see section 8 (the committed
light curves let you re-run its statistics without the JWST downloads).
Study 7 repeats the Study 5 fit per wavelength bin: see section 10 (its
committed bin curves plus Study 5's white curves reproduce every fit, again
without the downloads).

Script outputs land in `data/processed/` and `reports/` (both gitignored);
the committed files under `studies/*/results/` are curated copies of those
outputs.

## License and citation

Code: MIT. Result catalogs and figures: CC BY 4.0. To cite, see
[CITATION.cff](CITATION.cff).

## Acknowledgments

This research has made use of the NASA Exoplanet Archive, which is operated by
the California Institute of Technology, under contract with the National
Aeronautics and Space Administration under the Exoplanet Exploration Program.
The spectra themselves are the work of their original authors — every catalog
row carries the source bibcode.
