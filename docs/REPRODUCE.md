# Reproducing the results

Studies 1–2 run on a normal PC (no GPU needed; ~1 GB RAM; minutes, not hours).
Study 3 also runs on a normal PC but needs ~155 GB of disk for JWST detector
data and benefits from ≥16 GB RAM. Study 4 needs the retrieval codes and
their opacity data (section 7) and ~1–2 days of CPU time. Study 5 reuses the
Study 3 downloads — or, for its statistics, just the light curves committed
in this repository (section 8). Python ≥ 3.11.

Note on paths: the scripts write their outputs to `data/processed/` and
`reports/` (both gitignored). The files committed under `studies/*/results/`
and `data/spectra_summary.csv` are curated copies of exactly those outputs.

```bash
pip install -r requirements.txt
```

## 1. Get the data (~112 MB)

```bash
python src/download_spectra.py
```

Downloads all 1826 spectrum files (IPAC .tbl) from the NASA Exoplanet Archive
into `data/spectra/`. The script is idempotent (skips existing files) and
handles the archive's session-token workspace URLs automatically. The spectrum
index is fetched from the TAP service:

```bash
curl "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+*+from+spectra&format=csv" -o data/raw/spectra_index.csv
```

Note: the archive is a living database — a later index may contain more
spectra than the 1826 analyzed here (index snapshot: 2026-08-13).

## 2. Build the base catalog and pair list

```bash
python src/build_summary.py
```

Writes `data/processed/spectra_summary.csv` (one row per spectrum) and
`data/processed/pairs.csv` (3,991 candidate pairs).

## 3. Study 1 — epoch variability

```bash
python src/phase3_full_archive.py   # pair statistics + FDR + top-anomaly plots
python src/phase4_classify.py       # provenance classes from the archive note field
```

Outputs: `pair_results.csv`, `pair_results_classified.csv`, gold-sample CSVs,
`reports/phase3_summary.md`, overlay plots.

## 4. Study 2 — blind anomaly search

```bash
python src/phase5_features.py       # Tier A structure stats + Tier C point anomalies
python src/phase5_cohort.py         # Tier B cohort shape oddballs + plots
python src/phase5_report.py         # summary + cross-check vs Study 1
```

Outputs: `phase5_features.csv`, `phase5_point_anomalies.csv`,
`phase5_instrument_hotspots.csv`, `phase5_cohort_scores.csv`,
`reports/phase5_summary.md`, oddball/PCA plots.

## 5. JWST annex (optional, large downloads)

`src/jwst_whitelight.py` and `src/jwst_spectrum.py` reduce JWST/NIRSpec G395H
`rateints` files (public, from MAST) to white-light curves and transmission
spectra. They were used for the GJ 1132 b end-to-end validation in Study 1.
Expect multi-GB downloads per visit; paths are configured at the top of each
script.

## 6. Study 3 — JWST mini-survey (large downloads: ~155 GB)

```bash
python src/mast_survey_scout.py     # find public multi-visit G395H/M hosts (astroquery)
python src/survey_download.py       # fetch rateints/x1dints for the 5 chosen hosts
python src/survey_ephem.py          # transit ephemerides from the archive TAP service
python src/survey_analyze.py        # extraction + visit-pair comparison (all targets)
```

Requires `astroquery` and `astropy` (in `requirements.txt`). The downloader is
idempotent (verifies cached files by size) and retries on connection drops;
expect the MAST transfer to take many hours. `survey_analyze.py` accepts
target names as arguments to process a subset (e.g.
`python src/survey_analyze.py TOI-776`). Outputs land in `reports/survey/`;
the committed files under `studies/03-jwst-mini-survey/results/` are curated
copies of those outputs.

## 7. Study 4 — L 98-59 b SO₂ verification matrix

Extra dependencies (not in `requirements.txt`, study 4 only):

```bash
pip install taurex==3.3.2 platon==5.4 dynesty==3.1.0
```

Opacity data: TauREx needs ExoMol R = 15 000 cross-sections in TauREx HDF5
format for SO₂ (ExoAmes) and CO₂ (UCL-4000), placed in
`data/opacities/taurex_ready/` (download from exomol.com, "Cross sections
(TauREx)"; ~365 MB each). PLATON downloads its own opacity package on first
use.

Inputs: the Study 3 extraction must exist for L 98-59 (`reports/survey/`),
plus the two published spectra from the archive download (step 1).

```bash
python src/l9859b_prepare.py            # unify the 9 spectrum variants
python src/l9859b_retrieval.py --code taurex   # 27 retrievals (~1 day CPU)
python src/l9859b_retrieval.py --code platon   # 27 retrievals (~1 day CPU)
python src/l9859b_report.py             # sigma matrix + heatmap + fit plots
```

The retrieval script is idempotent (finished cells are skipped on restart;
`--max-runs N` limits a batch — useful for running in fresh processes) and
the sampler is seeded, so results are deterministic for the pinned package
versions. Outputs land in
`reports/l9859b_verification/`; the committed files under
`studies/04-l9859b-so2-verification/results/` are curated copies.

## 8. Study 5 — limb asymmetry in white-light transits

Extra dependency (escalation stage only): `pip install emcee`.

The full chain needs the Study 3 JWST downloads (section 6, ~155 GB — plus
the two GJ 1132 b visits of program 1981 via `src/jwst_whitelight.py` paths)
and the ephemerides from `src/survey_ephem.py`:

```bash
python src/limb_asym_run.py        # stage 2: catalog + per-fit diagnostics
python src/limb_asym_inject.py     # stage 3a: injection-recovery calibration
python src/limb_asym_mcmc.py       # stage 3b: emcee escalation of >2 sigma rows
python src/limb_asym_report.py     # FDR verdict, per-planet table, summary
```

Shortcut without the downloads: the 33 white-light curves are committed under
`studies/05-limb-asymmetry/results/lightcurves/`. Copy them together with
`limb_asym_catalog.csv` into `reports/limb_asymmetry/`, then run
`python src/survey_ephem.py` once — a network-only Exoplanet Archive query
that writes the `data/processed/survey_ephemerides.csv` the injection and
MCMC stages read. After that the injection, MCMC and report stages run as-is
(only `limb_asym_run.py` needs the raw frames). Outputs land in
`reports/limb_asymmetry/`; the committed files under
`studies/05-limb-asymmetry/results/` are curated copies.

## 9. Study 6 — transit timing and blind dip search

Needs the Study 3 downloads (section 6) and the ephemerides from
`src/survey_ephem.py`; no extra dependencies:

```bash
python src/survey_planet_hunt.py            # all survey targets
python src/survey_planet_hunt.py K2-18      # or a subset
```

Outputs land in `reports/planet_hunt/`; the committed files under
`studies/06-transit-timing/results/` are curated copies.

## Tools

`tools/plain_english.py` — the plain-English documentation pass used for this
repository: rewrites Markdown via a local LLM (Ollama) and **fails loudly if
any number changes** between input and output.
