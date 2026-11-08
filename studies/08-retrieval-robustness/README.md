# Study 8 — How robust are atmospheric retrieval conclusions to upstream choices? A 144-cell two-code sensitivity matrix

*8 science pairs of published or re-reduced transmission spectra of the same
planet — each pair differing in exactly one upstream choice (data reduction,
epoch, extraction pipeline, instrument) — plus one consistent-visit negative
control, fit with identical model classes in two independent radiative-transfer
codes (TauREx 3 and PLATON) by full nested sampling. Headline: swapping only
the reduction of the same HST data erases a capped-7σ H₂O detection to 0.0σ
in both codes (GJ 1132 b, the Swain/Mugnai dispute quantified); honest
error inflation for epoch-to-epoch disagreement drives every marginal
K2-18 b CH₄ cell to 0.0σ; strong detections survive everything.*

## TL;DR (for everyone)

When astronomers announce "water detected on planet X at 7 sigma," that
number quietly depends on decisions made long before the statistics run:
whose version of the data you start from, which night's observation you
use, which software extracted the spectrum. We built a controlled
experiment: pairs of spectra of the *same planet* that differ in exactly
one of those choices, pushed through the exact same analysis in two
independent codes. The result is sobering. For one famous small planet
(GJ 1132 b), the water detection is as strong as our scale allows in one
team's version of the data — and completely absent in another team's
version *of the same telescope data*. For K2-18 b, the methane hints
vanish entirely once you admit, in the error bars, that the planet's
spectrum does not repeat between visits. Genuinely strong detections (hot
giants like HD 189733 b) survive every swap. The lasting product is the
full grid of 144 model comparisons (120 nested-sampling runs, 118 fully
converged) — a ground-truth set anyone can use to test faster retrieval
methods.

## Summary (for the technical reader)

**Design.** 8 science pairs + 1 negative control (9 pairs in the matrix),
each differing along exactly one axis, drawn from our own catalogs
(Study 1 pair catalog, Study 3 mini-survey re-reductions, Study 4
variants):

| pair | axis | side a vs side b |
|---|---|---|
| GJ 1132 b | reduction | Swain et al. 2021 vs Mugnai et al. 2021 (same HST/WFC3 data) |
| WASP-12 b | reduction | Changeat et al. 2022 vs Swain et al. 2013 (HST/WFC3) |
| HD 189733 b | reduction | Fu et al. 2024 vs Zhang et al. 2025 (same JWST G395H observation) |
| L 98-59 b | reduction | FIREFLy vs Eureka published reductions (Study 4 runs reused) |
| TOI-5205 b | pipeline | same data, two reduction pipelines (Cañas et al. 2026) |
| TOI-5205 b | epoch | different visits, same pipeline (Cañas et al. 2026) |
| K2-18 b | epoch | our homogeneous reduction of visits jw02722001001 vs jw02372005001 |
| WASP-39 b | instrument | Feinstein et al. 2023 (NIRISS) vs Carter et al. 2024 (combined) |
| LHS 1140 b | **negative control** | two visits our pipeline found consistent (chi²_red 0.75, p = 0.76) |

Each side is fit with a flat model and single-gas models (H₂O / CH₄ / CO₂ /
SO₂ per the pair's science case, plus grey-cloud variants) in **TauREx 3
and PLATON**, by nested sampling (dynesty, budgeted and checkpointed).
K2-18 b adds an error-inflation variant ("infl"): uncertainties inflated by
the measured epoch disagreement from Study 3's catalog. Significance uses
the Benneke & Seager (2013) Bayes-factor mapping, **capped at 6.99σ**
(the p-value solver saturates near p ~ 1e-12; lnB is always reported
alongside). The matrix has 144 cells, 72 per code. They come from 120
distinct nested-sampling runs. The flat baseline does not depend on the
radiative-transfer code, so for 16 of the 18 flat baselines (spectrum ×
error variant) both codes have an identical flat run (same input and seed,
hence the same result); for LHS 1140 b a/b the TauREx flat run predates the
sampler switch (see `DATA_DICTIONARY.md`) and differs. TOI-5205 b side a is
the same spectrum in two pairs, and 9 of its 10 cells are identical runs
(PLATON H₂O was run separately). `run_id` in `results/evidences_all.csv`
marks every identical run; treat them as one, not as independent samples.
118 distinct runs ended on the dlogz = 0.1 stopping criterion (19 of them
predate the `capped` flag, but all stopped below every call budget ever in
force, so they ended on dlogz). 2 runs hit the call budget (100 000 at the
time) and are flagged `capped = 1`: TauREx CH₄ on both sides of
HD 189733 b. Their lnB vs flat is about −99 and −98, and even their best
likelihood found is about 81 below the flat model, so the 0σ verdict
there (R7) is very unlikely to depend on the unfinished part of the evidence
integral (the best likelihood found is a lower bound, not a strict limit).
*(Correction in v1.1: v1.0 stated "72 runs per code, 144 total, all
converged".)*
(`results/sigma_all.csv` holds 116 rows = 108 non-flat fits — 144 cells
minus 36 flat cells, which enter only as the lnB reference — plus 8
L 98-59 b rows imported from Study 4.)

**Results.**

* **R1 — a reduction swap alone can erase a detection.** GJ 1132 b H₂O:
  6.99σ (Swain) → 0.00σ (Mugnai) in both codes, with or without clouds;
  the posterior moves by z = 10.7 in planet radius and z = 5.6 in
  log-abundance (both PLATON; TauREx 3.8 and 0.4). Nothing warns you —
  each side is internally confident.
* **R2 — or merely maim it.** WASP-12 b H₂O: 6.99 → 4.10σ (TauREx),
  6.42 → 3.61σ (PLATON) between reductions.
* **R3 — the L 98-59 b SO₂ case reproduces inside the matrix**: FIREFLy
  1.33 / 0.00σ (TauREx/PLATON) vs Eureka 3.19 / 3.10σ — same answer as
  Study 4, now as one row of a homogeneous grid.
* **R4 — epochs + honest errors kill marginal claims.** K2-18 b CH₄ on raw
  errors is at most marginal (PLATON 0.87 / 1.19σ, peaking at 1.48σ with
  cloud; TauREx 0.00 everywhere); with epoch-inflated errors **every
  K2-18 b cell in both codes is 0.00σ**.
* **R5 — robust detections exist (with one asterisk each).** HD 189733 b
  H₂O/CO₂ holds 6.99 ↔ 6.99 across reductions in both codes — though its
  PLATON posteriors still shift by up to z = 21.6 in temperature.
  TOI-5205 b holds 6.99 ↔ 6.99 on the epoch axis in both codes and in
  TauREx on the pipeline axis; in PLATON the pipeline swap costs up to
  1.47σ (CO₂, 6.99 → 5.52), all cells staying above 5σ. The negative
  control is clean: LHS 1140 b = 0.00σ in all 12 cells (posterior z-shifts
  between consistent visits ≤ 3.8 — the noise baseline for the z values
  above).
* **R6 — the instrument axis moves parameters, not verdicts.** WASP-39 b
  H₂O is overwhelming on both sides (lnB vs flat 583 and 927 in PLATON,
  ~1030 in TauREx), but the posteriors disagree at z = 16.5 in planet
  radius. "Detected in both" ≠ "the same atmosphere in both".
* **R7 — the one qualitative code disagreement is a lesson.** HD 189733 b
  CH₄: TauREx 0.00σ vs PLATON 6.99σ on both sides (|ΔlnB| 328–407), yet
  both codes pin CH₄ at the prior floor — neither sees methane. The model
  classes differ physically: PLATON includes H₂–H₂/H₂–He collision-induced
  absorption and its continuum beats the flat model by Δloglike ≈ +319;
  our TauREx wrapper has no CIA and its saturated CH₄ line cores
  anti-correlate with the data, driving the sampler to the temperature
  prior floor (483.7 K). **Takeaway: "Nσ of gas X vs flat" on hot giants
  partly measures continuum physics, not the gas.** Within-code Δσ along
  an axis — the quantity this study reports — remains valid, because both
  sides of a pair share the same model class.

## For the expert

* Model classes are deliberately minimal (isothermal single gas + N₂
  filler, optional grey cloud) so that both sides of every pair, and both
  codes, share identical parameterizations; absolute abundances are not a
  product of this study — Δσ and posterior z-shifts along controlled axes
  are.
* The 6.99σ cap is a property of the Benneke & Seager p-value inversion
  (brentq saturates at p ≈ 1e-12); treat all 6.99 entries as "≥ 6.99" and
  use `lnB_vs_flat` for anything quantitative.
* The TauREx wrapper lacks a CIA contribution (see R7); a full CIA re-run
  would not change any within-code conclusion and is documented rather
  than fixed. The R7 probe diagnostics quoted in prose (continuum
  amplitude ~175 ppm, correlations r = 0.64 / −0.30) are internal probe
  numbers, not reproducible from the shipped CSVs; everything else in R7
  (prior-floor pinning, Δloglike, |ΔlnB|) is.
* The K2-18 b inflation factor comes from our Study 3 epoch catalog; other
  defensible inflation choices shift R4 quantitatively, not qualitatively.
* Prior art: catalogue-level pipeline comparisons (e.g. MNRAS 2024) and
  simulation-based challenges (ABC/Ariel). To our knowledge this is the
  first homogeneous two-code sensitivity grid along
  reduction/epoch/pipeline/instrument axes on real archival spectra of the
  same planets. The absolute dynesty lnZ of all 144 cells, with the exact
  input spectra and priors (added in v1.1), double as a ground-truth set
  for validating ML-accelerated retrievals (FMPE + importance sampling,
  transformer priors).

## Results & data

* `results/sigma_all.csv` — every (pair, side, code, model, err_variant):
  lnB vs flat, σ, best-fit and posterior parameters.
* `results/sensitivity_pairs.csv` — per pair-axis: σ_a, σ_b, Δσ and
  posterior z-shifts (radius, T, log-abundance, cloud pressure).
* `results/codes_agreement.csv` — TauREx vs PLATON per cell: Δσ, ΔlnB.
* `results/pairs.csv` — the pair definitions (files, axes, authors).
* `results/evidences_all.csv` *(v1.1)* — all 144 cells with the absolute
  lnZ, lnZ_err, max log-likelihood, likelihood calls, live points,
  `capped` flag and `run_id` (shared runs have the same id).
* `results/input_spectra/<pair>__<side>.csv` *(v1.1)* — the 16 spectra
  exactly as fitted: cut to the pair's common wavelength range, in ppm.
  For K2-18 b `infl` cells, multiply `err_ppm` by `err_scale` from
  `retrieval_setup.csv`.
* `results/retrieval_setup.csv` *(v1.1)* — per spectrum and error variant:
  system parameters, model wavelength range, model list and the exact
  uniform prior bounds.
* Validation of the v1.1 export (`src/gap7_export_release.py`): the flat
  model's lnZ has a closed form; recomputed from the spectra as fitted, for
  all 36 flat cells, it matches dynesty within 1.3 lnZ_err everywhere
  (max |difference| 0.23), and the same check from the published CSVs gives
  the same result. The script itself runs only in the working repository
  (it reads the intermediate retrieval files); its `flat_lnz()` can be
  reused on the published CSVs.
* `results/delta_sigma.png` + `results/fit__<pair>.png` — overview heatmap
  and per-pair best fits for the 8 newly-run pairs (L 98-59 b reuses
  Study 4 artifacts).
* Reproduce: `src/gap7_retrieval.py` (the matrix; expects the pair spectra
  listed in `results/pairs.csv`) and `src/gap7_report.py` (aggregation;
  imports `l9859b_retrieval.py` from this repository). Software versions:
  TauREx 3.3.2, PLATON 5.4, dynesty 3.1.0. TauREx cross-sections (R = 15000,
  TauREx HDF5 format): H₂O POKAZATEL, CH₄ MM, CO₂ UCL-4000, SO₂ ExoAmes
  (ExoMol line lists); PLATON uses its own bundled opacity data. Spectra
  sources:
  NASA Exoplanet Archive Atmospheric Spectroscopy table (archive pairs)
  and Study 3's homogeneous re-reductions (K2-18 b, LHS 1140 b; their
  `source_file` lives in `studies/03-jwst-mini-survey/results/spectra/`).

## Version history

* **v1.1 (2026-11-08)** — data-only release. Adds `evidences_all.csv`,
  `input_spectra/` and `retrieval_setup.csv`, so the matrix can be used as a
  benchmark without re-running anything. Corrects the run count and the
  convergence statement (120 distinct runs; 2 capped at the call budget, see
  Summary) and documents the sampler and budget history. No result, figure
  or conclusion changes.
* **v1.0 (2026-10-11)** — first release.
