# Data dictionary — Study 1

## `results/pair_results.csv` and `results/pair_results_classified.csv`

One row per tested pair of spectra (same planet, same spectrum type,
overlapping wavelengths). The `_classified` file adds `pair_class`.

| column | meaning |
|---|---|
| `pl_name` | planet name (NASA Exoplanet Archive convention) |
| `spec_type` | Transmission / Eclipse / Direct Imaging |
| `spec_id_a`, `spec_id_b` | row ids into `data/spectra_summary.csv` |
| `authors_a`, `authors_b` | source papers (may contain HTML entities, e.g. `&ntilde;`) |
| `instrument_a`, `instrument_b` | instruments as listed by the archive |
| `same_instrument` | both spectra from the same instrument (True/False) |
| `epoch_class` | eclipse pairs only, from per-point observation dates: `independent` / `same_obs` / `unknown` |
| `tested` | pair had ≥3 matched points |
| `n` | matched points used |
| `n_interp` | matched points obtained by interpolation fallback (rather than bandpass averaging) |
| `offset`, `offset_err` | fitted vertical offset B−A and its 1σ error (units of the measured value, i.e. depth in %) |
| `chi2`, `dof`, `p_value` | offset-model chi-squared, degrees of freedom (n−1), and tail probability; **small p = discrepant pair** |
| `chi2_red` | chi2 / dof |
| `slope`, `slope_err` | offset+slope model: fitted slope of B−A vs. wavelength (% per micron) |
| `p_slope`, `chi2_red_slope` | tail probability and reduced chi2 of the offset+slope model (dof = n−2); a pair with small `p_value` but large `p_slope` differs by a smooth tilt only |
| `shared_suspect` | flagged as likely re-reduction of the same data (p > 0.999 with n ≥ 8, or same-obs eclipse dates) |
| `significant` | rejected at Benjamini–Hochberg FDR 1% within its sample (same-instrument or cross-instrument) |
| `pair_class` | provenance (classified file only): `epoch_same_pipeline` (gold sample) / `same_data_diff_pipe` / `epoch_diff_pipeline` / `derived` / `same_visit_same_pipe` / `within_paper_other` / `cross_paper` |

## `results/pairs_epoch_same_pipeline.csv`

The 47-pair gold sample (same instrument, same reduction pipeline, different
visits) — subset of the classified file, same columns.

## `results/pairs_same_data_diff_pipe.csv`

The 36 pairs where two pipelines reduced the same observation — same columns.

## `../../data/spectra_summary.csv` (base catalog, shared by both studies)

One row per published spectrum in the archive snapshot (1826 rows).

| column | meaning |
|---|---|
| `spec_id` | stable row id (join key used by all result files) |
| `file` | local filename of the IPAC .tbl spectrum file |
| `pl_name`, `spec_type`, `authors`, `instrument`, `facility`, `bibcode` | metadata copied from the archive index |
| `n_points`, `n_usable` | total points / real measurements (finite value+error, not a limit) |
| `wave_min`, `wave_max` | wavelength coverage [µm] |
| `obs_date_min`, `obs_date_max` | per-point observation date range [JD], where available (mostly eclipse spectra) |
