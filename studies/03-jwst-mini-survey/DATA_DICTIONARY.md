# Data dictionary — Study 3

All times are BMJD_TDB (BJD_TDB − 2400000.5). All depths and errors are in
percent of stellar flux in the spectra files and in ppm in the pair table
(1% = 10000 ppm).

## `results/survey_visits.csv` (one row per extracted visit-spectrum, 22 rows)

| column | meaning |
|---|---|
| `target` | host-star directory name as downloaded from MAST (e.g. `TOI-776`) |
| `visit` | JWST visit identifier (`jwPPPPPOOOVVV`: program, observation, visit) |
| `planet` | the planet transiting in this visit (archive name) |
| `n_bins` | wavelength bins surviving the quality cut (see README, "For the expert") |
| `t0_bjd` | **visit start** (mid-time of the first integration), BMJD_TDB — NOT the transit mid-time |
| `file` | the spectrum file in `results/spectra/` |

## `results/spectra/*.csv` (one row per wavelength bin)

| column | meaning |
|---|---|
| `wave` | bin central wavelength [µm] |
| `dwave` | bin width [µm] |
| `depth_pct` | transit depth [% of stellar flux] |
| `err_pct` | 1σ white-noise uncertainty [%] |
| `planet` | planet name |
| `det` | detector (`nrs1` / `nrs2`) |

## `results/survey_pairs.csv` (one row per same-planet visit pair, 28 rows)

| column | meaning |
|---|---|
| `target`, `planet` | host and planet |
| `visit_a`, `visit_b` | the two compared visits |
| `n` | wavelength bins common to both spectra |
| `offset_ppm` | fitted constant offset B − A [ppm] |
| `chi2_red`, `p_value` | reduced chi2 and tail probability of the offset model (dof = n − 1); **small p = the two visits disagree** |
| `slope_ppm_um` | fitted slope of B − A vs. wavelength [ppm/µm] (offset+slope model) |
| `slope_sigma` | absolute slope over its 1σ error |
| `chi2_red_slope`, `p_slope` | reduced chi2 and tail probability of the offset+slope model (dof = n − 2); a pair discrepant in `p_value` but consistent here has a smooth, tilted difference (stellar-contamination signature) |
| `significant` | discrepant at Benjamini–Hochberg FDR 1% (on `p_value`, across the 28 pairs) |
