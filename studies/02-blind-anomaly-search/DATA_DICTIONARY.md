# Data dictionary — Study 2

`spec_id` in all files joins to `data/spectra_summary.csv` (see Study 1
dictionary for the base-catalog columns).

## `results/phase5_features.csv` (Tier A — one row per spectrum, 741 rows)

| column | meaning |
|---|---|
| `n_used` | usable points with finite positive errors |
| `chi2_red_flat`, `p_flat` | reduced chi2 and tail probability vs. a constant; **small p = the spectrum has structure** |
| `chi2_red_slope` | reduced chi2 vs. the offset+slope model |
| `slope_snr` | absolute fitted slope over its 1σ error |
| `chi2_red_smooth`, `poly_deg` | reduced chi2 vs. a weighted polynomial and its degree |
| `acf_lag1`, `acf_lag1_snr` | lag-1 autocorrelation of slope-model residuals (and ×√n); high values = coherent, band-like residuals |
| `max_abs_z_local` | largest single-point deviation from the running local median, in units of that point's σ |
| `n_pts_gt3sig` | number of points with local deviation > 3σ |
| `struct_amp` | excess scatter beyond the quoted errors, in data units (depth %) |
| `structured` | significant structure at Benjamini–Hochberg FDR 1% within its spectrum type |

## `results/phase5_cohort_scores.csv` (Tier B — one row per cohort member, 435 rows)

| column | meaning |
|---|---|
| `cohort` | cohort name (e.g. `T_WFC3_G141` = transmission, Hubble WFC3, 1.10–1.66 µm) |
| `coverage` | fraction of the cohort's wavelength bins covered by this spectrum (≥0.8 required; the rest imputed) |
| `amp_snr` | shape amplitude / median bin error; low values = the "shape" is mostly noise |
| `score_recon` | PCA reconstruction residual outside the top-k components |
| `score_mahal` | robust Mahalanobis distance in PC space |
| `score_knn` | mean distance to the 5 nearest neighbors in PC space |
| `oddball_rank` | median of the three score ranks within the cohort; **1 = most anomalous** |
| `pca_k` | number of principal components used for this cohort |

## `results/phase5_point_anomalies.csv` (Tier C — one row per anomalous point, 243 rows)

| column | meaning |
|---|---|
| `wave`, `value`, `sigma` | wavelength [µm], measured value (depth %), quoted 1σ |
| `z_local` | deviation from the spectrum's own running local median, in units of σ |
| `n_other_specs` | other spectra of the same planet & type that cover this wavelength (±2%) |
| `n_confirming` | of those, how many deviate the same way (same sign, \|z\| > 2) |
| `n_contradicting` | of those, how many deviate significantly the other way |

## `results/phase5_instrument_hotspots.csv` (12 rows)

Instrument × wavelength bins (1% wide, logarithmic) collecting ≥3 anomalies
from ≥2 different planets — candidate instrument/reduction systematics.

| column | meaning |
|---|---|
| `n_anomalies`, `n_planets`, `n_spectra` | anomaly count and how many distinct planets/spectra contribute |
| `wave_lo`, `wave_hi` | wavelength range of the contributing anomalies [µm] |
| `mean_z` | mean signed deviation (systematics tend to a consistent sign) |

Caveat: recurring *real* features also produce hotspots (e.g. the 1.36–1.38 µm
water band on WFC3) — separating physics from systematics needs a human.
