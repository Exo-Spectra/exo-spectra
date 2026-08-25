# Data dictionary — Study 7

Common keys: `target` (MAST directory), `planet`, `visit` (`jwPPPPPOOONNN`),
`det` (nrs1/nrs2), `bin` (0–7, per-detector wavelength bin index),
`wave_lo_um`/`wave_hi_um` (bin edges [µm]), `band` (CH4 / CO2 / CO / cont —
by bin center against CH₄ 3.20–3.45, CO₂ 4.20–4.45, CO 4.50–4.75 µm).
dd = egress − ingress transit depth [ppm], measured on the divide-white
corrected curve, i.e. the **chromatic contrast relative to the band
average** (see study README).

## `results/limb_chrom_fdr.csv` (239 calibrated bins)

| column | meaning |
|---|---|
| `dd_ppm`, `dd_err_boot_ppm` | per-bin dd and its cyclic-bootstrap error (raw, uncalibrated) |
| `k_in`, `k_eg` | fitted ingress/egress radius ratios |
| `depth_ppm` | symmetric-fit depth of the bin |
| `ramp_model` | baseline chosen by BIC (linear / quad / quad_exp) |
| `sigma_ppm` | per-point white scatter of the corrected bin curve |
| `chi2_red`, `n_pts` | asymmetric-fit reduced chi-squared and point count |
| `spot_flag` | 1 = the white-light fit was spot-flagged in Study 5 |
| `dd_white_ppm` | the white-light dd of the parent fit (Study 5 catalog) |
| `sigma_cal_ppm`, `sigma_resid_ppm` | injection-calibrated error (zero-amplitude scatter; pooled variant) |
| `bias_slope` | recovered-vs-injected dd slope |
| `z_cal`, `p_cal` | dd/sigma_cal and its two-sided p-value |
| `fdr_pass` | Benjamini–Hochberg FDR 1% verdict across the 239 bins |
| `ul95_ppm` | \|dd\| + 1.645·sigma_cal |

## `results/stack_planet_bin.csv` (89 rows: planet × detector × bin)

`n_fits` visits stacked; `dd_ppm`, `err_ppm` = inverse-variance weighted
mean and its error; `z`, `p`, `fdr_pass` as above.

## `results/band_contrast.csv` (19 rows: planet × molecular band)

`dd_band_ppm`/`err_band_ppm` and `dd_cont_ppm`/`err_cont_ppm` = weighted
means over the band and continuum bins of that planet; `contrast_ppm` =
band − continuum with `contrast_err_ppm`, `contrast_z`, `p`, `fdr_pass`.

## `results/chromatic_structure.csv` (33 rows: one per fit)

Chi-squared of dd(λ) against its own weighted mean: `n_bins`,
`dd_mean_ppm`, `chi2`, `dof`, `p_struct`, `fdr_pass`.

## `results/limb_chrom_calibration.csv` (239 rows)

Injection-recovery per bin: `n_fits` (recovery fits), `sigma_cal_ppm`,
`sigma_resid_ppm`, `bias_slope`, `dd_obs_ppm`, `cal_vs_boot`
(sigma_cal / bootstrap error), `ul95_abs_dd_ppm`.

## `results/limb_chrom_mcmc.csv` (27 escalated bins)

`ramp` (baseline of the fresh refit), `dd_lsq_ppm` (stage-2 value),
`dd_map_ppm` (fresh least-squares start), `dd_mcmc_ppm` /
`dd_mcmc_err_ppm` (posterior median and half the 68% interval; the
posterior width is NOT injection-calibrated — judge against
`sigma_cal_ppm`), `z_cal`, `z_mcmc`.

## `results/ld_sensitivity.csv`, `results/ramp_sensitivity.csv`

The two stress tests on the candidate planets' band bins:
`dd_pinned_ppm` vs `dd_freeLD_ppm` (with the fitted `q1_free`/`q2_free`),
and `dd_stage2_ppm`/`ramp_stage2` vs `dd_fullramp_ppm`/`ramp_full`;
`shift_sigma` = the dd change in units of sigma_cal.

## `results/bins_index.csv` (264 rows)

Bin bookkeeping: `n_cols` (detector columns in the bin), `n_ok` (finite
integrations), `oot_scatter_ppm` (out-of-transit rms per integration).

## `results/lightcurves/<target>__<visit>_<det>__bins.csv` (33 files)

`t_bmjd` (BJD_TDB − 2400000.5) plus `bin00`…`bin07`: the normalized,
segment-renormalized per-bin fluxes BEFORE the divide-white correction.
Combine with Study 5's committed white light curves (flux, model_sym,
ramp) to rebuild the corrected curves and reproduce every fit in this
study without any JWST downloads.
