# Data dictionary — Study 5

All times are BMJD_TDB (BJD_TDB − 2400000.5). All depths and depth
differences are in ppm of stellar flux. The asymmetry statistic is
**dd = egress depth − ingress depth** [ppm]; positive dd = the planet looks
bigger on the way out.

## `results/limb_asym_catalog.csv` (one row per candidate visit-detector-planet, 43 rows)

Stage-2 catalog: every transit found in the visits, including the ones that
failed the full-transit gate (no asymmetry fit, `dd_ppm` empty).

| column | meaning |
|---|---|
| `target` | host-star directory name as downloaded from MAST (e.g. `TOI-776`) |
| `planet` | the transiting planet (archive name) |
| `visit` | JWST visit identifier (`jwPPPPPOOOVVV`: program, observation, visit) |
| `det` | detector (`nrs1` / `nrs2`; LHS 1140 = G395M, `nrs1` only) |
| `n_int` | integrations in the white-light curve after clipping |
| `e`, `w_deg` | orbital eccentricity and argument of periastron used (archive values, fixed in the fit) |
| `eccen_assumed` | 1 = archive gives no eccentricity, e = 0 assumed |
| `partial_flag` | 1 = transit contacts + margin not fully inside the visit — row catalogued but **not fit** |
| `spot_flag`, `spot_z` | in-transit excursion of the 5-point moving average above 4× its own robust scatter (starspot-crossing suspect) and its worst z |
| `seg_step_flag`, `seg_step_z` | flux step at a segment boundary detected against the out-of-transit level |
| `t0_bmjd` | fitted mid-transit time (symmetric fit), BMJD_TDB |
| `depth_ppm` | symmetric-fit transit depth k² [ppm] |
| `k`, `aRs`, `b` | fitted radius ratio, scaled semi-major axis, impact parameter |
| `q1`, `q2` | fitted Kipping (2013) limb-darkening parameters |
| `ramp_model` | systematics model selected by BIC (`linear` / `quad` / `quad_exp`) |
| `t14_h` | transit duration T1→T4 from the fitted geometry [hours] |
| `sigma_ppm` | per-integration white scatter [ppm] |
| `chi2_red_sym` | reduced chi-squared of the symmetric fit |
| `bic_linear`, `bic_quad`, `bic_quad_exp` | BIC of the symmetric fit per ramp model |
| `box_depth_ppm` | model-free box depth used to locate the transit |
| `dd_ppm` | **Method A asymmetry statistic**: egress − ingress depth [ppm] |
| `dd_err_boot_ppm` | dd error from the cyclic-shift residual bootstrap (300 draws) |
| `beta` | binned-RMS red-noise inflation factor (Winn et al. 2008 style) |
| `dd_err_ppm` | conservative first-pass error = bootstrap × beta |
| `dd_sigma`, `p_value` | dd over `dd_err_ppm` and its two-sided p-value (the *conservative* significance) |
| `dbic_asym` | BIC(symmetric) − BIC(asymmetric); positive favors the asymmetric model |
| `k_in`, `k_eg` | fitted ingress-side and egress-side radius ratios |
| `mb_n`, `mb_z`, `mb_p_chi2`, `mb_p_sign` | **Method B** folded-residual test: points compared, mean-difference z-score, chi-squared p, sign-test p |
| `dd_ld_sens_ppm` | change in dd when limb darkening is frozen at the symmetric solution (sensitivity check; small = the asymmetry does not trade against limb darkening) |
| `dd_mcmc_ppm`, `dd_mcmc_err_ppm` | emcee posterior median and 1σ for dd (only rows escalated: \|dd\| > 2σ pre-calibration) |

## `results/limb_asym_catalog_fdr.csv` (one row per full fit, 33 rows)

The 33 rows of the catalog that passed the full-transit gate, extended with
the injection-calibrated significance. Same columns as above, plus:

| column | meaning |
|---|---|
| `sigma_cal_ppm` | **the calibrated dd error** from injection-recovery (see `upper_limits.csv`) |
| `sigma_cal_resid_ppm` | pooled-residual calibration variant (robustness check) |
| `ul95_abs_dd_ppm` | 95% upper limit on \|dd\| for this fit |
| `dd_err_cal_ppm` | error used for the final significance (= `sigma_cal_ppm`) |
| `dd_sigma_cal`, `p_value_cal` | dd over the calibrated error and its two-sided p-value |
| `significant` | discrepant from zero at Benjamini–Hochberg FDR 1% (on `p_value_cal`, across the 33 fits) — **False for every row** |

## `results/upper_limits.csv` (one row per full fit, 33 rows)

Injection-recovery sensitivity per fit: asymmetric transits at 9 known dd
amplitudes × 10 red-noise-preserving realizations = 90 refits.

| column | meaning |
|---|---|
| `target`, `planet`, `visit`, `det` | as in the catalog |
| `n_fits` | successful injection refits (of 90) |
| `sigma_cal_ppm` | scatter of recovered dd at zero injection = the calibrated error |
| `sigma_cal_resid_ppm` | variant: scatter around the fitted linear response, pooled over all amplitudes (~9× more samples) |
| `bias_slope`, `bias_intercept_ppm` | linear response of recovered vs injected dd (slope 1, intercept 0 = unbiased) |
| `dd_obs_ppm` | the observed dd of this fit (copied from the catalog) |
| `cal_vs_boot` | sigma_cal / bootstrap error (≈ 1 = the bootstrap was honest) |
| `cal_vs_beta` | sigma_cal / (bootstrap × beta) error (< 1 = beta double-counts red noise) |
| `ul95_abs_dd_ppm` | 95% upper limit on \|dd\|: \|dd_obs\| + 1.645 × sigma_cal |

## `results/per_planet_summary.csv` (one row per planet, 7 rows)

| column | meaning |
|---|---|
| `planet` | planet name |
| `n` | fits combined (spot-flagged fits excluded when ≥ 2 clean fits exist) |
| `dd_w_ppm`, `dd_w_err_ppm` | inverse-variance weighted mean dd and its error (calibrated errors) |
| `chi2`, `dof`, `p_consistency` | visit-to-visit consistency of dd around the mean; small p = the fits disagree with each other (systematics signature) |
| `dd_w_sigma` | weighted mean over its error |

## `results/lightcurves/*.csv` (one file per full fit, 33 files; one row per integration)

Named `<target>__<visit>_<det>__<planet>.csv`. Enough to re-run the
injection-recovery and MCMC stages without the raw JWST downloads.

| column | meaning |
|---|---|
| `t_bmjd` | integration mid-time, BMJD_TDB |
| `flux` | normalized white-light flux (clipped, segment-normalized) |
| `model_sym` | best symmetric model (transit × ramp) |
| `model_asym` | best asymmetric model (transit × ramp) |
| `ramp` | the ramp (systematics) component alone — noiseless, used to recover ramp coefficients exactly |
| `res_sym` | flux − model_sym |
