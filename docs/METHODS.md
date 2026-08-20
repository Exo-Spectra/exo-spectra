# Methods

The statistics of Studies 1–3 are model-free: no atmospheric models, no
assumed chemistry. Inputs are the published data points and their quoted
uncertainties. Study 4 is the exception by design: a targeted verification
that uses Bayesian model comparison with atmospheric forward models
(last section).

## Data

NASA Exoplanet Archive, TAP table `spectra` (Atmospheric Spectroscopy):
1826 spectra, 289 planets. Each spectrum is an IPAC table with central
wavelength, bandwidth, measured value (transit depth [%] / eclipse depth [%] /
F_lambda), asymmetric uncertainties, and a limit flag. We use only real
measurements (finite value, at least one finite error, limit flag = 0) and
symmetrize asymmetric errors by averaging the two sides. Spectrum types are
never mixed; direct-imaging spectra are excluded from population statistics
(flux units are not consistent between papers).

## Pair consistency test (Study 1)

For two spectra A, B of the same planet and type with overlapping wavelengths:

1. **Common grid.** The sparser spectrum defines the grid. For each of its
   points, the denser spectrum's points inside the bandpass
   [w − dw/2, w + dw/2] are averaged with inverse-variance weights; if the
   bandpass is empty or undefined, linear interpolation is used and flagged
   (`n_interp`). At least 3 matched points are required.
2. **Offset model.** chi2 = Σ w_i (B_i − A_i − c)^2 with w_i = 1/(σ_A² + σ_B²),
   c fitted analytically, dof = n − 1. p = P(chi2 ≥ observed).
3. **Offset + slope model.** B − A = c + m·(λ − ⟨λ⟩_w) with weighted-centered
   wavelength (orthogonal parameters), dof = n − 2. A pair that is discrepant
   under the offset model but consistent under offset+slope has a smooth,
   tilted difference — the signature of changing stellar contamination.
   A pair discrepant under both has structural (band-like) differences.
4. **Multiple testing.** Benjamini–Hochberg FDR at 1%, applied separately to
   same-instrument and cross-instrument samples.

**Pair provenance.** Pairs are classified before interpretation:
- *shared-data suspects*: p > 0.999 with n ≥ 8 (statistically too consistent),
  or overlapping per-point observation dates (eclipse spectra);
- the archive `note` field labels visits and reduction pipelines
  (Eureka!, ExoTiC-JEDI, Tiberius, FIREFLy) and derived products
  (co-adds, joint fits, averages), yielding classes:
  `epoch_same_pipeline` (gold sample), `same_data_diff_pipe`,
  `epoch_diff_pipeline`, `derived`, `cross_paper`.

Known limitation: the shared-data flag catches re-reductions that *agree*;
re-reductions that *disagree* (the interesting disputes) appear as anomalies
and must be resolved from the literature.

## Structure statistics (Study 2, Tier A)

Per spectrum with ≥5 usable points, weighted least squares against three nulls:
constant (dof n−1), offset+slope (dof n−2), and a weighted polynomial of degree
3/2/1 depending on n (capped by the number of distinct wavelengths; fit on a
centered, scaled wavelength variable). Diagnostics: lag-1 autocorrelation of
slope-model residuals (times √n for significance) and per-point deviations from
a centered running median (window 5), in units of the point's own σ.
BH-FDR at 1% on the flat-model p-values, per spectrum type.

## Cohort shape analysis (Study 2, Tier B)

Cohorts = instrument regex + wavelength window with a fixed bin grid.
Per member: inverse-variance bin means; ≥80% bin coverage required;
missing bins imputed with the cohort median shape (flagged via `coverage`).
Shape normalization: subtract the weighted mean, divide by the standard
deviation of the binned values — removing transit-depth scale and
scale-height amplitude. PCA via SVD on the column-centered shape matrix;
k components at 90% explained variance (max 8). Outlier scores:

1. reconstruction residual outside the top-k subspace,
2. Mahalanobis distance with iteratively trimmed (3 × 10%) center/covariance,
3. mean distance to the 5 nearest neighbors in PC space.

Final rank = median of the three score ranks. `amp_snr` = shape amplitude /
median bin error distinguishes signal-shaped from noise-shaped members.

## JWST detector-frame extraction (Study 3)

Inputs are the public `rateints` (per-integration slope images) and `x1dints`
(pipeline 1-D spectra, used only for the wavelength solution) files from MAST,
per visit and detector (NRS1/NRS2), long time series concatenated across
segments.

1. **Background.** Per column, the median of rows ≥ 10 rows away from the
   spectral trace is subtracted.
2. **Aperture.** Box extraction of ±6 rows around the trace.
3. **Bad pixels.** DQ-flagged pixels are replaced from a median-image template
   (inserted unscaled) — without this, box extraction is dominated by
   missing-pixel noise.
4. **Cleaning and detrending.** 5σ clipping on each band's time series; a
   linear trend fitted on the out-of-transit baseline is divided out.
5. **Transit window.** Predicted mid-times from the archive (`pscomppars`) are
   refined per visit with a matched box filter. In-transit:
   |t − c| < 0.375·T14; baseline: |t − c| > 0.6·T14; other planets of the same
   system are masked. Depth per bin = 1 − median(in-transit) / median(baseline)
   of the detrended band flux; uncertainty propagated from the out-of-transit
   scatter (white noise only).
6. **Binning and quality.** 22 wavelength bins per detector; bins with errors
   above max(200 ppm, 3× the visit's median bin error) dropped (box extraction
   degrades at the red end).

Practical notes: the detector-cutout x-offset needed to place the wavelength
solution is only available in the calibration-pipeline log inside the ASDF
extension ("Subarray x-extents"); and `INT_TIMES.int_mid_BJD_TDB` is
BJD_TDB − 2400000.5 while the archive publishes full JD.

Visit pairs of the same planet are then compared exactly as in Study 1
(offset model, offset+slope model, BH-FDR at 1%).

## Bayesian model comparison (Study 4)

A verification matrix for one published detection claim: every spectrum
variant × every atmosphere model × two independent retrieval codes, with
identical data, priors and sampler settings.

1. **Models** (isothermal, well-mixed): `flat` — one parameter, a constant
   transit depth; single-gas models — three parameters (planet radius,
   temperature, log surface pressure; uniform priors), 99.9% of the tested
   gas + 0.1% N₂, absorption + Rayleigh scattering, no clouds. System
   parameters (stellar radius and temperature, planet mass) are fixed at the
   published values. The setup deliberately mirrors the model comparison of
   the paper under test rather than adding physics.
2. **Codes.** TauREx 3 (ExoMol R = 15 000 cross-sections) and PLATON 5.4
   (native opacity grid, per-gas mixtures via `custom_abundances`). Forward
   models are averaged over each data point's bandpass — the same matching
   rule as the Study 1 pair statistic.
3. **Evidence.** Static nested sampling (dynesty), fixed live-point count and
   random seed — every cell deterministic. Model preference is
   ln Bayes factor vs `flat`, mapped to n-sigma via Benneke & Seager (2013).
4. **Controls.** A flat synthetic spectrum with the real error bars (the
   codes must not detect a molecule on it), and the `flat`-model lnZ — which
   is independent of the radiative-transfer code — must agree between codes
   for every spectrum (it does, to better than 1e-10).
5. **Error-inflation axis.** Where multiple visits of the same planet are
   averaged, a second variant inflates each bin's error by √max(1, χ²_red)
   of the visit-to-visit scatter — propagating measured epoch inconsistency
   into the averaged error bars.

## Point anomalies and hotspots (Study 2, Tier C)

Points with |z| > 4 against their spectrum's running local median enter the
catalog. Each anomaly is cross-checked for repeatability: other spectra of the
same planet and type covering that wavelength (±2%) either confirm
(same sign, |z| > 2) or contradict it. Anomalies are then aggregated per
instrument in 1%-wide logarithmic wavelength bins; bins collecting ≥3
anomalies from ≥2 different planets are "hotspots" — candidate instrument or
reduction systematics rather than astrophysics (though recurring real features,
like the 1.4 µm water band, also produce hotspots and must be separated by hand).
