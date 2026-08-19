# Study 2 — A blind, model-free anomaly search over the real archive

*Unsupervised anomaly detection applied, for the first time, to the real
(not simulated) archive of published exoplanet spectra.*

## TL;DR (for everyone)

We used a computer to find unusual features in 1826 published exoplanet spectra. We did not tell the computer what "strange" meant, and we did not include any physics or chemistry. The computer only looked at numbers and error bars. This blind search found three things: a famous real detection (hydrogen glowing in the ultra-hot planet KELT-9 b), a known instrument artifact (light from a background star leaking into a JWST detector), and the same disputed spectrum that Study 1 found using a different method. The other results are a ranked list of oddities for humans to check.

## Summary (for the technical reader)

Methods to find anomalies in exoplanet spectra exist in the literature (Matchev et al. 2022; "Hunting for Oddballs", 2026). However, those methods only worked on large sets of simulated spectra. Real archives are more complex because every instrument has different wavelength grids, resolutions, and error behaviors. We used a three-tier design to handle this:

- **Tier A — per-spectrum structure.** We tested each of the 741 spectra (≥5 usable points) against three model-free baselines: a constant, a straight line, and a smooth polynomial. **47% of transmission spectra and 66% of eclipse spectra show significant structure** beyond their error bars (FDR 1%). These are mostly real molecular features and serve as the baseline for other tiers. Sanity checks passed: WASP-39 b (deep JWST features) shows strong structure; the flat GJ 1214 b spectrum of Kreidberg et al. 2014 does not (chi2_red 1.02), while the disputed ground-based claims of structure for the same planet are flagged.
- **Tier B — shape oddballs.** We grouped 435 spectra into 7 homogeneous instrument cohorts. We normalized each spectrum to a pure shape and scored it with three outlier measures (PCA reconstruction error, robust Mahalanobis distance, k-nearest-neighbor distance). The top oddball in the Hubble WFC3 cohort is GJ 1132 b, Swain et al. 2021. This is the spectrum at the center of a published dispute about whether an atmosphere exists. Our independent method found it.
- **Tier C — point anomalies.** 243 individual points deviate by more than 4 sigma from their local surroundings. 149 of these points are confirmed by another spectrum of the same planet. Recurring anomalous wavelengths across different planets on the same instrument form 12 "hotspots." These are likely instrument systematics rather than astrophysics.

## For the expert

**Tier A.** Weighted least squares against three nulls: constant (weighted
mean, dof n−1), offset+slope (weighted-centered x, dof n−2), and a weighted
polynomial with degree scaled to n and capped by the number of distinct
wavelengths. Diagnostics: lag-1 autocorrelation of residuals (coherent bands
vs. white noise) and the largest single-point deviation from a running local
median. Benjamini–Hochberg FDR at 1% within each spectrum type.

**Tier B.** Cohorts are defined by instrument regex + wavelength window
(WFC3 G141, STIS, NIRSpec G395, NIRSpec PRISM, NIRISS SOSS, MIRI LRS for
transmission; WFC3 for eclipse). Members are resampled onto a fixed grid
(inverse-variance bin means, ≥80% bin coverage required, cohort-median
imputation for gaps), normalized to zero weighted mean and unit standard
deviation (removing transit-depth scale and scale-height amplitude), then
decomposed by SVD. Outlier scores: reconstruction residual outside the top-k
principal components (k at 90% variance, capped at 8), Mahalanobis distance
with iteratively trimmed center/covariance, and mean distance to the 5
nearest neighbors in PC space. Final rank = median of the three score ranks.
`amp_snr` (shape amplitude over median bin error) flags spectra whose "shape"
is mostly noise.

**Tier C spot checks (top 3 by |z|):** (1) KELT-9 b at 0.656 and 0.486 µm,
z = 97 and 68 — the Balmer lines H-alpha and H-beta; the blind scan
rediscovered the published detection of hydrogen absorption (Cauley et al.
2019). (2) WASP-17 b at 2.05–2.07 µm, z = −25 — consistent with the
zeroth-order background-star contamination that the JWST-TST DREAMS program
itself reports as the dominant systematic for this target. (3) WASP-43 b at
5.25 µm — the blue edge of MIRI LRS, repeating across all spectra of the same
paper: a detector-edge systematic.

**Cross-check against Study 1.** Spearman correlation between a spectrum's
oddball percentile (Tier B) and the maximum chi2_red over its Study-1 pairs is
only **rho = 0.10 (p = 0.058)**: shape oddity and epoch-to-epoch variability
are nearly independent axes of "anomalous" in this archive. The two catalogs
complement each other rather than duplicate.

**Caveats.** Teams that underestimate their uncertainties look artificially
"structured" in Tier A (we flag per-paper clusters). Hotspots mix recurring
real physics (the 1.36–1.38 µm water band on WFC3) with candidate systematics
(e.g. 3.31–3.34 µm on NIRSpec G395H) — they need manual separation. Cohort
imputation slightly damps anomalies in poorly covered bins.

## Results & data

| file | contents |
|---|---|
| `results/phase5_features.csv` | Tier A: per-spectrum structure statistics and flags (741 spectra) |
| `results/phase5_cohort_scores.csv` | Tier B: per-spectrum oddball scores and ranks (435 spectra, 7 cohorts) |
| `results/phase5_point_anomalies.csv` | Tier C: 243 point anomalies with repeatability cross-checks |
| `results/phase5_instrument_hotspots.csv` | 12 instrument × wavelength hotspots |
| `results/phase5_oddball_*.png` | top oddballs vs. their cohort median shape |
| `results/phase5_pca_*.png` | PC1/PC2 shape-space maps per cohort |
| `results/phase5_vs_phase3.png` | oddball rank vs. epoch discrepancy cross-check |
| `results/summary.md` | machine-generated run summary |

Column definitions: [DATA_DICTIONARY.md](DATA_DICTIONARY.md).
Reproduce: [docs/REPRODUCE.md](../../docs/REPRODUCE.md) — Tier A/C:
`python src/phase5_features.py`, Tier B: `python src/phase5_cohort.py`,
report: `python src/phase5_report.py`.
