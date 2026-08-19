# Study 1 — Do repeated spectra of the same exoplanet agree?

*An archive-wide consistency test of published exoplanet spectra across epochs, instruments, and reduction pipelines.*

## TL;DR (for everyone)

When a planet passes in front of its star, the starlight passes through the planet's atmosphere. This light shows a chemical fingerprint called a spectrum. Scientists have published hundreds of these spectra. However, most are single measurements. We took all 1826 published spectra from the NASA Exoplanet Archive. We asked one simple question: **when the same planet was measured twice, do the two measurements agree?** Often they do not. These differences happen mostly around small, active stars. This means we cannot fully trust a single observation of these planets.

## Summary (for the technical reader)

We created the first consistent catalog for published exoplanet spectra. From 1826 spectra of 289 planets, we found 3,991 candidate pairs (same planet, same spectrum type, overlapping wavelengths). Of these, 1462 were testable. For each pair, we adjusted one free vertical offset to remove calibration differences. Then, we tested if the residuals matched the reported error bars (chi-squared test, false-discovery-rate control at 1%). We also classified every pair by provenance: independent epochs vs. re-reductions of the same observation vs. derived products.

Key findings:

1. In the **gold sample** — 47 pairs that differ *only* in observation epoch (same instrument, same reduction pipeline) — **8 pairs (17%) disagree beyond their stated uncertainties**, and every one of them orbits an M-dwarf star (TOI-5205 b, TOI-260 b, TOI-776 b and c, GJ 1132 b). The likely cause is changing star spots, not planetary weather.
2. In 36 pairs where two teams reduced the **same raw data with different pipelines**, results usually agree closely (median chi2_red 0.17) — but for the hardest targets the choice of pipeline alone shifts the spectrum by more than the stated errors.
3. A two-parameter test (offset + slope) splits the discrepancies into two mechanisms: smooth tilts (stellar contamination) and genuine spectral structure differences.
4. Blind validation: with no knowledge of the literature, the pipeline independently recovered famous published controversies (GJ 1132 b, HD 209458 b).

## For the expert

**Method.** Pair statistic: chi-squared of B−A after fitting a free constant
offset (dof = n−1), with an optional offset+slope model (weighted-centered x,
orthogonal parameters, dof = n−2). Wavelength matching: the sparser spectrum
defines the grid; the denser one is averaged inside each bandpass
(inverse-variance weights), with linear interpolation as a flagged fallback.
Asymmetric errors are symmetrized. Multiple testing: Benjamini–Hochberg at
FDR 1%, applied separately to the same-instrument ("clean") and
cross-instrument samples. Full details: [docs/METHODS.md](../../docs/METHODS.md).

**Provenance classification.** A major pitfall in the archive is that many
"pairs" are re-reductions of the same observation, not independent epochs.
We flag statistically-too-consistent pairs (p > 0.999 with n ≥ 8), use per-point
observation dates where available (eclipse spectra), and parse the archive's
`note` field, which labels visits and reduction pipelines. This yields the
pair classes in `pair_results_classified.csv`: epoch_same_pipeline (47, the
gold sample), same_data_diff_pipe (36), epoch_diff_pipeline (169), derived
(153), cross_paper (791), and others.

**Headline numbers.** Clean same-instrument test: 621 pairs → 42 significant
at FDR 1% (15 planets). Gold sample: 8/47 discrepant, all M-dwarf hosts.
With the offset+slope model the gold-sample discrepancies drop from 8 to 4:
for TOI-260 b and TOI-776 b/c a smooth slope absorbs the disagreement
(the textbook signature of unocculted-spot contamination), while GJ 1132 b
and TOI-5205 b keep structural, band-like differences.

**Blind validations.** The anomaly ranking independently surfaced:
GJ 1132 b Swain et al. 2021 vs. Mugnai et al. 2021 (chi2_red 4.5; the
published dispute over an atmosphere detection from the same HST data) and
HD 209458 b Diamond-Lowe et al. 2014 vs. Knutson et al. 2008 (chi2_red 13.8;
the classic thermal-inversion dispute).

**JWST annex.** As an end-to-end check we re-reduced both JWST/NIRSpec G395H
visits of GJ 1132 b ("Double Trouble", May et al. 2023) from `rateints` files
with our own extraction code. Our visit spectra agree with the published
Eureka! reductions (chi2_red 1.22 / p = 0.169 and 1.32 / p = 0.097), and our
visit 1 vs. visit 2 comparison is consistent (n = 35, offset −43 ppm,
chi2_red 0.92, p = 0.60): at R~100, point-by-point with a free offset, the
"Double Trouble" disagreement is a dispute between atmospheric model fits,
not between the spectra themselves — consistent with Bennett et al. 2025.

**Caveats.** We trust the authors' quoted uncertainties; underestimated
correlated noise would show up as discrepancy. Cross-paper pairs have unknown
provenance unless the note field says otherwise. Direct-imaging spectra are
excluded (inconsistent flux units between papers).

## Results & data

| file | contents |
|---|---|
| `results/pair_results.csv` | all 1462 tested pairs: offset, slope, chi2, p-values, flags |
| `results/pair_results_classified.csv` | the same + provenance class per pair |
| `results/pairs_epoch_same_pipeline.csv` | the 47-pair gold sample |
| `results/pairs_same_data_diff_pipe.csv` | the 36 pipeline-comparison pairs |
| `results/anomaly_*.png` | overlay plots of the top discrepant pairs |
| `results/gj1132b_*` | JWST annex: our own GJ 1132 b reductions and comparisons |

Column definitions: [DATA_DICTIONARY.md](DATA_DICTIONARY.md).
Reproduce from scratch: [docs/REPRODUCE.md](../../docs/REPRODUCE.md).
