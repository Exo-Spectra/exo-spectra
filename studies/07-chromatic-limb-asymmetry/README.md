# Study 7 — Colors of the terminator: a chromatic limb-asymmetry search in JWST transits of M-dwarf planets

*The 33 white-light transits of Study 5, re-fitted in 8 wavelength bins per
detector with injection-calibrated errors: no morning/evening color
difference passes FDR 1%, two ~3.3σ band candidates survive every
systematic test, and the by-product is the first per-band upper-limit
catalog for small planets around M dwarfs (median 95% limit 170 ppm).*

## TL;DR (for everyone)

Study 5 asked whether a planet looks bigger leaving its star than entering
it — averaged over all colors — and found nothing. But theory says the
effect should hide in specific colors: the wavelengths where molecules like
methane absorb. Averaging all colors together dilutes it. So we cut the
same 33 transit recordings into 8 color channels each and asked the
question again, color by color: is the morning side different from the
evening side *inside the methane band*? *Inside the carbon-dioxide band?*
The strict answer is again **no** — no color passes our significance bar.
But two signals came close and refused to die under every stress test we
threw at them: on **GJ 1132 b** the methane band looks different from the
neighboring colors, and on **TOI-776 c** the carbon-monoxide band does. Both
sit just under the bar — follow-up targets, not discoveries. And as always,
the lasting product is the limit table: for every planet and color channel
we publish the largest asymmetry that could have hidden in the noise.

## Summary (for the technical reader)

Input: the 33 full-transit (visit, detector, planet) white-light fits of
Study 5, each cut into 8 uniform wavelength bins per detector (264 bin
light curves from the same rateints extraction as Study 3; NRS1
2.86–3.72 µm, NRS2 3.82–5.18 µm, the G395M target LHS 1140 spans
2.85–5.19 µm on one detector). Each bin curve is **divide-white** corrected,
so the measured per-bin dd = egress − ingress depth [ppm] is the **chromatic
contrast relative to the band average** — and the band average itself is
null (Study 5). Geometry and limb darkening are pinned to the white
solution; the per-bin baseline is chosen by BIC from linear / quadratic /
quadratic × settling-exponential. Errors are calibrated per bin by
injection-recovery: 239 of 264 bins are calibratable (the 25 dropped are
the low-flux NRS2 red end beyond ~4.8 µm), median sigma_cal = 68 ppm,
recovery unbiased (median response slope 0.996).

**Result: 0 of 19 molecular-band-vs-continuum contrasts significant at
Benjamini–Hochberg FDR 1%** (bands: CH₄ 3.20–3.45, CO₂ 4.20–4.45, CO
4.50–4.75 µm). Two candidates end just below the bar and survive every
systematic test (limb-darkening release, ramp-model menu, MCMC posteriors,
spot-flag exclusion, visit-to-visit sign consistency):

| candidate | contrast | consistency |
|---|---|---|
| GJ 1132 b, CH₄ band | −96 ± 30 ppm (−3.3σ) | same sign in all 4 fits (2 visits × 2 adjacent CH₄ bins — the band lies entirely on NRS1, so there is no cross-detector check here) |
| TOI-776 c, CO band | +148 ± 45 ppm (+3.3σ) | same sign in both visits (CO₂: +94 ± 42, 2.2σ) |

A physical note on the first: the photosphere of a ~3250 K M dwarf carries
CO and H₂O lines but essentially no CH₄, so chromatic stellar contamination
is a poor explanation for a signal localized in the methane band.

Every *individually* significant bin fails repeatability: 4 of the 5
per-bin FDR passers are TOI-1685 b bins whose signs flip between visits and
detectors (its familiar systematics from Studies 3 and 5), and the fifth
(TOI-776 b, CO₂ band, −69 ± 18 ppm) flips to +47 ± 37 ppm in the only other
visit. The one FDR-passing planet-stacked bin (TOI-1685, 3.82–3.99 µm) is
driven by 2 of 5 internally inconsistent visits — the per-visit
consistency tests exist precisely to catch this.

**Product:** the first uniform chromatic limb-asymmetry upper-limit catalog
for small M-dwarf planets — per-bin 95% limits (median 170 ppm), per-planet
stacked dd(λ), and band contrasts (`results/`).

## For the expert

**Extraction** (`src/limb_chrom_extract.py`). Per (visit, detector): the
Study 3 per-column aperture photometry, columns grouped into 8 uniform
wavelength bins; per bin 3×5σ clip, median normalization, and the Study 5
per-segment out-of-transit renormalization with the OOT mask built from the
Study 5 catalog (t0, T14, 0.6·T14 padding) over every planet transiting in
the visit.

**Fits** (`src/limb_chrom_fit.py`). Common mode cm = white_flux ×
white_ramp / white_sym_model from the Study 5 light curves (committed under
`../05-limb-asymmetry/results/lightcurves/`); the bin curve divided by cm
keeps its own transit and loses the shared systematics — including any
band-averaged asymmetry, hence the contrast semantics (a gray,
wavelength-independent asymmetry is invisible here by construction; that is
Study 5's domain). The Study 5 asymmetric model is refit with t0, a/R*, b,
q1, q2 pinned to the white solution by tight Gaussian priors and k_in/k_eg
free; baseline per bin by BIC over the full ramp menu — restricting the
menu to linear/quad biased dd by up to 1.4 sigma_cal in exactly the
interesting bins (the settling-exponential/ingress-depth degeneracy of
Study 5), which is why the released catalog uses the full menu everywhere
(`results/ramp_sensitivity.csv` documents the effect).

**Calibration** (`src/limb_chrom_inject.py`). Per bin: rebuild the corrected
curve, refit sym+asym to recover model, ramp and residuals (the catalog
stores no ramp coefficients — the Study 5 lesson), inject mean-depth
preserving asymmetries at ±(1–3)× the bootstrap error plus 12 zero-amplitude
realizations with cyclically shifted real residuals, refit. sigma_cal =
scatter of recovered dd at zero injection; median response slope 0.996;
sigma_cal/bootstrap median 0.95 (the cyclic bootstrap is honest here).
Bins with bootstrap error ≥ 500 ppm (NRS2 red end, box-extraction failure
mode known from Study 3) are excluded up front.

**Verification ladder.** (1) LD sensitivity (`src/limb_chrom_ldtest.py`):
re-fitting the 26 band bins of the three planets flagged by the first-pass
contrasts (GJ 1132 b, TOI-776 c, LHS 1140 b) with free q1/q2 moves
dd by median 0.02, max 0.29 sigma_cal — pinned LD explains nothing.
(2) Ramp sensitivity (`src/limb_chrom_ramptest.py`): documented above,
folded into the final run. (3) MCMC escalation (`src/limb_chrom_mcmc.py`):
every bin above 2σ calibrated (27 bins) gets an emcee posterior of the
asymmetric model with free white-noise jitter; the posteriors reproduce the
least-squares dd essentially exactly (only two noisy red-end TOI-1685 bins
shrink) — unlike Study 5's white-light escalation, nothing here rested on
a fragile least-squares solution. All quoted significances use sigma_cal.

**Statistics** (`src/limb_chrom_report.py`). Two-sided p-values from
dd/sigma_cal with BH-FDR 1% at four levels: per bin (5/239 pass — all
non-repeatable, above), chromatic structure per fit (dd(λ) vs its weighted
mean: 3/33, all TOI-1685), per-planet stacked bins (1/89, the inconsistent
TOI-1685 stack), and molecular-band vs continuum contrasts per planet
(0/19). Upper limit per bin: |dd| + 1.645·sigma_cal.

**Caveats.** Divide-white pins the flux-weighted sum of the bins near the
white dd, so only wavelength-*dependent* asymmetry is measurable; sigma_cal
assumes noise stationarity within a visit; the two candidates are below the
pre-registered threshold — the honest headline is the upper-limit catalog.

## Results & data

| file | contents |
|---|---|
| `results/limb_chrom_fdr.csv` | all 239 calibrated bins: dd, sigma_cal, calibrated z/p, band, FDR verdict, UL95 |
| `results/stack_planet_bin.csv` | per-planet inverse-variance stacks per (detector, bin) |
| `results/band_contrast.csv` | the 19 band-vs-continuum contrasts with errors and FDR verdicts |
| `results/chromatic_structure.csv` | per-fit dd(λ)-vs-constant chi-squared tests |
| `results/limb_chrom_calibration.csv` | injection-recovery per bin: sigma_cal, bias slope, UL95 |
| `results/limb_chrom_mcmc.csv` | the 27 escalated bins: LSQ vs posterior dd |
| `results/ld_sensitivity.csv`, `results/ramp_sensitivity.csv` | the two systematic stress tests |
| `results/bins_index.csv` | bin wavelength edges, column counts, OOT scatter |
| `results/lightcurves/*.csv` | the 264 per-bin light curves (33 files × 8 bins) — with Study 5's committed white curves these reproduce every fit without the ~160 GB of JWST data |
| `results/dd_lambda__*.png` | stacked dd(λ) per planet with band shading |

Column definitions: [DATA_DICTIONARY.md](DATA_DICTIONARY.md).
Reproduce from scratch: [docs/REPRODUCE.md](../../docs/REPRODUCE.md),
section 10.

**Context reading:** chromatic/limb-resolved asymmetries have so far been
published only for giants — WASP-39 b (Espinoza et al. 2024,
https://arxiv.org/abs/2407.10294), WASP-107 b (Murphy et al. 2024,
https://www.nature.com/articles/s41550-024-02367-9); this study is their
small-planet counterpart.
