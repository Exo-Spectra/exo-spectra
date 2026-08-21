# Study 5 — Morning vs evening: a search for limb asymmetry in JWST white-light transits

*A separate transit depth fitted to the first and second half of 33 JWST/NIRSpec white-light transit light curves of 7 planets around M dwarfs — with injection-calibrated errors — finds zero significant asymmetries and produces the first uniform upper-limit table for this sample.*

## TL;DR (for everyone)

When a planet crosses its star, starlight filters through the ring of
atmosphere around the planet's edge. That ring has two halves with different
weather: a cooler "morning" side and a warmer "evening" side. A warmer side is
more puffed up, so the planet can look slightly bigger when it leaves the star
than when it enters. JWST has seen this effect on a giant, very hot planet
(WASP-39 b). Nobody had checked it systematically for small planets around
small red stars — the targets JWST observes most often. We tested 33 light
curves of 7 such planets, fitting a separate size to the entry and the exit of
each transit. The result: **nothing**. No planet shows a significant
morning-evening difference. Every candidate signal died under closer
inspection — it either shrank with a better statistical treatment or
contradicted itself between two detectors that watched the same transit. A
null result is still a result: we publish, for each light curve, the largest
asymmetry that could have hidden in the noise (median 246 ppm). Anyone
modeling these atmospheres can use those limits.

## Summary (for the technical reader)

We took the 24 JWST/NIRSpec G395H/G395M transit visits already reduced
uniformly in Studies 1 and 3 (the 22 mini-survey visits plus the two GJ 1132 b
pilot visits) and fitted each white-light curve with a limb-darkened transit
model in which the ingress-side and egress-side planet radii are independent.
The statistic is the depth difference **dd = egress − ingress** [ppm]. Visits
that do not cover the full transit with margin are excluded (10 of 43
visit-detector rows, all of L 98-59 b among them), leaving **33 fits of 19
visits, 7 planets**: GJ 1132 b, K2-18 b, LHS 1140 b and c, TOI-1685 b,
TOI-776 b and c.

Result: **0 of 33 fits significant at Benjamini–Hochberg FDR 1%** — on
injection-calibrated errors and on conservative beta-inflated errors alike.

| planet | fits | weighted mean dd [ppm] | note |
|---|---|---|---|
| GJ 1132 b | 4 | +94 ± 32 | largest mean (2.9σ), below the FDR threshold; its strongest single fit on calibrated errors is 2.65σ (p = 0.008), and the one fit escalated to MCMC drops from 3.0σ pre-calibration to 1.7σ |
| K2-18 b | 6 | −52 ± 55 | consistent with zero |
| LHS 1140 b | 3 | +27 ± 99 | consistent with zero |
| LHS 1140 c | 2 | +15 ± 49 | consistent with zero |
| TOI-1685 b | 9 | +41 ± 24 | single fits reach 3.5σ on calibrated errors (at most −2.9σ after MCMC), but the signs flip between visits (consistency p = 0.004) and between detectors of the same visit (up to 3.1σ apart) — systematics, not the planet |
| TOI-776 b | 4 | −1 ± 28 | consistent with zero |
| TOI-776 c | 3 | −0 ± 19 | consistent with zero |

*(The "fits" column sums to 31 of the 33: the 2 spot-flagged fits are
excluded from the per-planet means — see Method A below.)*

Three checks make the null trustworthy. First, errors are calibrated by
**injection-recovery**: we inject known asymmetries into the real light curves
(preserving each visit's red noise) and measure what comes back — recovery is
unbiased (median response slope 1.00) and yields a median per-fit error of
89 ppm. Second, every candidate above 2σ in the first-pass fit (pre-calibration) was
escalated to a full **MCMC posterior**, which killed all of them. Third, a real planetary asymmetry must
repeat across visits and agree between the two detectors (NRS1/NRS2) that
observe the same transit — the only planet with interesting single fits
(TOI-1685 b) fails both tests.

The product is the first uniform upper-limit table on white-light limb
asymmetry for planets around M dwarfs: **median 95% upper limit
|dd| = 246 ppm** per fit, with the best single fits reaching 31–84 ppm
(`results/upper_limits.csv`).

## For the expert

**Sample and light curves.** The 24 visits are the Study 3 mini-survey set
(programs jw02372/jw02722 for K2-18, jw03942 for L 98-59, jw07073 for
LHS 1140, jw03263/jw04195 for TOI-1685, jw02512 for TOI-776) plus the two
GJ 1132 b visits of program jw01981 from the Study 1 JWST annex. White-light
curves come from the same extraction machinery as Study 3
(`src/survey_analyze.py`), per visit and detector (NRS1/NRS2; LHS 1140 is
G395M, NRS1 only), segments concatenated with per-segment out-of-transit
normalization and a segment-boundary step check. Full-transit gate: both
contacts plus a margin of max(0.3 · T14, fixed floor) must lie inside the
visit; the 10 failing rows (all four L 98-59 b visits on both detectors and
K2-18 visit jw02722001001) are catalogued with `partial_flag = 1` and not fit.

**Model** (`src/limb_asym_model.py`). Pure-numpy quadratic-limb-darkening
transit: the occulted flux is integrated radially over the annulus covered by
the planet disk (validated to < 10 ppm against pytransit). Keplerian orbit
with e and ω fixed at archive values, so the genuine ingress/egress *duration*
asymmetry of an eccentric orbit is modelled and cannot leak into the fitted
radius asymmetry (`eccen_assumed = 1` marks planets where the archive gives no
eccentricity and e = 0 is assumed). Limb darkening in the Kipping (2013) q1/q2
parametrization, free. Multiplicative systematics ramp — linear, quadratic, or
quadratic × settling exponential — selected per light curve by BIC on the
symmetric fit. Gaussian priors on a/R* and b from the archive; t0 free
(a pure timing offset cannot masquerade as asymmetry).

**Method A (the headline statistic).** Refit with independent ingress/egress
radii split at the fitted t0, everything else shared;
dd = (k_eg² − k_in²) × 10⁶. First-pass uncertainty: cyclic-shift residual
bootstrap (300 draws; preserves red noise) times a binned-RMS beta factor
(Winn et al. 2008). **Method B (model-independent check):** fold the
symmetric-fit residuals about t0 inside [T1, T4], mirror-interpolate the
egress side onto the ingress side, and test the difference (mean z-score,
chi-squared, sign test). A spot-crossing flag marks in-transit excursions of
the 5-point moving average above 4× its own robust scatter (2 of the 33 fits;
per-planet means are computed without them).

**Injection-recovery calibration** (`src/limb_asym_inject.py`). Per fit:
replace the fitted transit with an asymmetric one at a known dd (mean-depth
preserving radius split), re-add the real residuals cyclically shifted, refit
with the cheap pipeline — 9 amplitudes × 10 noise realizations = 90 fits per
row. Products: `sigma_cal` (the calibrated per-fit error: the scatter of
recovered dd at zero injection), the bias slope, and a 95% upper
limit on |dd|. Calibration verdicts: the cyclic bootstrap is honest (median
sigma_cal/bootstrap = 0.91), the beta inflation double-counts red noise the
bootstrap already carries (median sigma_cal/beta-inflated = 0.59) — so all
reported significances use sigma_cal. A pooled-residual variant changes
sigma_cal by 4% (median), confirming the estimate is not an artifact of
recycling each visit's own residuals.

**MCMC escalation** (`src/limb_asym_mcmc.py`). Every row above 2σ
pre-calibration gets an emcee posterior of the asymmetric model (40 walkers,
1500 burn-in + 3000 steps, free white-noise jitter). Reproducibility gotcha
worth recording: the catalog does not store ramp coefficients, and refitting
the ramp from zero drifts into local minima — the settling-exponential
timescale is partially degenerate with the ingress depth, which produced a
spurious −5.6σ on TOI-776 c in the first pass (visible even on dd ≈ 0 control
rows). Fix: recover the ramp coefficients from the noiseless `ramp` column
persisted with each light curve, with a multi-start over the timescale.
Post-escalation: GJ 1132 b jw01981022001 +1.7σ, TOI-776 c jw02512003001
−2.8σ (also spot-flagged), TOI-1685 b at most −2.9σ (jw03263001001, also
spot-flagged) with inconsistent signs.

**Significance.** Two-sided normal p-values from dd/sigma_cal,
Benjamini–Hochberg FDR at 1% across the 33 fits: zero pass. The same test on
the conservative beta-inflated errors also yields zero.

**Caveats.** White light dilutes the expected signal: limb asymmetry is
chromatic (strongest in molecular bands), so per-wavelength-bin asymmetry —
not done here — is the natural sequel and these limits apply to the
band-integrated G395 (≈2.9–5.1 µm) depth. The stellar limb darkening is symmetric by
construction; an asymmetric stellar disk (spots, gravity darkening) would be
absorbed into dd, which is one more reason the repeatability and NRS1-vs-NRS2
tests matter. Errors inherit the Study 3 white-noise treatment; sigma_cal
absorbs each visit's red noise empirically but assumes it is stationary
within the visit.

## Results & data

| file | contents |
|---|---|
| `results/limb_asym_catalog.csv` | all 43 candidate (visit, detector, planet) rows: symmetric-fit parameters, gates, dd with bootstrap/beta errors, Method B stats, MCMC columns |
| `results/limb_asym_catalog_fdr.csv` | the 33 full fits with injection-calibrated errors, calibrated p-values and the FDR verdict |
| `results/upper_limits.csv` | per-fit sensitivity: sigma_cal, bias slope, 95% upper limit on \|dd\| |
| `results/per_planet_summary.csv` | per-planet inverse-variance weighted mean dd and visit-to-visit consistency |
| `results/lightcurves/*.csv` | the 33 white-light curves (time, flux, fitted models, ramp, residuals) — enough to re-run the injection and MCMC stages without the ~160 GB of JWST downloads |
| `results/summary_dd.png` | all 33 dd measurements with calibrated errors, by planet and detector |
| `results/GJ-1132__jw01981022001_nrs1__GJ_1132_b.png` | example per-fit diagnostic (light curve, detrended transit, residuals, folded ingress vs mirrored egress) |
| `results/inject__GJ-1132__jw01981022001_nrs1__GJ_1132_b.png` | example injection-recovery calibration (recovered vs injected dd) |

Column definitions: [DATA_DICTIONARY.md](DATA_DICTIONARY.md).
Reproduce from scratch: [docs/REPRODUCE.md](../../docs/REPRODUCE.md),
section 8.

**Context reading:** the WASP-39 b limb-asymmetry detection this study is the
small-planet counterpart to: Espinoza et al. 2024,
https://arxiv.org/abs/2407.10294. The two-semicircle transit parametrization
idea: catwoman, Jones & Espinoza 2022,
https://joss.theoj.org/papers/10.21105/joss.02382. Limb-darkening
parametrization: Kipping 2013, https://arxiv.org/abs/1308.0009. Beta red-noise
factor: Winn et al. 2008, https://arxiv.org/abs/0804.4475.
