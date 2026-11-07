# Erratum (2026-11-07): matched-filter mid-times were biased by ramp systematics

## What was wrong

The O−C values in the original release were measured with the Study 3
matched filter: the transit center was refined by maximizing the in/out
median difference on the *raw* white-light curve, i.e. **before** the ramp
systematics were removed. On curves with a visit-long ramp this estimator is
biased: the ramp tilts the local baseline and drags the recovered center
away from the true mid-transit. A direct comparison against full transit
model fits (ramp fitted explicitly, same `rateints` extraction) shows the
matched-filter centers for TOI-776 were shifted by −13 to −20 minutes,
consistently in all four visits, while for the other targets the shift
stays within ~6 minutes.

## Corrected timing table

Mid-times re-measured with full transit-model fits (quadratic/exponential
ramp selected by BIC, per-detector fits combined by weighted mean; the same
measurements we provided to the ExoClock project). O−C is still quoted
against the NASA Exoplanet Archive default ephemeris retrieved 2026-08,
positive = transit later than predicted.

| planet | visits | mean O−C [min] (corrected) | scatter [min] | original value |
|---|---|---|---|---|
| K2-18 b | 4 | **+49.6** | 5.5 | +51 |
| L 98-59 b | 4 | +0.4 | 0.3 | +4.7 |
| LHS 1140 b | 3 | +2.8 | 0.1 | +7.0 |
| LHS 1140 c | 2 | −3.8 | 0.2 | −2.9 |
| TOI-1685 b | 5 | −0.2 | 0.5 | +2.3 |
| TOI-776 b | 2 | **−1.5** | 1.0 | +12.2 |
| TOI-776 c | 2 | **−2.6** | 0.6 | +16.2 |

## What changes

- **Retracted:** "TOI-776 b and c are 12–16 minutes late." With model fits
  both planets are consistent with the archive ephemeris (−1.5 and −2.6
  minutes, within a few minutes like the rest of the sample). The published
  offset was an artifact of the matched-filter center estimator, not of the
  data.
- **Unchanged:** the headline result. K2-18 b remains ~50 minutes late with
  a stable offset across 2.3 years of visits (+49.6 ± 5.5 min from model
  fits vs +51 min published) — a stale ephemeris, not a timing variation.
- **Unchanged:** the dip search and its conclusion (zero credible new
  transit candidates); it did not use the matched-filter centers.
- The remaining planets shift by ≤ 4 minutes, within the quoted empirical
  precision; no conclusion about them was drawn in the original note.

## Lesson recorded

Transit centers must be measured at the final ephemeris with the ramp
modeled explicitly; an in/out median estimator on a ramp-tilted curve is a
biased clock. The corrected mid-times (per visit, with formal errors from
the fit covariance inflated by the Winn 2008 beta factor) are published in
`results/` of this study as `jwst_midtimes_modelfit.csv`.
