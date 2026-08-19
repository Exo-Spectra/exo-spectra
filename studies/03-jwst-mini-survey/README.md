# Study 3 — One code, five stars: re-reducing JWST detector data to test visit-to-visit variability

*A uniform re-reduction of 22 public JWST/NIRSpec transit observations of 7 planets, from per-integration detector frames to spectra, with a single extraction code — so that any disagreement between visits cannot come from different reduction choices.*

## TL;DR (for everyone)

Study 1 compared *published* spectra. It found that repeated measurements of the same planet often disagree. However, those published spectra were made by different teams using different software. We wanted to see if the differences were caused by the software or the stars. To test this, we downloaded minimally processed camera frames from the James Webb Space Telescope for 5 well-observed stars (about 155 GB). We processed every observation with **exactly the same code, written by us**. The result: **half of the repeat-visit comparisons still disagree** (14 of 28). These disagreements happen mostly around small, active stars. A quiet star we used as a control showed no disagreements. This means the effect is real, and it comes mostly from the stars (starspots changing between visits), not from differences in reduction software.

## Summary (for the technical reader)

We selected 5 exoplanet hosts with the most public multi-visit JWST/NIRSpec G395H (G395M for LHS 1140) time-series observations: K2-18, L 98-59, LHS 1140, TOI-1685 and TOI-776 — 22 usable transit visits of 7 planets. For each visit, we reduced the public `rateints` files from MAST with one shared extraction code (`src/survey_analyze.py`): column-wise background subtraction, fixed aperture, bad-pixel in-fill from a median template, white-light transit detection, and a low-resolution transmission spectrum (22 wavelength bins per detector). We then compared visit pairs of the same planet using the same statistics as Study 1 (free vertical offset, chi-squared, then an offset+slope model; Benjamini–Hochberg FDR at 1%).

Result: **28 same-planet visit pairs, 14 of them (50%) discrepant beyond the propagated errors.**

| planet | discrepant pairs | max chi2_red | note |
|---|---|---|---|
| L 98-59 b | **6/6** | 4.60 | strongest variability; white-light depths span 815–995 ppm across 4 visits |
| TOI-776 b | **1/1** | 3.25 | significant slope of 98 ppm/µm (6.0σ) — the starspot-contamination signature — but the slope does not absorb all of it |
| TOI-776 c | **1/1** | 2.71 | structural difference; a slope does not help |
| K2-18 b | **4/6** | 2.05 | the planet of the DMS debate: its 4 public visits are not fully consistent point-by-point |
| TOI-1685 b | 2/10 | 1.73 | marginal; most pairs agree |
| LHS 1140 b, c | **0/4** | — | negative control: quiet star, all visits agree (including a pair ~1 year apart) |

Two cross-checks confirm the result. First, TOI-776 b and c were flagged as epoch-variable in Study 1 using *published* spectra. Here, the same planets are discrepant in our independent reduction of the raw frames. Both methods agree. Second, the quiet-star control (LHS 1140) shows no false alarms. Our white-light transit depths also match the published values for these planets.

## For the expert

**Target selection.** MAST query for public NIRSpec BOTS time series
(`instrument_name = "NIRSPEC/SLIT"`, filtered to per-visit `obs_id` entries),
curated to known transiting-planet hosts with ≥2 public G395H/G395M visits
(`src/mast_survey_scout.py`). The five richest hosts were downloaded in full
(`src/survey_download.py`): programs jw02372/jw02722 (K2-18), jw03942
(L 98-59), jw07073 (LHS 1140), jw03263/jw04195 (TOI-1685), jw02512 (TOI-776).

**Extraction** (`src/survey_analyze.py`), per visit and detector (NRS1/NRS2),
segments concatenated: per-column background from rows ≥10 rows away from the
trace; box aperture of ±6 rows around the trace; DQ-flagged pixels in-filled
from a median-image template, inserted unscaled (this step is what makes
box extraction usable); 5σ outlier clipping on each band's time series; linear
detrend fitted on the out-of-transit baseline. Wavelength solution taken from
the public `x1dints` products; the detector-cutout x-offset is recovered from
the calibration-pipeline log embedded in the ASDF extension ("Subarray
x-extents"), since it is not present in the FITS keywords. Timing gotcha worth
recording: `INT_TIMES.int_mid_BJD_TDB` is BJD_TDB − 2400000.5, while the
Exoplanet Archive publishes full JD — mixing them silently misses every
transit window.

**Transit handling.** Predicted mid-times from the archive (`pscomppars`)
refined per visit with a matched box filter (stale ephemerides shift by up to
tens of minutes); in-transit window |t − c| < 0.375·T14, baseline
|t − c| > 0.6·T14; other planets of the same system masked in multi-planet
visits. Depth per wavelength bin = 1 − median(in-transit) / median(baseline)
of the detrended band flux; the per-bin uncertainty is the white-noise error
propagated from the out-of-transit scatter. Bins with errors above
max(200 ppm, 3× the visit's median bin error) are dropped (box extraction
degrades at the red, low-flux end).

**Comparison statistic** identical to Study 1 (docs/METHODS.md): free-offset
chi-squared, then offset+slope; BH-FDR at 1% across the 28 pairs. Headline
pairs: L 98-59 b jw03942002001 vs jw03942004001 chi2_red 4.60
(p = 6.7e-20, n = 40); TOI-776 b jw02512006001 vs jw02512005001 chi2_red 3.25
(p = 6.0e-11) with a 6.0σ slope of 98 ppm/µm that reduces but does not remove
the discrepancy (chi2_red 2.35 with slope); TOI-776 c jw02512003001 vs
jw02512004001 chi2_red 2.71 (p = 7.0e-8), slope-insensitive.

**Sanity anchors.** White-light depths: K2-18 b ≈ 3000 ppm, LHS 1140 b
≈ 5800 ppm, TOI-776 b 1124–1195 ppm, TOI-776 c 1146–1154 ppm — consistent
with published values. The visit flagged most discrepant for L 98-59 b
(jw03942004001) is discrepant against all three other visits independently.

**Caveats.** Uncertainties are white-noise only (no correlated-noise/systematics
model), so chi2_red values are upper bounds on significance — but the LHS 1140
control suggests the inflation is mild. No limb-darkening model (box transit
window; a constant-depth approximation inside a conservative window). Spectral
resolution (22 bins per detector) is below the published R~100. A single free
offset per pair absorbs any visit-level normalization difference, so the
flagged discrepancies are shape changes, not depth-scale changes.

## Results & data

| file | contents |
|---|---|
| `results/survey_visits.csv` | the 22 extracted visit-spectra (target, visit, planet, bins, visit start time) |
| `results/survey_pairs.csv` | all 28 visit pairs: offset, slope, chi2, p-values, FDR flag |
| `results/spectra/*.csv` | each visit's transmission spectrum (wave, dwave, depth_pct, err_pct) |
| `results/<target>__<planet>.png` | per-planet overlay of all visit spectra |

Column definitions: [DATA_DICTIONARY.md](DATA_DICTIONARY.md).
Reproduce from scratch: [docs/REPRODUCE.md](../../docs/REPRODUCE.md)
(warning: ~155 GB of JWST detector-frame downloads).
