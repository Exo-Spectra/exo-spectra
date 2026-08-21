# Data dictionary — Study 6

## `results/known_transit_timing_oc.csv` (22 rows: one per known transit per visit)

| column | unit | meaning |
|---|---|---|
| `target` | — | MAST target directory name (host star) |
| `visit` | — | JWST visit id (`jwPPPPPOOONNN`) |
| `planet` | — | planet whose transit was timed |
| `pred_bmjd` | BMJD (BJD_TDB − 2400000.5) | predicted mid-transit from the NASA Exoplanet Archive default ephemeris (retrieved 2026-08 via `src/survey_ephem.py`) |
| `obs_bmjd` | BMJD | observed mid-transit: matched-filter refinement of the NRS1 white-light curve (same filter as Study 3) |
| `o_minus_c_min` | minutes | (obs − pred) × 1440; positive = transit arrived later than predicted. No formal per-row error — use the per-planet visit-to-visit scatter (0.5–4.1 min) as the empirical precision |
| `white_depth_ppm` | ppm | matched-filter white-light transit depth (diagnostic; not the calibrated depth of Studies 3/5) |

## `results/dip_candidates.csv` (every SNR ≥ 5 dip outside known-transit windows)

| column | unit | meaning |
|---|---|---|
| `center_bmjd` | BMJD | center of the sliding box at maximum SNR |
| `dur_h` | hours | box duration (0.5, 1.0 or 2.0) |
| `depth_ppm` | ppm | median(out-of-box baseline) − median(in-box), on the quadratically detrended baseline |
| `snr` | — | depth / error, error from baseline point scatter scaled by in/out counts |
| `target`, `visit` | — | where the event was found |

The single row in this catalog (TOI-1685, jw03263001001) is classified as an
edge artifact — see the study README.

## `results/*_lightcurve.png` (one per visit, 22 files)

Detrended NRS1 white light (gray points), 200-bin median curve (blue),
known-planet transit windows (green bands), dip candidates (red dashed
verticals). X axis: hours since visit start.
