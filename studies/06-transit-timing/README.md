# Study 6 — Clock check: 22 fresh JWST transit timings and a blind hunt for extra planets

*A by-product sweep of the 22 uniformly re-reduced JWST/NIRSpec visits: a
matched-filter mid-transit time for every known transit (7 planets — the
archive ephemeris of K2-18 b runs ~51 minutes late), and a blind SNR ≥ 5 dip
search over the out-of-transit baseline that finds no credible new planets.*

## TL;DR (for everyone)

Our survey light curves (Studies 3 and 5) contain more information than the
spectra we built them for. This study harvests two free by-products. First, a
**clock check**: every transit we recorded has a precise arrival time, and we
compared it with the prediction of the public ephemeris that observers use to
schedule telescopes. Most planets arrive within a few minutes of schedule —
but **K2-18 b, one of the most famous small planets in the sky, arrives about
51 minutes late**, consistently across four visits spanning over two years. Its
public ephemeris simply has not been refreshed since the K2 era; anyone
scheduling expensive telescope time from it should know. Second, a **planet
hunt**: we scanned the hours of flat light curve around each known transit
for unexplained dips — the signature of an unknown planet crossing the star.
Result: nothing credible. One formal blip passes the statistical bar but sits
at the very edge of its observation window, where such blips are a known
instrumental artifact — and it never repeats in the other four visits of the
same star. Not a planet.

## Summary (for the technical reader)

Input: the 22 JWST/NIRSpec G395H/M visits of 5 M-dwarf systems (K2-18,
L 98-59, LHS 1140, TOI-1685, TOI-776) extracted uniformly in Study 3
(NRS1 white light). For every transit of a known planet predicted inside a
visit we refine the mid-transit time with the Study 3 matched filter and
report **O−C against the NASA Exoplanet Archive default ephemeris**
(retrieved 2026-08); positive O−C = transit later than predicted. The
out-of-transit baseline (known windows masked with padding) is then
quadratically detrended and box-scanned at 0.5/1/2 h durations for
unexplained dips at SNR ≥ 5.

**Timing (22 mid-times, 7 planets):**

| planet | visits | mean O−C [min] | visit-to-visit scatter [min] |
|---|---|---|---|
| K2-18 b | 4 | +51 | 3.1 |
| L 98-59 b | 4 | +4.7 | 1.8 |
| LHS 1140 b | 3 | +7.0 | 3.0 |
| LHS 1140 c | 2 | −2.9 | 4.1 |
| TOI-1685 b | 5 | +2.3 | 2.5 |
| TOI-776 b | 2 | +12.2 | 2.0 |
| TOI-776 c | 2 | +16.2 | 0.5 |

The K2-18 b offset is stable (no drift) across 2.3 years of visits — a stale
ephemeris, not a timing variation we could attribute to an unseen perturber.
TOI-776 b and c are 12–16 minutes late. The rest are within ~10 minutes.

**Dip search:** one formal candidate survives the SNR ≥ 5 cut in 22 visits
(TOI-1685, visit jw03263001001: 213 ppm, 0.5 h, SNR 5.2). It sits in the
final 30 minutes of the ~19 h visit, at the very edge of the data — where the
sliding box loses its trailing baseline and edge systematics accumulate — and
it does not repeat in the four other TOI-1685 visits. We classify it as an
instrumental artifact, not a planet. **Zero credible new transit
candidates.**

## For the expert

**Pipeline** (`src/survey_planet_hunt.py`, reusing the Study 3 extraction
from `src/survey_analyze.py`). Per visit, NRS1 `rateints` segments only —
white light needs no second detector. Known-planet windows are predicted
from `src/survey_ephem.py` output (NASA Exoplanet Archive TAP, default
ephemeris per planet, retrieved 2026-08) and refined with the same matched
filter used by Study 3 to place its transit windows; the refined center and
matched-filter white depth per (visit, planet) are the timing catalog.

**O−C caveats.** The matched filter reports no formal per-time uncertainty;
use the per-planet visit-to-visit scatter (0.5–4.1 min) as the empirical
precision. Two rows have O−C ≈ 0 to floating-point precision — the refinement
grid landed exactly on the prediction; they are genuine measurements at the
grid step, not copies. O−C mixes ephemeris staleness with any real timing
variation; within our baselines the offsets are constant per planet, so
staleness is the parsimonious reading. The K2-18 b default archive ephemeris
(Sarkis et al. 2018) is anchored to a 2015 Spitzer transit epoch (Benneke et
al. 2017); over ~11 years even a small period error accumulates to tens of
minutes, which is what we see. These mid-times are suitable for ephemeris refresh
(e.g. an ExoClock submission) but are **not** a TTV search: one white-light
time per visit, minutes-level precision.

**Dip search.** Baseline = points outside every gap-padded known window;
quadratic detrend in time; box widths 0.5/1/2 h, each slid in steps of a
quarter of its own duration;
depth = median(out) − median(in), error from the point scatter of the
detrended baseline scaled by the in/out counts; SNR ≥ 5 keeps the strongest
non-overlapping events. Sensitivity is therefore to single-visit,
hour-scale transits of roughly ≥ 150–250 ppm depending on the visit — an
Earth-sized planet around these M dwarfs would produce ~280–1800 ppm
depending on the host, at or above that threshold; the search
is blind to periods, it only sees what crossed during the ~3–19 h visits.

**Relation to Studies 3 and 5.** Study 5's asymmetric fits also produce a
free t0 per (visit, detector), but only for the 33 full transits it fits
(L 98-59 b excluded there); this catalog covers all 22 visits (L 98-59 b and
the partial K2-18 visit included), one time per visit, with a single
consistent method — use this one for timing, Study 5's for shape.

## Results & data

| file | contents |
|---|---|
| `results/known_transit_timing_oc.csv` | the 22 timings: predicted and observed mid-transit [BMJD], O−C [min], matched-filter white depth [ppm] |
| `results/dip_candidates.csv` | every dip passing SNR ≥ 5 outside known windows (1 row — the TOI-1685 artifact) |
| `results/*_lightcurve.png` | per-visit diagnostic: detrended white light, binned curve, known-transit windows (green), dip candidates (red dashed) |

Column definitions: [DATA_DICTIONARY.md](DATA_DICTIONARY.md).
Reproduce from scratch: [docs/REPRODUCE.md](../../docs/REPRODUCE.md),
section 9.

**Context reading:** the ExoClock project (ephemeris maintenance for Ariel),
Kokori et al. 2022, https://arxiv.org/abs/2110.13863. K2-18 b discovery
and validation: Montet et al. 2015, https://arxiv.org/abs/1503.07866; the
Spitzer transit anchoring its archive ephemeris: Benneke et al. 2017,
https://arxiv.org/abs/1610.07249.
