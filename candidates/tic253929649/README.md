# TIC 253929649 -- WITHDRAWN candidate (not supported by the data)

> **Status (2026-11-09): withdrawn by the author.** Do not use the candidate parameters in
> the original report below, and do not treat this star as a planet host on the basis of
> this package. The accompanying Research Note (RNAAS submission AAS82200) has not been
> published; a request to withdraw it was sent to the journal on 2026-11-09.

## Why it was withdrawn

A re-check in November 2026 tested the statistics the note relied on against null
distributions, which the original analysis did not do. Findings 1 and 2 together remove
the claim: the signal is not above noise either with per-series time shifts or at a fixed
ephemeris. Finding 3 explains where the strongest apparent series came from, and finding 4
shows that a blind box search does not single out this period either.

1. **The per-series timing statistic is indistinguishable from noise.** The note's
   procedure searched each light-curve series for the transit-time offset (within ±2 h)
   that maximised the depth SNR. It was run on all 8 archival series (SPOC 2-min sector
   28; TESS-SPOC sectors 1, 28 and 68; QLP sectors 1, 28, 68 and 95) plus our own
   sector-106 TESScut photometry, and kept only the series with SNR >= 2: five archival
   series and sector 106 (`recheck_2026-11-09/note_timing_run.txt`). The re-check repeats
   this on the same 8 archival series, with preprocessing set up to match the note
   (symmetric 4-sigma clip, no QLP flag mask); it reproduces the note's per-series values
   (e.g. QLP sector 95: SNR 4.3 at -38 min) and the same five series with SNR >= 2. At 400
   random orbital phases (0.1 to 0.9 of the period away from the claimed transit), the
   median number of series reaching SNR >= 2 is also 5 (p = 0.67 for at least as many).
   As a second summary statistic, the sum of the best SNR over the 8 series is 18.6,
   against a null distribution with mean 17.6 and standard deviation 2.7 (empirical
   p = 0.39). The maximum SNR over a time window is biased high even on pure noise; without
   a null run it is not evidence of a transit. The sector-106 value (SNR 3.56 at +110 min)
   was not null-tested separately; the sector-106 alias test is discussed below.
2. **At the published linear ephemeris the transit depth is consistent with zero:**
   61 ± 47 ppm (weighted mean over the same 8 series, with cadences flagged 2^30, 2048 or
   4096 in the QLP QUALITY column removed; 49 ± 49 ppm if SPOC and TESS-SPOC cadences close
   in time to those flagged QLP cadences are also removed). The published 432 ± 56 ppm was
   not measured at a fixed ephemeris: it is the inverse-variance weighted mean of the 6
   series (out of 9 searched) that reached SNR >= 2, each taken at its own best time
   shift (column `dt_min` in `timing_per_series.csv`). An independent re-analysis by a
   separate AI agent, which did not use this project's code, finds 36 to 63 ppm for all 8
   series pooled across 7 detrending and quality-mask configurations, with a
   random-ephemeris null width of 62 to 71 ppm, so a 432 ppm transit at this ephemeris
   would be about 6 to 7 times the null width. If the transit times really scattered by
   up to about 80 min, as the note reported, a fixed-ephemeris depth would be diluted; finding 1
   shows that allowing per-series time shifts does not produce a signal above noise
   either.
3. **The strongest series rested on flagged cadences.** In the QLP sector 95 light curve,
   40 of the 99 cadences within ±0.5 transit durations (±33 min) of the note's sector-95
   transit times carry flag 2^30 (bit 31). This is QLP's own flag for poor-quality
   cadences, used since the extended mission EM2 (sectors 56 onward) according to the QLP
   documentation (tess.mit.edu/qlp, "Changes to MAST QLP light curves for EM2",
   2023-11-21); lightkurve's default quality mask does not remove it. Re-measured at that
   series' own transit times, the sector-95 depth is 516 ± 224 ppm with no flags removed
   and 17 ± 176 ppm with cadences flagged 2^30, 2048 or 4096 removed (columns
   `depth_bits0` and `depth_new`; flag 4096 does not occur in this sector; this is a different
   estimator from the note's: fixed transit times instead of an offset search, the median
   within ±0.4 transit durations, a baseline at 1 to 3 transit durations, and an error
   multiplied by 1.25, hence the larger error bars). Adding only this mask to our re-implementation of the note's
   procedure (variant `A_old` to `D_sym_mask`), the sector-95 value drops from SNR 4.3 at
   -38 min to SNR 1.0 at -20 min. The note's data package (`timing_per_series.csv`, row
   QLP-s95) gives 467 ± 108 ppm for this series.
4. **A blind box search does not single out this period.** The original report relied on
   a stitched SNR of about 8. In the independent re-analysis, in three settings (one
   series per sector; QLP only with a strict mask and the agent's own detrending; QLP only
   with lightkurve's default mask and QLP's detrended flux), the highest box-search peak
   within ±0.003 d of the claimed period, at any epoch, reaches SNR 5.2 to 5.5, while 0.2
   to 2.2% of all trial periods score higher and the strongest dip reaches SNR 6.6 to 8.7.

The out-of-sample alias test in sector 106 (item 6 of the original report) does not hold up
either. That test compared five candidate periods; the one the note favoured is
P = 4.814564 d (from `alias_test_s106.csv`, slightly different from the published
4.81460 d). The note reported 368 ± 144 ppm (SNR 2.55) for it, with SNR = depth / error. In
the independent re-analysis this period reaches an empirical z (depth divided by the
spread of depths at shifted phases) of 0.27 to 1.88 across 7 variants: combinations of
3x3 and 5x5 pixel apertures centred on the brightest pixel within the central 5x5 of the
cut-out, a 3x3 aperture shifted by one pixel to the cut-out centre, running-median
detrending windows of 0.75 or 1.5 d, and a strict or a looser quality mask (the looser
one keeps flag 2048, among others). The spread at shifted phases is 1.4 to 1.7 times
the formal error in that re-analysis, which is why an SNR based on formal errors looks
more significant than the data support.
On noise alone, the chance that at least one of five independent trials does as well as
P = 4.814564 d, computed as 1 - (1 - p)^5, is 14 to 91%, and in 3 of the 7 variants the
best of the five periods is not 4.814564 d but 4.789319 d.

The file `depth_by_pipeline.csv` in this folder was written by the project's vetting script
at the BLS peak period 4.878611 d, which was later rejected as an alias, not at the
published period (see `recheck_2026-11-09/vetting_run_ephemeris.txt`).

Supporting files are in `recheck_2026-11-09/`; its README says which rows and columns
each number above comes from and how the files were produced.

The original vetting report is kept below unchanged, except that its heading is relabelled
and three file names that had been corrupted by control characters are repaired. Some PNG
files it lists were never included in this public folder.

---

# Original report (superseded): TIC 253929649 -- transiting super-Earth candidate

> Vetting data package accompanying the RNAAS note submitted in October 2026
> (Godala 2026, *Research Notes of the AAS*, in review). Candidate parameters:
> P = 4.81460 +/- 0.00013 d, depth 432 +/- 56 ppm, Rp ~ 1.6 R_Earth.
> Files: figure1.png (note figure), alias_test_s106.csv (out-of-sample alias
> test), timing_per_series.csv (per-series mid-transit offsets),
> depth_by_pipeline.csv (cross-pipeline depths). Full vetting report below.

---

# TIC 253929649 — transiting super-Earth candidate (vetting report)

**Prepared 2026-10-12 by Mateusz Godala (Exo-Spectra project). Contact:
mt.godala@gmail.com. All analysis from public data (TESS QLP / TESS-SPOC /
SPOC light curves, TESScut FFIs, Gaia DR3) on a home PC.**

## Candidate parameters

| quantity | value |
|---|---|
| Period | **4.81460 ± 0.00013 d** (per-sector transit timing, 8-yr baseline; formal linear-fit error 2e-5 d inflated by sqrt(chi2_red)=6.3 — see Timing) |
| T0 (BTJD) | 1329.651 ± 0.010 (formal; ±0.06 with the same inflation) |
| Depth | **432 ± 56 ppm** (inverse-variance mean over 6 independent series) |
| Duration (T14) | ≈1.1 h, weakly constrained (BLS and trapezoid fits at the final ephemeris span 1.0–1.3 h) |
| T23/T14 (shape) | not constrained — trapezoid fits degenerate at per-sector SNR ~3 (see correction in Section "Vetting", item 5) |
| Rp | **1.64 R_Earth** (for TIC R* = 0.722 R_Sun) |
| Host | Tmag 10.68, Teff 4183 K, R* 0.72 R_Sun, M* 0.65 M_Sun, d = 70 pc |

## Detection history

Found in an independent BLS survey of 12,000 bright K/M dwarfs without
TOIs/CTOIs (QLP FFI light curves, 30-min binning, period grid 0.9–15 d).
The star shows the dip in **sectors 1, 28, 68, 95 (QLP) — four epochs
spanning 2018–2025** — with cross-sector period agreement at Tier-1 stage.

## Vetting performed

1. **Photometric tests (stitched QLP, noise-weighted):** odd/even depth
   difference z = 1.0; no secondary eclipse (z = 0.7 at phase 0.5);
   per-sector depth consistency chi2_red = 2.2.
2. **Contamination:** TIC contamination ratio 0.013. Gaia DR3: target
   RUWE = 0.88 (no astrometric binarity signal), parallax 14.24 mas
   consistent with TIC distance; nearest neighbour G = 18.4 at 22.6"
   (7.0 mag fainter); no comparable-brightness source within the aperture.
   (Gaia cannot exclude companions inside ~1".)
3. **Difference-image centroids** (TESScut 11x11, out-minus-in): on-target
   in both tested sectors (offsets 0.65 and 0.68 px; our discriminator
   threshold 1.0 px = 21").
4. **Cross-pipeline depth check:** the dip is present at consistent depth
   in SPOC (2-min), TESS-SPOC and QLP reductions and in our own aperture
   photometry of the sector-106 FFIs (depths 363–520 ppm across 6 series).
5. **Transit shape — CORRECTION (2026-10-12):** an earlier version of this
   report quoted a flat-bottomed trapezoid (517 ppm, T14 = 1.10 h,
   T23/T14 = 0.79) from the SPOC 2-min fold. That fold used the
   later-rejected alias period 4.8786 d; refolding at the final ephemeris
   shows the shape is NOT meaningfully constrained at per-sector SNR ~3
   (trapezoid fits are degenerate, with depth estimates 180–1100 ppm
   depending on method). The robust depth is the median in/out statistic
   per series (363–520 ppm; weighted mean 432 ± 56 ppm). The 2-min
   odd/even values previously quoted here (634 vs 577 ppm) were computed
   at the rejected alias and are withdrawn; the stitched odd/even test
   (z = 1.0, item 1) stands.
6. **Out-of-sample period confirmation:** sector 106 (observed
   May–June 2026) has no HLSP light curves yet; from our own TESScut
   aperture photometry, among 5 alias periods (pre-registered from the
   2018–2025 joint fit) **only P = 4.8146 d predicts the sector-106
   transits** (+368 ± 144 ppm in the predicted windows; the other aliases
   give 0 or negative depth). Prediction spans ~600 cycles from a
   2018-anchored epoch.

## Timing

Per-(pipeline, sector) local mid-transit offsets against the reference
ephemeris, searched within ±2 h, error floor T14/(2 SNR); 6 series yield
detections (SNR 2.1–4.3). Weighted linear fit:
P = 4.814602 ± 0.000020 d, T0(BTJD) = 1329.6512 ± 0.0098, but
**chi2_red = 39.5** — the O−C scatter (±80 min) exceeds white-noise
errors. We attribute this to FFI detrending systematics at per-sector
SNR ~ 3 and therefore quote the inflated errors above; the present data
do not justify a TTV claim.

## Known caveats (stated plainly)

* Per-sector SNR is ~3: this candidate lives at the detection floor of
  FFI photometry for a Tmag 10.7 star; the global stitched SNR is ~8.
* **Duration tension:** T14 ≈ 1.1 h at P = 4.81 d around a 0.72 R_Sun star
  implies a near-grazing geometry (b ≈ 0.9) for a circular orbit; the
  transit shape is too poorly constrained at this SNR to test this, and a
  smaller true stellar radius, eccentricity, or duration bias could
  ease it.
* O−C scatter above formal errors (see Timing).
* No radial-velocity or high-resolution-imaging constraints; a <1"
  companion is not excluded.

## Files

`spoc_fold_trapez.png` (SPOC 2-min fold + trapezoid), `timing_oc.png`
(O−C), `alias_test_s106.{csv,png}` (out-of-sample alias test),
`depth_by_pipeline.csv`, `timing_per_series.csv`,
`centroid_tic253929649_s{1,28}.png` (difference-image centroids).
