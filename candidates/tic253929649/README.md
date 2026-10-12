# TIC 253929649 -- transiting super-Earth candidate

> Vetting data package accompanying the RNAAS note submitted in October 2026
> (Godala 2026, *Research Notes of the AAS*, in review). Candidate parameters:
> P = 4.81460 +/- 0.00013 d, depth 432 +/- 56 ppm, Rp ~ 1.6 R_Earth.
> Files: igure1.png (note figure), lias_test_s106.csv (out-of-sample alias
> test), 	iming_per_series.csv (per-series mid-transit offsets),
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
