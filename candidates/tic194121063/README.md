# TIC 194121063 -- transiting sub-Neptune candidate

> Vetting data package accompanying the RNAAS note submitted in November 2026
> (Godala 2026, *Research Notes of the AAS*, in review, AAS83333). Candidate
> parameters: P = 4.732307 +/- 0.000047 d, depth 1504 +/- 104 ppm (QLP),
> Rp ~ 2.6-3.1 R_Earth. Files: `figure1.png` (note figure),
> `depth_by_pipeline.csv` (cross-pipeline depths at the final ephemeris),
> `ffi_depth.csv` (per-event depths from independent aperture photometry of
> the raw TESScut images), `centroid_neighbor.csv` (difference-image centroid
> projections onto the target--neighbour axis, per sector).
> Full vetting report below.

---

# TIC 194121063 — transiting sub-Neptune candidate (vetting report)

**Prepared 2026-11-07 by Mateusz Godala (Exo-Spectra project). Contact:
mt.godala@gmail.com. All analysis from public data (TESS QLP / TESS-SPOC
light curves, TESScut FFIs, Gaia DR3) on a home PC.**

## Candidate parameters

| quantity | value |
|---|---|
| Period | **4.732307 ± 0.000047 d** (per-series transit timing, 8 series 2019–2024, chi2_red(O−C) = 3.1, residuals within ±25 min) |
| T0 (BTJD) | 1714.9131 ± 0.0111 |
| Depth | **1504 ± 104 ppm** (QLP weighted, chi2_nu = 2.0 across sectors); raw-FFI joint 1637 ± 149 ppm; trapezoid floor 2190 ± 230 ppm (window-dilution spread) |
| Duration (T14) | 2.7 ± 0.2 h (trapezoid, bootstrap over events); central-transit prediction ≈2.2 h (excess consistent with mixed 200–1800 s FFI cadence smearing) |
| T23/T14 (shape) | **0.40 ± 0.17 — flat-bottomed** |
| Rp | **~2.6–3.1 R_Earth** (for TIC R* = 0.605 R_Sun) |
| Host | Tmag 11.16, Teff 4220 K, R* 0.61 R_Sun, M* 0.66 M_Sun, d = 75 pc; a ≈ 0.048 au, ~45 S_Earth |

## Detection history

Found in an independent BLS survey of bright K/M dwarfs without TOIs/CTOIs
(QLP FFI light curves). The star has no 2-min SPOC photometry; the signal is
recovered in **all 8 available (pipeline, sector) series** — TESS-SPOC
(sectors 15, 41, 55) and QLP (15, 41, 55, 75, 82) — spanning 2019–2024.

## Vetting performed

1. **Independent raw-pixel photometry:** because the two pipelines disagreed
   in sector 15, we performed aperture photometry of the raw 11×11 TESScut
   cutouts with no pipeline detrending: the dip is present independently in
   all five sectors (joint 1637 ± 149 ppm), significance against an empirical
   null of 300 random-phase folds z ≈ +5, and the raw images side with QLP in
   sector 15 — the TESS-SPOC PDCSAP detrending suppressed the signal there.
   Per-event depths: `ffi_depth.csv`.
2. **Eclipsing-binary tests:** odd/even depth difference 1.7σ; no secondary
   eclipse at phase 0.5; per-sector depth consistency chi2_nu = 2.0. The
   folded raw-FFI phase profile shows no sinusoidal (rotational) structure
   at P.
3. **Transit shape:** trapezoid fit to the QLP fold (bootstrap over events):
   floor 2190 ± 230 ppm, T14 = 2.7 ± 0.2 h, T23/T14 = 0.40 ± 0.17 —
   flat-bottomed.
4. **Multiplicity:** Gaia DR3 RUWE = 1.02, parallax 13.27 mas (consistent
   with the TIC distance), no non-single-star solution.
5. **Crowding — THE PRINCIPAL CAVEAT:** the field is crowded (TIC
   contamination ratio 0.33). The nearest significant neighbour, ΔG = 3.1 at
   9.8" (0.47 px), is a Gaia DR3 SOLAR_LIKE photometric variable with a 2.9%
   peak-to-peak G range; attributing the signal to it would require ~2.2%
   dips — comparable to its total variability — so Gaia photometry alone
   cannot exclude the blend. Difference-image centroids across all five
   sectors, projected onto the target–neighbour axis (0 = target,
   1 = neighbour) with bootstrap uncertainties, give **+0.19 ± 0.45**:
   consistent with the target, but rejecting the neighbour at only 1.8σ,
   since the 0.47-px separation is at the resolution limit of TESS
   difference imaging. Per-sector projections: `centroid_neighbor.csv`.

## What would settle it

A single ground-based transit observation (~1500–2200 ppm, T14 = 2.7 h,
Tmag = 11.2) would both test the ephemeris and trivially resolve the 9.8"
pair, removing the blend ambiguity. The next TESS coverage of this field is
expected in sectors 120–121 (2027), enabling a pre-registered out-of-sample
ephemeris test.
