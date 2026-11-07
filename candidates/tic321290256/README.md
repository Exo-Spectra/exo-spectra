# TIC 321290256 -- transiting super-Earth-sized candidate

> Vetting data package accompanying the RNAAS note submitted in November 2026
> (Godala 2026, *Research Notes of the AAS*, in review). Candidate parameters:
> P = 2.331884 +/- 0.000013 d, depth 413-608 ppm (method-dependent),
> Rp ~ 1.8-2.2 R_Earth. Files: `figure1.png` (note figure),
> `timing_per_series.csv` (per-series mid-transit times and O-C),
> `depth_by_pipeline.csv` (cross-pipeline depths at the final ephemeris),
> `ffi_depth.csv` (51 per-event depths from independent aperture photometry
> of the raw TESScut images), `transit_windows.csv` (49 predicted transit
> windows, Nov 2026 - Feb 2027, for ground-based follow-up).
> Full vetting report below.

---

# TIC 321290256 — transiting super-Earth-sized candidate (vetting report)

**Prepared 2026-11-07 by Mateusz Godala (Exo-Spectra project). Contact:
mt.godala@gmail.com. All analysis from public data (TESS QLP / TESS-SPOC
light curves, TESScut FFIs, Gaia DR3) on a home PC.**

## Candidate parameters

| quantity | value |
|---|---|
| Period | **2.331884 ± 0.000013 d** (per-series transit timing on 5 QLP series, 2019–2026 baseline; uncertainties scaled by sqrt(chi2_red) = 3.1) |
| T0 (BTJD) | 1545.2201 ± 0.0091 |
| Depth | **413 ± 20 ppm** (QLP weighted); 511 ± 26 ppm (raw-FFI joint); 608 ± 38 ppm (trapezoid floor) — the spread is a window-dilution effect |
| Duration (T14) | 1.90 ± 0.07 h (trapezoid, bootstrap over events); central-transit prediction 2.16 h |
| T23/T14 (shape) | **0.50 ± 0.10 — flat-bottomed** |
| Rp | **~1.8–2.2 R_Earth** (for TIC R* = 0.82 R_Sun; range tracks the depth-method spread) |
| Host | Tmag 9.91, Teff 4787 K, R* 0.82 R_Sun, M* 0.77 M_Sun, d = 76 pc (Gaia parallax 13.14 mas); a ≈ 0.032 au, ~320 S_Earth |

## Detection history

Found in an independent BLS survey of ~30,000 bright K/M dwarfs without
TOIs/CTOIs (QLP FFI light curves). The dip is recovered in **all five QLP
sectors (9, 35, 62, 89, 100), spanning 2019–2026**, with cross-sector period
agreement at the Tier-1 stage. No SPOC 2-min photometry exists for this star.

## Vetting performed

1. **Independent raw-pixel photometry:** TESS-SPOC PDCSAP suppresses the
   signal (weighted −11 ± 75 ppm over its 3 sectors), so we performed our own
   aperture photometry of the raw 11×11 TESScut cutouts with no pipeline
   detrending: the dip is present in every sector (joint +511 ± 26 ppm),
   its significance against an empirical null of 300 random-phase folds is
   z = +5.6, and in sector 35 the raw images give +685 ± 59 ppm where
   TESS-SPOC shows −76 ± 85 ppm — a detrending artifact, resolved in favour
   of the signal. Per-event depths: `ffi_depth.csv`.
2. **Eclipsing-binary tests:** no secondary eclipse at phase 0.5. Odd/even
   depth difference 2.3σ on the QLP fold (360 ± 31 vs 465 ± 32 ppm) and 2.8σ
   per event (42 transits); the per-sector pattern is inconsistent (3 of 4
   sectors < 1σ; the aggregate is driven by sector 100), as expected for red
   noise rather than an eclipsing binary at 2P. Reported as a caveat to
   monitor.
3. **Transit shape:** trapezoid fit to the QLP fold (bootstrap over events):
   floor 608 ± 38 ppm, T14 = 1.90 ± 0.07 h, T23/T14 = 0.50 ± 0.10 —
   flat-bottomed. (An earlier internal fit gave T23/T14 = 0.26, an artifact
   of a mis-centred ephemeris dragged by the excluded TESS-SPOC s62 point;
   see `timing_per_series.csv`.)
4. **Multiplicity:** Gaia DR3 RUWE = 1.17, parallax consistent with the TIC
   distance, no non-single-star solution.
5. **Crowding — THE PRINCIPAL CAVEAT:** the nearest significant neighbour
   sits at 2.1" (ΔG = 2.9, 0.10 TESS px). Difference-image centroids across
   the five sectors show no consistent offset (which constrains more distant
   blends) but **cannot test a 0.10-px pair**. Attributing the signal to the
   neighbour would require ~0.7% eclipses on it, and Gaia DR3 carries no
   photometric-variability information for that star (two-parameter
   astrometry, no BP/RP photometry — its windows are blended with the 7×
   brighter target). **The blend scenario is untested from space data.**

## What would settle it

A single resolved ground-based transit observation — the pair separates
under sub-arcsecond seeing ($T_{14}$ = 1.9 h, ~500 ppm on a bright
Tmag = 9.9 star) — or speckle imaging plus an on/off-target photometric
series. Predicted windows for Nov 2026 – Feb 2027: `transit_windows.csv`
(mid-transit BTJD and UTC, computed from the ephemeris above).
