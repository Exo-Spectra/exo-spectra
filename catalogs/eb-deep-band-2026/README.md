# 1360 Uncatalogued Eclipsing Binary Candidates from TESS Full-Frame Images (the "deep band")

**EXO-SPECTRA project — catalog release, 2026-11-08**

## What this is

A catalog of **1360 eclipsing binary (EB) candidates around bright K/M dwarfs
(Tmag < 12), none of which appear in VSX, the ASAS-SN variable star catalog,
or the TESS Eclipsing Binary Catalog** (positional cross-match, 21 arcsec =
1 TESS pixel). They were found as a by-product of a transit-search re-scan of
30 000 TIC targets in TESS FFI photometry (QLP light curves), after fixing a
sigma-clipping step that had been silently removing deep eclipse signals
("the deep band": depths above ~2-5 ppt, where symmetric outlier clipping
deletes the eclipses themselves).

Why they are likely new: this depth-and-period regime is a blind spot of the
major variability catalogs. Gaia DR3 samples each star only ~70 times in
34 months and classifies these stars as spotted rotators (we verified this:
0 of our cross-matched targets appear in `vari_eclipsing_binary`, while 23%
appear in `vari_summary` as SOLAR_LIKE/RS); the TESS EB Catalog is built
primarily from 2-min cadence targets (only 12 of our objects match); ground
surveys (VSX 31%, ASAS-SN 4%) catch the shallower/brighter part.

## Method (short)

1. BLS transit search on QLP FFI light curves (asymmetric outlier clipping:
   upper 4 sigma only), multi-sector period clustering.
2. EB identification via odd/even eclipse-depth alternation
   (`EB_ODD_EVEN` verdict: alternating deep/shallow events means the BLS
   period is half the true orbital period).
3. **Per-target verification** (this catalog includes only targets that
   passed): local depth measured for every individual eclipse event across
   all sectors; weighted odd vs even depths; quality cut = at least 6
   measured eclipses and event-level SNR >= 5. Median: 13 eclipses,
   SNR 34 per target.

## Columns of `eb_catalog.csv`

| column | meaning |
|---|---|
| `tic` | TIC identifier |
| `ra_deg`, `dec_deg` | position (TIC, J2000, degrees) |
| `tmag`, `teff_k` | TESS magnitude and effective temperature (TIC) |
| `period_bls_d` | period found by BLS [d] |
| `period_adopted_d` | adopted orbital period [d]: 2x BLS period when the odd/even test is significant (abs(z) >= 3; 1310/1360 targets), else the BLS period |
| `t0_btjd` | reference eclipse epoch (BTJD = BJD - 2457000) |
| `depth_primary_ppm`, `depth_secondary_ppm` | weighted depths of the deeper/shallower eclipse family [ppm] |
| `oddeven_z` | significance of the primary-secondary depth difference |
| `n_eclipses_measured`, `n_sectors` | number of individually measured eclipse events / QLP sectors used |
| `snr` | combined event-level detection SNR |

## Caveats (read before using)

- **Depths are from QLP light curves and are NOT corrected for dilution**
  (TESS pixels are 21 arcsec; QLP applies its own crowding correction to the
  target star). For some targets the eclipsing source may be a nearby blended
  star; no per-target centroid vetting was performed at catalog scale.
- `period_adopted_d` doubles the BLS period based on the odd/even test; for
  the 50 targets with abs(z) < 3 the true period may still be twice the
  listed value (near-equal eclipse depths).
- Ephemerides come from a BLS fit to sectors spanning up to ~5 years;
  accumulated phase drift of minutes-to-tens-of-minutes is possible.
- QLP light curves in recent sectors fill the orbital-gap (downlink) window
  with scattered-light cadences flagged by quality bit 2^30, which the
  default lightkurve mask does not reject. This can fabricate single
  transit-like events, but cannot fabricate the alternating, periodic
  eclipse pattern required here: 99.9% of catalog entries are confirmed in
  >=2 sectors (median 4), none of the entries has its reference epoch inside
  such a window, and the per-event quality cut (>=6 eclipses, SNR>=5)
  removes targets dominated by single-window artifacts.
- Positional novelty was assessed on 2026-11-08 against VSX, ASAS-SN
  (II/366, v2021-10-08 via the AVSD snapshot of 2024-09-10) and the TESS EB
  Catalog (Prsa et al. 2022, J/ApJS/258/16); catalogs grow, re-check before
  claiming individual novelty.

## Files

- `eb_catalog.csv` — the catalog (1360 rows).
- `figure_period_depth.png` — adopted period vs primary depth, colored by Teff.
- `eb_lightcurves_png.zip` (GitHub release asset of tag `eb-catalog-v1.0`,
  72 MB) — a phase-folded light curve plot for every catalog entry
  (evidence for database submissions).

## License and citation

Data: see repository `LICENSE-DATA`. Please cite this repository release
(Zenodo DOI of tag `eb-catalog-v1.0`) and the EXO-SPECTRA project.

## Provenance

Produced with TESS QLP light curves (MAST), the TESS Input Catalog, VSX
(AAVSO), ASAS-SN and the TESS EB Catalog via the VizieR mirror at CfA
(Harvard). This research made use of lightkurve, astropy, pandas and the
NASA/MAST archive. Analysis pipeline co-developed with AI assistance
(Claude, Anthropic); all numbers in the catalog come from deterministic
scripts (BLS search, per-event depth measurement, positional cross-match)
run on public QLP light curves; the methodology is fully described above
and the scripts can be shared on request.
