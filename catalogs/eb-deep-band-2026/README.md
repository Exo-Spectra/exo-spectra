# WITHDRAWN (v1.1, 2026-11-09): "1360 Uncatalogued Eclipsing Binary Candidates from TESS Full-Frame Images"

**EXO-SPECTRA project. Catalog release v1.1, which withdraws 1340 of the 1360 entries of v1.0
and withdraws the claim that the remaining ones are uncatalogued.**

## Withdrawal notice

I withdraw this catalog. An independent re-check on the night after the v1.0 release
(2026-11-08) showed that the evidence given in v1.0 does not hold. For about 98.5% of the
entries the listed period and epoch do not predict repeating eclipses. The large majority of
these entries are instrumental artifacts. A few are real eclipsing binaries for which v1.0
gave a wrong period. Only 20 entries pass the repetition test at the listed ephemeris, and
a second check showed that for 8 of these 20 a known VSX eclipsing binary with a matching
period lies nearby, so the eclipses most likely belong to that object.

Please do not use the v1.0 list (Zenodo DOI 10.5281/zenodo.24058988) and do not cite it as
a source of new eclipsing binaries.

### What went wrong

1. **A quality flag was not rejected.** QLP light curves mark cadences that the QLP
   pipeline itself considers bad ("Low precision points" in the QLP data product
   description). In sectors 1 to 55 this flag is bit 13 of the QUALITY column (value
   2^12 = 4096). From sector 56 on it is bit 31 (value 2^30). The change is described on the
   MIT QLP page (https://tess.mit.edu/qlp/, entry "Changes to MAST QLP light curves for EM2
   11-21-23"). My mask rejected 2^30 and 2048 but not 4096, and the default lightkurve mask
   does not reject either QLP flag.
2. **The outlier clipping was one-sided.** To keep deep eclipses, the method clipped only
   upward outliers. Flagged cadences with flux drops of up to about 50% therefore stayed in
   the light curves and entered the period search.
3. **The verification counted measured events, not significant ones.** The criterion
   "at least 6 measured eclipses and combined SNR >= 5" is passed by one huge drop plus five
   epochs with no signal. The odd/even test then reported a large depth difference, which I
   read as a primary/secondary pair.
4. **I did not test for shared event times.** The same bad cadences affect many stars, so
   the "eclipses" of different stars fall at the same times.
5. **The novelty check used too small a radius.** I matched against VSX within 21 arcsec
   (one TESS pixel). The QLP aperture is 2.5 to 3 pixels in radius and the wings of the
   point spread function of a bright star reach further, so an eclipsing binary tens of
   arcseconds away, in one case 117 arcsec, shows up in the light curve of the target.

### The numbers

For every entry I re-measured each predicted eclipse of the deeper eclipse family (adopted
period, epoch moved to the deeper eclipse), once with the quality mask used for v1.0 and
once with bit 4096 rejected as well. An event is significant when its depth is positive and
its SNR is above 5. An entry is coherent when it has at least 3 significant events, these
are at least half of the measured events, and their median depth is at least 0.3 of the
maximum.

| | v1.0 quality mask | bit 4096 rejected as well |
|---|---|---|
| entries with 0 significant events | 50 | 1227 |
| with exactly 1 | 851 | 86 |
| with exactly 2 | 362 | 10 |
| with 3 or more | 97 | 37 |
| coherent (retained) | 25 | 20 |
| predicted event times shared by 3 or more catalog stars | 165 | 1 |

- With bit 4096 rejected, 1178 entries lose every significant event they had.
- Event times are compared after rounding to 0.01 d. With the v1.0 mask the most frequent
  time (BTJD 2081.20) is shared by 76 catalog stars.
- 17 further entries keep 3 or more significant events but are not coherent at the
  catalog ephemeris.

This test only looks at the epochs predicted by the v1.0 period. It cannot find a real
eclipsing binary whose v1.0 period is wrong. A search without the v1.0 ephemeris found
at least 7 withdrawn entries with real, repeating, multi-cadence eclipses at a different
period (table below). There may be more.

### What is in v1.1

`eb_catalog.csv` still has 1360 rows, so that every v1.0 entry can be looked up. The 14
columns of v1.0 are unchanged. The new columns are:

| column | meaning |
|---|---|
| `status` | `withdrawn` (1340 rows) or `retained` (20 rows). The status records the result of the repetition test at the v1.0 ephemeris (the "coherent" row of the table above). `retained` does not mean that the object is new or that the eclipses are on this star. |
| `n_events_measured_masked` | predicted eclipses of the deeper family (v1.0 adopted period) that could be measured with bit 4096 rejected |
| `n_significant_events_masked` | how many of them are significant (SNR > 5, depth > 0) |
| `period_best_d`, `t0_best_btjd` | re-derived period and the epoch of an observed eclipse of the family that the fit locked on, normally the deeper one; filled for the retained rows and for the withdrawn rows of the second table |
| `ephemeris_source` | `bls_refit` (box least squares fit within 0.5% of the starting period, then a straight-line fit to the eclipse times of the individual sectors) or `v1.0_t0_shifted` (v1.0 adopted period with the v1.0 epoch, moved by one BLS period where that was needed to land on the deeper eclipse; used for two rows where the refit was not better) |
| `n_events_best`, `n_significant_best` | measured and significant eclipses at `period_best_d` |
| `depth_sector_min_ppm`, `depth_sector_max_ppm` | smallest and largest per-sector median eclipse depth [ppm]; retained rows only |
| `vsx_name`, `vsx_type`, `vsx_period_d`, `vsx_sep_arcsec` | VSX object within 150 arcsec whose period equals the BLS, adopted or re-derived period, or half or twice one of them, within 0.4% |
| `note` | short remarks |

### The 20 retained entries

All 20 show real, multi-cadence, repeating eclipses. That is all the word "retained" says.

| TIC | best period [d] | t0 [BTJD] | significant / measured eclipses | depth, shallowest to deepest sector [ppm] | VSX object within 150 arcsec with matching period | note |
|---|---|---|---|---|---|---|
| 96795030 | 3.41888 | 1470.043 | 36 / 36 | 59065 to 75707 | none |  |
| 118832092 | 1.85540 | 1518.742 | 65 / 65 | 53330 to 114583 | ASASSN-V J090003.48-224133.8 (EA, 3.7108 d, 39 arcsec) | known VSX eclipsing binary nearby with 2.0 x this period |
| 220052771 | 2.23329 | 1470.083 | 25 / 25 | 50558 to 67364 | none |  |
| 232378691 | 2.27181 | 2201.885 | 21 / 21 | 29768 to 47611 | ASASSN-V J061122.65+030407.9 (EA, 2.2719 d, 62 arcsec) | known VSX eclipsing binary nearby with the same period |
| 95247597 | 14.37689 | 1496.987 | 7 / 7 | 60416 to 65767 | none |  |
| 130441077 | 4.15814 | 1494.369 | 21 / 21 | 24594 to 35509 | ASAS J075945-2651.7 (ESD/ED, 2.079 d, 95 arcsec) | known VSX eclipsing binary nearby with 0.5 x this period |
| 12929185 | 4.74236 | 1685.258 | 31 / 31 | 11396 to 17852 | none |  |
| 342061072 | 10.30138 | 2011.224 | 11 / 12 | 47791 to 51879 | none |  |
| 141812312 | 10.84475 | 1902.463 | 5 / 5 | 49657 to 53818 | none |  |
| 124181500 | 8.36967 | 2202.018 | 9 / 11 | 326 to 18924 | ASASSN-V J060828.66-215840.6 (EA, 8.3688 d, 35 arcsec) | v1.0 period was a wrong doubling; eclipses nearly vanish in sector 98, so they are probably on a neighbouring star; known VSX eclipsing binary nearby with the same period |
| 367195056 | 6.14360 | 2503.993 | 13 / 13 | 14236 to 16047 | none |  |
| 104705505 | 2.17602 | 3776.253 | 17 / 17 | 53171 to 62854 | none |  |
| 272646087 | 1.19176 | 1571.694 | 88 / 148 | 18 to 21982 | FW Vel (EA/DM, 2.3835 d, 117 arcsec) | eclipses nearly vanish in sectors 36, 37, 64, 99, 100, so they are probably on a neighbouring star; known VSX eclipsing binary nearby with 2.0 x this period |
| 381729821 | 2.09540 | 1792.395 | 29 / 29 | 61396 to 63342 | ASASSN-V J024419.41+392030.0 (EA, 2.0955 d, 39 arcsec) | v1.0 ephemeris misses the eclipse by 2.6 h; known VSX eclipsing binary nearby with the same period |
| 395477762 | 6.37143 | 2771.417 | 11 / 11 | 4046 to 4940 | Gaia DR3 4314106430556684160 (E, 6.3713 d, 28 arcsec) | v1.0 ephemeris misses the eclipse by 4.6 h; known VSX eclipsing binary nearby with the same period |
| 19691949 | 4.08167 | 1494.703 | 20 / 24 | 434 to 5571 | none | eclipses nearly vanish in sector 99, so they are probably on a neighbouring star |
| 345503390 | 7.43503 | 1717.224 | 16 / 16 | 2001 to 2951 | none |  |
| 231888176 | 4.96382 | 1626.746 | 16 / 28 | 494 to 4194 | none | eclipses nearly vanish in sector 104, so they are probably on a neighbouring star |
| 368035954 | 6.31364 | 1716.570 | 10 / 13 | 1670 to 2643 | Gaia DR3 1858610096216321664 (E, 6.3139 d, 75 arcsec) | known VSX eclipsing binary nearby with the same period |
| 417754962 | 8.75107 | 1843.312 | 6 / 6 | 3718 to 5243 | none |  |

What the table shows:

- **8 of the 20 retained entries have a known VSX eclipsing binary 28 to 117 arcsec away
  with the same period, or half or twice it.** They are not new, and the eclipses most
  likely belong to the VSX object, not to the listed star. The match uses a local VSX
  snapshot downloaded from VizieR on 2026-11-08.
- The other 12 retained entries have no such VSX object within 150 arcsec. I have not checked
  them against Gaia DR3, ASAS-SN as a separate catalog, or the FFI extension of the TESS
  Eclipsing Binary Catalog that is in preparation, so I do not claim that they are new.
- For 4 retained entries the eclipse depth in one sector is less than 15% of the depth in
  another sector. Eclipses that come and go with the pointing of the spacecraft or with the
  photometry method are a sign that they happen on a neighbouring star. No centroid test
  was done for any entry.
- v1.0 had wrong ephemerides for three of them: for TIC 124181500 the adopted period was a
  wrong doubling, and for TIC 395477762 and TIC 381729821 the epoch misses the eclipse by
  hours. The hours quoted in the notes are counted from the v1.0 epoch after it was moved
  to the deeper eclipse, not from the raw `t0_btjd` column. For TIC 118832092, 130441077 and 272646087 the VSX period differs from mine by a
  factor of two, and my data do not decide which is right.
- The re-derived periods have an uncertainty of about 1e-4 d, and the epochs of about half
  an hour. They are good enough to find the eclipses, not for timing work. For
  TIC 130441077 the eclipses at phase 0 and at phase 0.5 are equally deep, so "the deeper
  family" is not defined there.
- **The v1.0 column `t0_btjd` must not be used.** It came from the fit at the undoubled
  period, and for 8 of the 20 retained rows it points to the shallower eclipse, about half a
  period away from `t0_best_btjd`. The other v1.0 columns of the retained rows (depths, SNR,
  number of eclipses) were measured without rejecting bit 4096 and can be biased.
- `n_significant_best` can be smaller than `n_significant_events_masked` for the same row,
  because the two measurements use different eclipse durations.
- The depths are QLP depths, not corrected for dilution.

### Withdrawn entries that are real eclipsing binaries with a wrong v1.0 period

| TIC | v1.0 adopted period [d] | best period [d] | t0 [BTJD] | significant / measured eclipses | VSX object within 150 arcsec with matching period |
|---|---|---|---|---|---|
| 307361219 | 21.994986 | 23.82827 | 2089.299 | 24 / 24 | none |
| 114300563 | 12.147610 | 23.12038 | 2243.068 | 3 / 3 | none |
| 16778117 | 11.107010 | 10.67843 | 1845.656 | 4 / 4 | none |
| 365699980 | 2.379286 | 14.27563 | 2836.285 | 10 / 10 | none |
| 231708358 | 14.183038 | 12.98385 | 1333.916 | 8 / 8 | Gaia DR3 6504137089668690560 (E, 6.4919 d, 53 arcsec) |
| 46629702 | 4.333707 | 4.33023 | 1791.830 | 32 / 34 | ASASSN-V J033557.81+273454.8 (EA, 4.3297 d, 57 arcsec) |
| 269971671 | 3.242931 | 3.24828 | 1493.992 | 33 / 33 | ASASSN-V J073802.11+104702.2 (EA, 3.2482 d, 74 arcsec) |

These rows stay `withdrawn`, because the status records the test at the v1.0 ephemeris and
they did not pass it. The list is a lower limit: it comes from one search without the v1.0
ephemeris by an independent reviewer, which I confirmed with my own measurement for the
rows above. For some of them only one of two eclipse families is counted here (TIC 16778117
has a second family of similar depth near phase 0.52).

### Further caveats

- The audit columns were measured with QUALITY bits 2^30, 2048 and 4096 rejected. The
  re-derived ephemerides also reject bit 30 (value 2^29), which appears in QLP files from
  sector 94 on and is not documented.
- QLP time stamps in sectors 99 to 104 have a known barycentric correction error. The
  corrected files have not yet been redelivered to MAST (https://archive.stsci.edu/hlsp/qlp).
  Several retained entries use these sectors.
- The light-curve plots in `eb_lightcurves_png.zip` (asset of the v1.0 GitHub release)
  show the artifacts for the withdrawn entries and should not be used as evidence.

### Where the claim was sent, and what I did

The v1.0 manuscript submitted to OEJV on 2026-11-08 was withdrawn on 2026-11-09, before
review. The offer of the sample to the TESS Eclipsing Binary Catalog team was corrected
the same day. Nothing was submitted to VSX or VizieR.

---

## Original description of v1.0 (kept for the record; its claims are withdrawn)

The text below is the v1.0 README, unchanged apart from heading levels. Its statements
about verification, robustness against scattered-light cadences and novelty do not hold.

### 1360 Uncatalogued Eclipsing Binary Candidates from TESS Full-Frame Images (the "deep band")

**EXO-SPECTRA project — catalog release, 2026-11-08**

#### What this is

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

#### Method (short)

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

#### Columns of `eb_catalog.csv`

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

#### Caveats (read before using)

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
  (II/366, v2021-10-08 via the AVSD snapshot of 2024-09-10), the TESS EB
  Catalog (Prsa et al. 2022, J/ApJS/258/16) and the TESS Ten Thousand
  Catalog (Kostov et al. 2025, J/ApJS/279/50: zero matches by TIC or
  within 21 arcsec -- their vetted sample is much fainter, median T=13.8);
  catalogs grow, re-check before claiming individual novelty. Note that
  Kostov et al. also published ~900k UNVETTED neural-network candidates;
  overlap with that list was not checked (it is not a vetted catalog).

#### Files

- `eb_catalog.csv` — the catalog (1360 rows).
- `figure_period_depth.png` — adopted period vs primary depth, colored by Teff.
- `eb_lightcurves_png.zip` (GitHub release asset of tag `eb-catalog-v1.0`,
  72 MB) — a phase-folded light curve plot for every catalog entry
  (evidence for database submissions).

#### License and citation

Data: see repository `LICENSE-DATA`. Please cite this repository release
(Zenodo DOI of tag `eb-catalog-v1.0`) and the EXO-SPECTRA project.

#### Provenance

Produced with TESS QLP light curves (MAST), the TESS Input Catalog, VSX
(AAVSO), ASAS-SN and the TESS EB Catalog via the VizieR mirror at CfA
(Harvard). This research made use of lightkurve, astropy, pandas and the
NASA/MAST archive. Analysis pipeline co-developed with AI assistance
(Claude, Anthropic); all numbers in the catalog come from deterministic
scripts (BLS search, per-event depth measurement, positional cross-match)
run on public QLP light curves; the methodology is fully described above
and the scripts can be shared on request.
