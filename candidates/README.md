# Planet candidates from the TESS search: final status (2026-11-21)

**None of the planet candidates from this project is confirmed, and none should be treated as a planet host on the
basis of this repository.** This page lists every candidate that was announced publicly or reached working or hold
status in our vetting, with the reason it was closed. Many further signals were examined and rejected during
vetting (eclipsing neighbours, stellar variability, detrending artifacts) and are not listed. The project ended in November 2026; the candidates are left here so that anyone with better data can
pick them up, and as a record of which false-positive mechanisms caught them.

| TIC | period (d) | final status | why |
|---|---|---|---|
| 253929649 | 4.8146 | withdrawn (2026-11-09) | on re-check the signal does not stand out from noise in the statistic its note relied on; see `tic253929649/` |
| 142825877 | 7.5211 | withdrawn (2026-11-09) | after the full QLP quality mask the signal is at about 2σ (p ≈ 0.05), not the S/N 8.9 stated in its submitted RNAAS note (AAS82681); a withdrawal request was sent to the journal on 2026-11-09 |
| 194121063 | 4.732307 | **withdrawn (2026-11-17): blend** | the dips come from a previously uncatalogued eccentric eclipsing binary, TIC 194121086, 39″ away (now in `catalogs/new-eclipsing-binaries-2026/`); see `tic194121063/` |
| 321290256 | 2.331884 | unconfirmed, closed | empirical false-alarm probability p ≈ 0.05–0.06 when the period is chosen on all data; a neighbour at 2.1″ cannot be excluded without high-resolution imaging; the RNAAS note was declined on 2026-11-18 under the journal's new policy for single-object transit notes; see `tic321290256/` |
| 92134828 | 5.308735 | unconfirmed, closed | an independent re-reduction of sector 107 gives p = 0.025–0.08 instead of 0.002; no second reduction of that sector exists |
| 217972650 | 4.362166 | unconfirmed, closed | an out-of-sample test on sectors 103/104 (+2.4σ) did not survive the QLP quality and camera common-mode masks described below (−0.4σ with the masks) |
| 85379511 | 9.474657 | likely false positive, closed | real, consistent signal in three pipelines and raw pixels, but V-shaped (grazing) and Gaia RUWE = 3.87: the textbook case of an eclipsing binary inside an unresolved double star; only speckle imaging or radial velocities can decide |
| 136497248 | 4.1567 | closed | depth grew from about 180 to 490 ppm between 2019 and 2024 and the signal vanishes in sector 99: either a precessing eclipsing binary or a QLP systematic; not followed up |
| 326481416 | 8.3442 | unconfirmed, closed | depth consistent in all 9 series from 3 pipelines (~330 ± 53 ppm), but the combined S/N of 6.2 is a formal value without a null calibration and is below our announcement threshold of 7; the centroid check that looked on target used an older method we later found unreliable, and the QLP part of this analysis predates the quality mask below |
| 305776955 | 6.4438 | closed | the QLP depth (~1500 ppm) is not reproduced elsewhere: raw-pixel photometry gives 257 ± 101 ppm and SPOC about zero; only a weak 2.9σ raw-pixel trace remains (frozen, not refuted) |
| 12129646 | ~7.9 | closed | only one sector of data (period from that sector only); never vetted across sectors |

## The lesson that mattered most

From about sector 56 onward, QLP fills the orbital-gap downlink windows with cadences affected by scattered light
(quality bit 2^30, often together with 2048), and in sectors 1–55 it flagged its own bad cadences with bit 4096. The
default `lightkurve` quality mask rejects none of these bits. In our search, keeping them produced false transit
detections; masking them (bits 2^30, 2048, 4096, the undocumented 2^29 from sector 94, and SPOC bits 8192 and 32768)
removed, among others, the withdrawn 1360-entry eclipsing-binary catalog and the TIC 217972650 out-of-sample signal.
If you search QLP light curves for shallow transits, mask these bits first.
