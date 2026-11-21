# dasch_host_coverage.csv — century-long DASCH DR7 coverage of exoplanet host stars

Released 2026-11 as part of the final exo-spectra release. Built by the project's own scripts (not included in
this repository); the selection, cuts and columns are described below in enough detail to rebuild the table.

Sample: all hosts in the NASA Exoplanet Archive `pscomppars` table (15.11.2026) with B ≤ 14; 2246 rows with a result.
Photometry: DASCH DR7 via `daschlab` 1.1.0 (version DOI 10.5281/zenodo.23574176), B-like plate magnitudes calibrated
to APASS (refcat `apass`).

| column | meaning |
|---|---|
| hostname | host name as in the NASA Exoplanet Archive |
| ra_deg, dec_deg | archive position (deg, ICRS at the archive catalog epoch; NOT propagated to plate epochs — for a handful of high-proper-motion hosts the plate-era position differs, so their match can be wrong) |
| bmag_nasa | B magnitude from the archive (`sy_bmag`) |
| dasch_refcat_mag | `stdmag` of the matched DASCH reference-catalog source (refcat `apass`; for a minority of sources the catalog magnitude does not come from APASS) |
| dasch_match_sep_arcsec | separation of that source from the host position |
| ident_ok | False when \|dasch_refcat_mag − bmag_nasa\| > 1 mag or separation > 10″ (match not trusted) |
| n_detections | all DASCH detections of the source |
| n_std_quality | detections passing OUR approximation of the standard DASCH rejection set (flags HIGH_BACKGROUND, LARGE_ISO_RMS, LARGE_LOCAL_SMOOTH_RMS, CLOSE_TO_LIMITING, BIN_DRAD_UNKNOWN) — `daschlab`'s own `apply_standard_rejections` is described by its authors as "in development", so this is not an official product cut |
| n_good | detections passing our conservative cuts (a point is rejected if any of these is set: AFLAGS HIGH_BACKGROUND, BAD_PLATE_QUALITY, MULT_EXP_UNMATCHED, UNCERTAIN_DATE, MULT_EXP_BLEND, LARGE_ISO_RMS, LARGE_LOCAL_SMOOTH_RMS, CLOSE_TO_LIMITING, BIN_DRAD_UNKNOWN, CASE_B_BLEND, CASE_C_BLEND, CASE_BC_BLEND, LARGE_DRAD, PICKERING_WEDGE, SUSPECTED_DEFECT, SXT_BLEND, REJECTED_BLEND, LARGE_SMOOTHING_CORRECTION, TOO_BRIGHT, LOW_ALTITUDE; BFLAGS BLEND, SATURATED, NEAR_BOUNDARY, APERTURE_INCOMPLETE, ISOPHOT_INCOMPLETE, DEBLEND_OVERFLOW, EXTRACTION_OVERFLOW, PSF_SATURATED, TOO_BRIGHT, LARGE_JD_UNCERT; a non-zero `local_bin_reject_flag`; or a non-zero `reject`); used for everything below |
| year_first_good, year_last_good | time span of good points |
| median_mag_good, rms_mag_good | median and robust scatter (1.4826 · MAD, mag) of good points — not a plain r.m.s. |
| (note) | year_first_good, year_last_good, median_mag_good and rms_mag_good are left empty when n_good < 10 (117 hosts have 1–9 good points) |
| rms_mag_std | robust scatter (1.4826 · MAD, mag) of standard-quality points |
| n_comparison_stars | comparison stars of similar brightness (\|ΔB\| < 0.4 mag, widened to 1.0 then 2.0 mag when fewer than 5 found) in an annulus 3′–60′ from the host (actual range in this run: 3.0–19.7′); max 5 |
| evaluated | ident_ok and ≥ 2 comparison stars (the hosts that entered an internal population test, not published) |
| known_variable_5arcsec | VSX or ASAS-SN variable within 5″ (local copies, 8.11.2026); most matches are the planet hosts themselves (VSX type EP, transiting-planet host), and some matches have no type |
| known_variable_type | VSX/ASAS-SN type, if any |

Long-term trends are deliberately not provided: plate systematics reach ≲ 0.2 mag per century (Hippke et al. 2016).
About half of the hosts brighter than B = 7 (54 %) have n_good = 0 because of saturation flags; at B = 7–10 it is 9–13 %.

Required acknowledgment and citations (dasch.cfa.harvard.edu/citing/): NSF grants AST-0407380, AST-0909073,
AST-1313370; Grindlay et al. 2012 (2012IAUS..285...29G), Laycock et al. 2010 (2010AJ....140.1062L), Simcoe et al. 2006 (2006SPIE.6312E..17S), Tang et al. 2013
(2013PASP..125..857T, "Improved Photometry for the DASCH Pipeline" — verified in ADS 17.11.2026),
daschlab via its version DOI. Also to acknowledge: NASA Exoplanet Archive (`pscomppars`, queried 15.11.2026;
operated by Caltech/IPAC under contract with NASA), APASS (AAVSO Photometric All-Sky Survey, Henden et al.),
VSX (Watson et al. 2006) and ASAS-SN variable-star catalogs (local copies 8.11.2026) used for the
`known_variable_*` columns. Closest prior work: Tang, Sasselov & Grindlay 2013
(2013PASP..125..793T, "100-year DASCH Light Curves of Kepler Planet-Candidate Host Stars").
Data licence: DASCH does not state an explicit data licence; this table is a derived product shared with the
acknowledgement and citations above.

## What this table is and is not

It is a coverage table: for each exoplanet host with B ≤ 14 it says how many usable DASCH plate measurements exist,
over which years, and with what scatter. It is meant for anyone who wants to look at the century-long photometric
history of a host star.

It is not an event catalog. An earlier draft flagged short "dimming events"; an independent check showed that these
were single outlying plate measurements, and that the population test built on them depended on the choice of bins.
The `event_flag` column was therefore removed.
