# Re-check of TIC 253929649 (2026-11-09): supporting files

These files support the withdrawal notice in `../README.md`. Notation: ± and +/- both mean
one standard error. "QLP flag 2^30" is bit 31 of the QUALITY column.

## Provenance

The three CSV files are the rows for TIC 253929649 extracted from tables in the project's
private working repository, produced on 2026-11-09 by the scripts in `scripts/` and
committed in commit `301432b`. The scripts import project modules that are not part of
this public repository (the QLP flag mask, and for two of them also the light-curve
preparation code), so they are included for reference, not as a stand-alone pipeline.

Which QLP mask each file used:

- `rnaas_timing_variants.csv` and `rnaas_recheck_summary.csv` were computed in the morning
  of 2026-11-09 with the QLP mask 2^30 | 2048 | 4096, before flag 2^29 was added to the
  project's default mask later the same day (in commit `301432b`). Flag 2^29 (bit 30) is
  not described by QLP, whose documentation lists bits 16 to 30 as unused. Re-running the scripts
  at that commit would also remove 2^29 cadences, which occur in QLP sector 95 (856
  cadences, none within the transit window; see `qlp_s95_flag_counts.txt`), so the
  out-of-transit baselines and random phases would change slightly. Variants `A_old` and
  `B_asym` use no mask and are not affected.
- `rnaas_mask4096.csv` sets its masks explicitly (see below) and is not affected.

## `rnaas_timing_variants.csv` (script `scripts/probe_notes_timing_variants.py`)

Reproduces the note's statistic: in each series, the transit-time offset within ±2 h that
maximises the depth SNR; a series counts as a detection at SNR >= 2. Offsets are measured
from the reference ephemeris used for the note's timing table, P = 4.8146 d,
T0 = 1329.63239 BTJD (about 27 min from the published linear ephemeris 4.814602 d,
1329.651 BTJD, which is used in `rnaas_recheck_summary.csv`).

- `variant`: `A_old` = preparation set up to match the note (symmetric 4-sigma clip, no
  QLP flag mask; a re-implementation inside the script, not the old code itself), which
  reproduces the note's per-series values within 0.1 in SNR and 1 min in offset.
  `B_asym` = clip of upward outliers only, no mask. `C_current` = upward clip plus QLP
  mask. `D_sym_mask` = symmetric clip plus QLP mask.
- `group`: which series are combined (`ALL` = the 8 archival series).
- `sum_snr`: sum of the best SNR over the series. `null_mean`, `null_std`: the same sum at
  400 random orbital phases (0.1 to 0.9 of the period away from the claimed transit).
  `z_emp`, `p_emp`: position of the real value in that null distribution.
- `n_det_snr2`: series with SNR >= 2; `null_median_n_det`: its median over random phases;
  `p_n_det`: fraction of random phases with at least as many detections.
- `per_series`: best SNR and the offset in minutes for each series.

The notice uses the row `A_old` / `ALL`, and for sector 95 also `D_sym_mask` / `ALL`
(same clip as `A_old`, only the QLP mask added).

## `rnaas_recheck_summary.csv` (script `scripts/probe_notes_recheck.py`)

Depth at the published linear ephemeris, inverse-variance weighted over series. Per series,
the depth is the median baseline flux (1.2 to 6 transit durations from mid-transit) minus
the mean
flux within ±0.6 transit durations; because this window is 20% wider than the transit, a
real box-shaped transit would read about 17% shallower.

- `variant`: `mask_qlp_only` = QLP cadences flagged 2^30, 2048 or 4096 in the QLP QUALITY
  column removed (2048 is the SPOC stray-light bit, carried into the QLP files);
  `mask_all_pipelines` = in addition, SPOC and TESS-SPOC cadences within 0.6 QLP exposure
  times of a flagged QLP time stamp are removed.
- `depth_w_ppm`, `err_w_ppm`: weighted depth and formal error. `null_depth_std_ppm`,
  `z_empirical_depth`, `p_null_depth`: against 300 random phases.

The notice uses the rows `mask_qlp_only` / `ALL` (61 ± 47 ppm) and
`mask_all_pipelines` / `ALL` (49 ± 49 ppm).

## `rnaas_mask4096.csv` (script `scripts/probe_notes_mask4096.py`)

QLP depth per sector under three explicit quality masks: `*_bits0` = no QLP flags removed,
`*_old` = 2^30 and 2048 removed, `*_new` = 2^30, 2048 and 4096 removed.

- `tc_source`: `linear` = the note's linear ephemeris; `T_obs` = that series' own transit
  times from the note (used for sectors 68 and 95).
- `n_in_raw`: raw cadences within ±0.5 transit durations (T14 = 1.1 h, so ±33 min) of the
  transit times, before any mask. `n_in_bit30`, `n_in_4096`: how many of them carry flag
  2^30 or 4096.
- `depth_*`, `err_*`: depth and error measured after light-curve preparation, with
  in-transit cadences defined more narrowly, within ±0.4 transit durations; `nin_*` is the
  number of such cadences. So `nin_*` is not a subset count of `n_in_raw`.

The file name refers to flag 4096, which matters for sectors 1 to 55; it does not occur in
sector 95. The notice uses the sector-95 row.

## Independent re-analysis: `independent_*.txt` and `scripts/independent/`

Outputs of an independent re-analysis done on 2026-11-09 by a separate AI agent (Claude).
It had the note's manuscript and this data package, but not the project's code, working
notes or conclusions, and it read the public light curves and TESScut images directly.
Its main scripts are in `scripts/independent/` (`core.py` holds loading, masks and
detrending; `tesscut_phot.py` its own aperture photometry).

- `independent_depths.txt` (`s04_depths.py`): per-series and pooled depths at the note's
  published linear ephemeris in 7 configurations. Quality mask "lkdef" = lightkurve
  default; "strict" = QUALITY == 0 in every series, which in practice changes only the QLP
  series. Detrending: running median of 0.4, 0.75 or 1.5 d, a spline, or the pipeline's
  own detrended flux followed by a 0.75 d running median. `phase-null SD` is the spread of
  the same statistic at random phases, `random-ephemeris null SD` at random ephemerides.
  The notice uses the `POOLED all series` lines.
- `independent_bls.txt` (`s06_bls.py`): blind box searches in three settings. The notice
  uses the lines starting `>> claimed P=4.814602` and `>> claimed-period peak`.
- `independent_s106_alias_test.txt` (`s11_s106_alias_authors.py`): the note's sector-106
  test of five candidate periods, repeated on the agent's own TESScut photometry and
  anchored at T0 = 1329.63239 BTJD (the epoch implied by the note's timing table). The
  alias-test period is P = 4.814564 d. Variants: 3x3 and 5x5 pixel apertures centred on
  the brightest pixel near the star, a 3x3 aperture shifted by one pixel to the cut-out
  centre, a strict mask (QUALITY == 0) or a looser one that removes only flag values 1, 2, 4,
  8, 16, 32, 128 and 16384 (and so keeps 64, 256, 512, 1024, 2048, 4096 and others), and running-median
  detrending windows ("runmed") of 0.75 or 1.5 d; the 7 variants are 7 combinations of
  these, not a full grid. `z_emp` is the depth divided by the
  spread of depths on a regular grid of shifted phases (steps of T14/4, excluding
  ±2 T14 around the predicted transit). "best of five ... reaches this" is the simulated
  chance that the best of five random-phase periods reaches the best observed value.
  "'any of five' p" is 1 - (1 - p_single)^5, the chance that at least one of five
  independent trials does as well as P = 4.814564 d; the notice quotes this one.

## Other files

- `note_timing_run.txt`: log of the note's own timing run (2026-10-11), showing that all
  8 archival series and the sector-106 photometry were searched and the three archival
  series with SNR < 2 were dropped before the fit and the 432 ± 56 ppm mean.
- `qlp_s95_flag_counts.txt`: number of cadences in the QLP sector-95 file, and counts of
  flag 2^29 (in total and within the transit window) and flag 2^30 (within the window).
  The file has 10626 cadences; `rnaas_mask4096.csv` lists 10614 (`n_cad`) because the
  script reads the file with lightkurve, which drops the 12 cadences that carry flags of
  lightkurve's default quality mask.

## Which period `depth_by_pipeline.csv` was computed at

- `vetting_run_ephemeris.txt`: excerpt from the log of the project's deep-vetting run of
  2026-10-11 that also wrote `depth_by_pipeline.csv` (both files carry the same
  timestamp). That run selected P = 4.878611 d, T0 = 1329.63239 BTJD as its ephemeris,
  and the per-series depths it logged are identical to `depth_by_pipeline.csv`. The
  period was later rejected as an alias in favour of 4.8146 d.
- `independent_extra_depth_by_pipeline.txt`: the independent agent's own per-series depths
  at that BLS ephemeris, next to the values in `depth_by_pipeline.csv`. In the "lkdef"
  column they agree within the errors for most series, but QLP sector 95 differs
  (327 ± 120 ppm against 621 ± 111 ppm).
