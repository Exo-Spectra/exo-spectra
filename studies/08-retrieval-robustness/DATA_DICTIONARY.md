# Study 8 — data dictionary

## results/sigma_all.csv (116 rows)

One row per non-flat retrieval cell: 108 newly-run fits (144 cells minus 36
flat cells, which enter only as the lnB reference) plus 8 L 98-59 b
rows imported from Study 4's evidences.

| column | meaning |
|---|---|
| `pair` | pair key (see `pairs.csv`): gj1132b, wasp12b, hd189733b, l9859b, toi5205b_pipe, toi5205b_epoch, k218b, wasp39b, lhs1140b |
| `planet` | planet name |
| `axis` | the one upstream choice the pair varies: reduction / pipeline / epoch / instrument / control_consistent |
| `side` | a or b (mapping to concrete spectra in `pairs.csv`) |
| `code` | taurex or platon |
| `model` | gas model vs the flat baseline: h2o, ch4, co2, so2, with `_cloud` = + grey cloud deck |
| `err_variant` | raw = published/propagated errors; infl = errors inflated by the measured epoch-to-epoch disagreement (K2-18 b only) |
| `lnB_vs_flat` | ln Bayes factor of the gas model against the flat model (same code, same spectrum) |
| `sigma` | Benneke & Seager (2013) mapping of lnB to n-sigma, **capped at 6.99** (solver saturates at p ~ 1e-12); use lnB for anything quantitative |
| `best_params` | maximum-likelihood parameter vector `[rp_re, T, logx]` (+ `logpc` for cloud models): planet radius in Earth radii, isothermal temperature in K, log10 mixing ratio, log10 cloud-top pressure in Pa |
| `post_mean`, `post_std` | posterior means and standard deviations, same order |

## results/sensitivity_pairs.csv (58 rows)

One row per (pair, code, model, err_variant): the side-a vs side-b
comparison.

| column | meaning |
|---|---|
| `sigma_a`, `sigma_b` | detection significance on each side (capped at 6.99) |
| `dsigma` | sigma_b − sigma_a: the cost of swapping the upstream choice |
| `z_rp_re`, `z_T`, `z_logx`, `z_logpc` | posterior shift between sides in units of the combined error: \|mean_a − mean_b\| / sqrt(std_a² + std_b²). The LHS 1140 b control sets the noise baseline (≤ 3.8). |

## results/codes_agreement.csv (58 rows)

One row per (pair, side, model, err_variant): TauREx vs PLATON on the
identical spectrum.

| column | meaning |
|---|---|
| `sigma_taurex`, `sigma_platon` | per-code significance (capped) |
| `dsigma_codes` | sigma_platon − sigma_taurex |
| `dlnB_codes` | lnB_platon − lnB_taurex. Large values on hot giants partly reflect continuum physics differences between model classes (CIA in PLATON, none in our TauREx wrapper) — see README R7 before interpreting. |

## results/pairs.csv

The pair definitions: planet, axis, source catalog, spectrum identifiers
and files for sides a/b, authors, wavelength range, point counts. Archive
`file_a`/`file_b` are NASA Exoplanet Archive Atmospheric Spectroscopy
`.tbl` files; `survey` rows point to Study 3's homogeneous re-reductions.

## results/delta_sigma.png, results/fit__<pair>.png

Overview heatmap of dsigma across the matrix, and per-pair best-fit
overlays (data + both codes' best models, both sides) for the 8 newly-run
pairs.

## results/evidences_all.csv (144 rows, v1.1)

One row per matrix cell (pair, side, code, model, err_variant), including
the 36 flat baselines. L 98-59 b is not included: its runs belong to Study 4.

| column | meaning |
|---|---|
| `run_id` | id of the nested-sampling run behind the cell. Cells with the same id are identical runs (same spectrum, prior and seed 42, hence the same result): flat baselines of both codes except LHS 1140 b a/b, and 9 of the 10 TOI-5205 b side-a cells shared by the pipeline and epoch pairs. Treat them as one run, not as independent samples. 120 distinct ids. |
| `spectrum` | `<pair>__<side>`; the input file is `input_spectra/<spectrum>.csv` |
| `pair`, `planet`, `axis`, `side`, `code`, `model`, `err_variant` | as in `sigma_all.csv` |
| `n_points` | number of spectral points fitted |
| `lnZ`, `lnZ_err` | dynesty log-evidence and its error (natural log). Includes the Gaussian normalisation of the likelihood, so lnZ values of different spectra are not comparable; differences within one spectrum and error variant are the ln Bayes factors. |
| `max_lnL` | highest log-likelihood found by the sampler |
| `ncall` | likelihood calls used |
| `nlive_used` | live points of the kept run: 400 nominal, 200 after a budget retry (55 of the 120 distinct runs; their lnZ is noisier, see `lnZ_err`); empty = not recorded (19 runs made before the column existed) |
| `capped` | 1 = the kept run hit the call budget before dlogz = 0.1 (2 cells: TauREx CH₄, HD 189733 b); 0 = ended on dlogz = 0.1; empty = not recorded (19 runs; all stopped below every budget in force, so they ended on dlogz). Budget history: none for the earliest runs, then 300 000 calls, then 100 000 from 2026-10-11; flat runs are never budgeted. |
| `best_params`, `post_mean`, `post_std` | as in `sigma_all.csv` (flat model: `[depth_ppm]`) |

Full posterior samples were not kept.

## results/input_spectra/<pair>__<side>.csv (16 files, v1.1)

The spectra exactly as passed to the likelihood: both sides of a pair are
cut to the pair's common wavelength range; archive values converted to ppm;
asymmetric errors symmetrised as the mean of the upper and lower error.

| column | meaning |
|---|---|
| `wave_um` | central wavelength [µm] |
| `dwave_um` | bin width [µm] as given by the source. The model is the boxcar mean of its native-grid samples in [λ − Δλ/2, λ + Δλ/2]; where `dwave_um` is empty (all of `gj1132b__a` and `wasp12b__b`: no bandwidth in the archive) or the bin holds no model sample, the model is interpolated at λ. |
| `depth_ppm` | transit depth [ppm] |
| `err_ppm` | 1σ uncertainty [ppm], raw variant |

## results/retrieval_setup.csv (18 rows, v1.1)

One row per spectrum and error variant.

| column | meaning |
|---|---|
| `spectrum`, `pair`, `planet`, `side`, `err_variant` | identifiers |
| `err_scale` | factor applied to `err_ppm` for this variant (1 for raw; sqrt(2.048883) for K2-18 b infl) |
| `source_file` | original file name (`pairs.csv`) |
| `n_points`, `wave_min_um`, `wave_max_um` | fitted data |
| `model_lam_lo_um`, `model_lam_hi_um` | wavelength range of the forward model (data range widened by a factor 1.05) |
| `mp_me`, `rp_re`, `rs_sun`, `teff_k`, `teq_k` | planet mass [M⊕] and radius [R⊕], stellar radius [R☉] and temperature [K], planet equilibrium temperature [K] (NASA Exoplanet Archive) |
| `models` | models fitted to this spectrum |
| `prior_flat` | uniform prior on the flat depth [ppm]: [min(depth) − 5 max(err), max(depth) + 5 max(err)] with this variant's errors |
| `prior_rp_re`, `prior_T_k`, `prior_logx`, `prior_logpc_pa` | uniform prior bounds for the gas models: radius [R⊕], isothermal temperature [K], log10 mixing ratio, log10 cloud-top pressure [Pa] (cloud models only) |

Model physics: isothermal atmosphere, one trace gas in an H₂/He fill
(He/H₂ = 0.17), optional grey cloud deck. TauREx 3.3.2 wrapper: absorption +
Rayleigh + grey cloud, **no collision-induced absorption (CIA)**, 50 layers
from 1e-1 to 1e6 Pa, cross-sections R = 15000 (H₂O POKAZATEL, CH₄ MM,
CO₂ UCL-4000, SO₂ ExoAmes). PLATON 5.4 with its default settings and opacity
data (includes CIA, see README R7). lnZ values are therefore comparable only
within one code.

Sampler: dynesty 3.1.0, nlive 400 (200 after a budget retry), seed 42,
dlogz = 0.1, default bounding, `sample="rwalk"` for all runs except 4 flat
runs (run_id 1, 6, 62, 70) made earlier with dynesty's automatic choice
(`unif`). Code: `src/gap7_retrieval.py`.
