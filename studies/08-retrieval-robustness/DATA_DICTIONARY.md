# Study 8 — data dictionary

## results/sigma_all.csv (116 rows)

One row per non-flat retrieval cell: 108 newly-run fits (144 runs minus 36
flat baselines, which enter only as the lnB reference) plus 8 L 98-59 b
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
