# Data dictionary — Study 4

Spectrum variant names used throughout: `pub_firefly`, `pub_eureka` (the two
published reductions of Bello-Arufe et al. 2025), `own_v1`…`own_v4` (our
Study 3 re-extraction of the 4 public G395H visits, in visit order
jw03942001001…jw03942004001), `own_avg` (weighted mean of the 4 own visits),
`own_avg_infl` (the same mean with errors inflated by the measured
visit-to-visit inconsistency), `synth_flat` (flat synthetic negative control).

## `results/spectra/*.csv` (one row per wavelength bin)

| column | meaning |
|---|---|
| `wave_um` | bin central wavelength [µm] |
| `dwave_um` | bin width [µm] |
| `depth_ppm` | transit depth [ppm] |
| `err_ppm` | 1σ uncertainty [ppm] |

## `results/evidences_taurex.csv`, `results/evidences_platon.csv` (one row per retrieval run, 27 rows each)

| column | meaning |
|---|---|
| `spectrum` | spectrum variant (see above) |
| `code` | retrieval code (`taurex` / `platon`) |
| `model` | atmosphere model: `flat` (constant depth), `so2`, `co2` |
| `n_points` | data points in the spectrum |
| `lnZ`, `lnZ_err` | ln evidence from nested sampling and its estimated error |
| `best_params` | maximum-likelihood sample: `[depth_ppm]` for `flat`; `[Rp (R⊕), T (K), log10 P_surf (Pa)]` for gas models |
| `max_lnL` | maximum log-likelihood found |
| `ncall` | number of likelihood evaluations |
| `runtime_s` | wall-clock runtime of the run [s] |

## `results/sigma_matrix.csv` (one row per model-vs-flat comparison, 36 rows)

| column | meaning |
|---|---|
| `spectrum`, `code`, `model` | as above (`flat` itself does not appear — it is the reference) |
| `lnB_vs_flat` | ln Bayes factor: lnZ(model) − lnZ(flat) for the same spectrum and code |
| `sigma` | lnB mapped to n-sigma via Benneke & Seager (2013); lnB ≤ 0 → 0.0 |
| `best_params` | maximum-likelihood sample of the gas model (as in the evidences files) |
