# Study 4 — How robust is the SO₂ detection on L 98-59 b? A 54-run verification matrix

*An independent verification of the reported ~3σ evidence for sulfur dioxide on
the sub-Earth L 98-59 b ([Bello-Arufe et al. 2025](https://doi.org/10.3847/2041-8213/adaf22)):
the same Bayesian model comparison, run on 9 spectrum variants × 3 atmosphere
models × 2 independent retrieval codes.*

## TL;DR (for everyone)

L 98-59 b is a rocky planet smaller than Earth. In 2025, a science team
reported a possible sign of sulfur dioxide (SO₂) in its atmosphere — a hint of
volcanic activity — at about 3σ confidence. That level means "interesting, but
not yet certain". We tested how stable this result is. We used every available
version of the planet's spectrum: the two published versions of the same
JWST observations, plus our own re-extraction of the 4 public transits from
Study 3. We then asked two independent retrieval codes the same question:
does an SO₂ atmosphere fit the data better than a flat line?
The answer depends on which version of the data you use. On the published
Eureka! reduction, we confirm ~3.1–3.2σ for SO₂. On the other published
reduction (FIREFLy) — **the same observations, different software** — the
signal drops to 1.3σ or less: no real evidence. In our own re-extraction, the
signal mostly comes from a single visit that we already marked as
inconsistent with the other three. The two retrieval codes agree everywhere,
so the instability comes from the data processing, not from the analysis
software. Our conclusion: the SO₂ detection is fragile. It may still be real,
but the current data cannot confirm it at the ~3σ level.

## Summary (for the technical reader)

Bello-Arufe et al. (2025) reported evidence for SO₂ on L 98-59 b from 4
JWST/NIRSpec G395H transits, at 2.2–3.6σ depending on the retrieval code, and
published two reductions of the dataset: FIREFLy (80 points) and Eureka!
(218 points; the reduction behind the headline result). We ran a uniform
verification matrix: 9 spectrum variants × {flat, SO₂, CO₂} × {TauREx 3,
PLATON 5.4} = 54 nested-sampling retrievals with identical data, priors and
sampler settings (dynesty, nlive = 400, fixed seed). The 9 variants: the two
published reductions, our own uniform re-extraction of the same 4 public
visits from [Study 3](../03-jwst-mini-survey/) (per visit, their weighted
mean, and the mean with errors inflated by the measured visit-to-visit
inconsistency), and a flat synthetic spectrum as a negative control.

Significance of each atmosphere model vs a flat line (Bayes factor mapped to
σ following Benneke & Seager 2013, as in the original paper):

| spectrum | SO₂ TauREx | SO₂ PLATON | CO₂ TauREx | CO₂ PLATON |
|---|---|---|---|---|
| published Eureka! (218 pts) | **3.19σ** | **3.10σ** | 0.0 | 0.0 |
| published FIREFLy (80 pts) | **1.33σ** | **0.0σ** | 0.0 | 0.0 |
| own mean of 4 visits | 3.20σ | 3.37σ | 0.0 | 0.0 |
| own mean, errors ×√χ² | 2.67σ | 2.77σ | 2.09σ | 2.20σ |
| own visit 1 | 2.07σ | 1.93σ | 0.0 | 0.0 |
| own visit 2 | 0.0σ | 0.0σ | 0.0 | 0.0 |
| own visit 3 | 1.82σ | 1.71σ | 2.06σ | 2.15σ |
| own visit 4 | **4.64σ** | **5.15σ** | 2.37σ | 2.58σ |
| flat synthetic (control) | 0.0σ | 0.0σ | 0.0σ | 0.0σ |

Four findings:

1. **Reproduction: yes.** On the Eureka! spectrum both codes independently
   reproduce the published detection (3.19σ / 3.10σ).
2. **The detection does not survive a change of reduction.** The same
   observations reduced with FIREFLy give 1.33σ / 0.0σ — the choice of
   reduction pipeline flips the scientific conclusion. This is the
   single-planet counterpart of the population-level result in Study 1
   (same-data, different-pipeline pairs).
3. **In our own re-extraction, the averaged data reproduce ~3σ, but the
   signal is fragile.** The 4 visits are mutually inconsistent (6 of 6 pairs
   discrepant in Study 3); propagating that inconsistency into the error bars
   cuts SO₂ to 2.67σ / 2.77σ while CO₂ rises to ~2.1–2.2σ — the data no longer
   single out SO₂. Per visit, the preference is driven largely by visit 4
   (4.64σ / 5.15σ) — exactly the visit Study 3 flagged as the epoch outlier.
4. **The two retrieval codes agree.** No cell of the matrix changes its
   conclusion between TauREx and PLATON. The instability comes from the data
   (reduction, epoch), not from the retrieval software.

**Conclusion:** the SO₂ detection on L 98-59 b is fragile — reproducible for
only one of the two published reductions, cut below 3σ by honest propagation
of the measured epoch-to-epoch inconsistency, and driven in large part by one
anomalous visit. We do not claim SO₂ is absent; we claim the current data
cannot establish it at the ~3σ level, and the reported significance is the
upper envelope over the choice of reduction.

## For the expert

**Data variants** (`src/l9859b_prepare.py`, all normalized to wave_um /
dwave_um / depth_ppm / err_ppm). Published: archive spectra of Bello-Arufe
et al. (2025), FIREFLy (80 pts) and Eureka! (218 pts). Own: the Study 3
re-extraction of the 4 public G395H visits (jw03942001001…jw03942004001,
41 bins each after the quality cut). `own_avg` = inverse-variance weighted
mean on the common bin grid; `own_avg_infl` = the same with per-bin errors
inflated by √max(1, χ²_red) of the 4-visit scatter — this carries the
measured epoch inconsistency into the averaged error bars (published joint
analyses do not). `synth_flat` = flat truth at the median depth + Gaussian
noise at the `own_avg` error bars (fixed seed).

**Models** (`src/l9859b_retrieval.py`), deliberately matched to the paper's
model-comparison setup rather than full physics: isothermal, well-mixed,
`flat` (1 parameter: constant depth), `so2` and `co2` (3 parameters:
Rp [0.70–1.00 R⊕], T [300–1000 K], log₁₀ P_surf [1–7, Pa]; uniform priors;
99.9% of the tested gas + 0.1% N₂; absorption + Rayleigh). Fixed system
parameters: Mp = 0.47 M⊕, Rs = 0.3155 R☉, Teff = 3415 K (Bello-Arufe 2025 /
Cadieux 2025). TauREx 3.3.2 with ExoMol R = 15 000 cross-sections
(SO₂: ExoAmes; CO₂: UCL-4000); PLATON 5.4 with its native opacity grid via
`custom_abundances` (PLATON 5.4 has no per-gas VMR API). Forward models are
bandpass-averaged onto each data grid exactly as in the Study 1 pair
statistic.

**Sampling and significance.** dynesty 3.1.0 static nested sampler,
nlive = 400, dlogz = 0.1, rstate seeded (42) — every cell of the matrix is
deterministic given the pinned package versions. σ from lnB via the Bayes-factor → p-value → σ
mapping of Benneke & Seager (2013), the one used in the paper; lnB ≤ 0 is
reported as 0.0σ.

**Cross-code consistency.** lnZ of the `flat` model (independent of the
radiative-transfer code) agrees between TauREx and PLATON to better than
1e-10 for all 9 spectra. On the negative control, lnB is negative in both
codes for both gases (|ΔlnB| ≤ 0.05). In the decision regime the codes agree
to |ΔlnB| ≲ 0.75 (Δσ ≤ ~0.2); the largest differences sit where they cannot
change the conclusion: own-visit-4 SO₂ (lnB 9.66 vs 12.09 — both "strong
preference") and FIREFLy CO₂ (lnB −2.84 vs −4.25 — both "strongly
disfavored"). The apparent FIREFLy SO₂ discrepancy (1.33σ vs 0.0σ) is a
threshold artifact of the σ mapping at lnB ≈ 0 (0.52 vs −0.03) — both codes
say "no preference".

**Caveats.** Our own-extraction errors are white out-of-transit scatter only
(no correlated-noise model), so the own-spectrum σ values are upper bounds.
Models are 3-parameter isothermal pure-gas atmospheres without clouds or
hazes — matched to the paper's comparison, not exhaustive. Only SO₂ and CO₂
are tested (the molecules discussed by the paper). System parameters are
fixed (their uncertainties are not propagated). The paper also fitted
detector offsets and other visit combinations; our re-extraction treats
visits separately.

## Results & data

| file | contents |
|---|---|
| `results/sigma_matrix.csv` | the 36 model-vs-flat comparisons: lnB, σ, best-fit parameters |
| `results/sigma_matrix.png` | the σ matrix as a heatmap |
| `results/evidences_taurex.csv`, `results/evidences_platon.csv` | all 54 retrieval runs: lnZ ± error, best-fit parameters, likelihood calls, runtime |
| `results/fit_<spectrum>.png` | each spectrum with its best-fit flat/SO₂/CO₂ models (gas curves: TauREx best fits) |
| `results/spectra/*.csv` | the 9 input spectrum variants used by the matrix |

Column definitions: [DATA_DICTIONARY.md](DATA_DICTIONARY.md).
Reproduce from scratch: [docs/REPRODUCE.md](../../docs/REPRODUCE.md),
section 7 (the matrix takes ~1–2 days of CPU time on a desktop PC).
