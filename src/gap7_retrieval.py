"""Gap #7 (retrieval robustness), stage 2: Bayesian model comparison on PAIRS
of spectra flagged as structurally discordant (or as a consistent control) by
our catalogs — how much do atmospheric conclusions change when you swap the
reduction / epoch / pipeline / instrument of the same planet?

Input pairs: data/processed/gap7_selected_pairs.csv (gap7_select_targets.py).
System parameters: data/processed/retrieval_targets_tap.csv (NASA TAP).

Per (pair, side, code, model, err_variant): nested sampling (dynesty) -> lnZ.
sigma(model vs flat) via Benneke & Seager 2013, as in l9859b_retrieval.py.

Models (isothermal, trace gas in an H2/He fill, He/H2 = 0.17 — uniform across
all targets; rocky-planet caveat goes in the write-up):
  flat         1 param : constant depth
  <gas>        3 params: Rp [R_e], T [K], log10 X_gas          (no clouds)
  <gas>_cloud  4 params: + log10 P_cloud [Pa]                  (gray cloud)
First molecule listed for the pair = primary (gets the _cloud variant),
the rest are secondary (clear-sky only). Both spectra of a pair are CUT to
the pair's common wavelength range so each retrieval sees the same
spectral information.

Error-axis variants (err_variant): "raw" for everything; "infl" additionally
for pairs listed in INFLATE (errors x sqrt(chi2_red of the pair) — the honest
propagation of the measured epoch discordance, as own_avg_infl in study #1).

Usage:
  python src/gap7_retrieval.py --code taurex [--pairs gj1132b,k218b]
      [--models all] [--nlive 400] [--max-runs N] [--redo] [--list]
Results appended per code to reports/retrieval_robustness/evidences_<code>.csv
(idempotent: existing (spectrum, code, model, err_variant) rows are skipped).
"""
import argparse
import json
import os
import sys
import time
import warnings
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import spectra_io

ROOT = Path(__file__).resolve().parents[1]
# marathon guard: max likelihood calls per nested run (non-flat models).
# 2026-08-25: tightened 300k -> 100k after taurex burned 3 interrupted
# attempts on one K2-18 cell -- its forward model runs only ~5-10 calls/s,
# so 300k was an 8-16 h ceiling, never reached before the next interruption.
# 100k = ~3-5 h worst case at taurex speed, minutes-to-1h for platon.
MAXCALL = 100_000
# dynesty checkpoint cadence [s] for non-flat runs (2026-09-03); worst case a
# stop loses ~5 min of a marathon instead of the whole 25-50 h cell
CHECKPOINT_EVERY = 300
PAIRS_CSV = ROOT / "data" / "processed" / "gap7_selected_pairs.csv"
TAP_CSV = ROOT / "data" / "processed" / "retrieval_targets_tap.csv"
SPECTRA_DIR = ROOT / "data" / "spectra"
XSEC = ROOT / "data" / "opacities" / "taurex_ready"
OUT = ROOT / "reports" / "retrieval_robustness"

RJUP_RE = 0.0892147
MJUP_ME = 1 / 317.83
HE_H2 = 0.17
LAM_MARGIN = 1.05       # model range = data range extended by this factor

# priors (uniform across targets; scaled by measured Rp / T_eq)
PRIOR_RP_FACT = (0.5, 1.7)
PRIOR_LOGX = (-9.0, -0.3)
PRIOR_LOGPC = (0.0, 7.0)     # log10 cloud-top pressure [Pa]
T_FACT = (0.4, 1.6)

# short pair keys (CLI + spectrum keys in the CSV)
PAIR_KEY = {
    ("L 98-59 b", "reduction"): "l9859b",
    ("GJ 1132 b", "reduction"): "gj1132b",
    ("WASP-12 b", "reduction"): "wasp12b",
    ("WASP-39 b", "instrument"): "wasp39b",
    ("TOI-5205 b", "pipeline"): "toi5205b_pipe",
    ("TOI-5205 b", "epoch"): "toi5205b_epoch",
    ("HD 189733 b", "reduction"): "hd189733b",
    ("K2-18 b", "epoch"): "k218b",
    ("LHS 1140 b", "control_consistent"): "lhs1140b",
}
# pairs that additionally run with inflated errors: key -> chi2_red of the pair
INFLATE = {"k218b": 2.048883}
# pairs fully covered by study #1 — never re-run here, results imported in
# the analysis stage (l9859b: pub_firefly vs pub_eureka, all models, both codes)
SKIP_DONE = {"l9859b"}


# ------------------------------------------------------------------ data ---
def load_any(path_str):
    """Load one spectrum (archive .tbl / survey CSV / study#1 CSV) -> ppm."""
    p = Path(path_str)
    if not p.is_absolute():
        p = ROOT / p
    if p.suffix == ".tbl":
        df = spectra_io.load_spectrum(str(SPECTRA_DIR / p.name), "Transmission")
        err = 0.5 * (df.err_hi + df.err_lo)
        ok = (df.lim == 0) & np.isfinite(df.value) & np.isfinite(err) & (err > 0)
        df, err = df[ok], err[ok]
        wave, dwave = df.wave.to_numpy(), df.dwave.to_numpy()
        depth, err = df.value.to_numpy() * 1e4, err.to_numpy() * 1e4
    else:
        df = pd.read_csv(p)
        if "depth_pct" in df.columns:       # survey CSV (percent units)
            wave, dwave = df.wave.to_numpy(), df.dwave.to_numpy()
            depth, err = df.depth_pct.to_numpy() * 1e4, df.err_pct.to_numpy() * 1e4
        else:                               # study#1 CSV (already ppm)
            wave, dwave = df.wave_um.to_numpy(), df.dwave_um.to_numpy()
            depth, err = df.depth_ppm.to_numpy(), df.err_ppm.to_numpy()
    order = np.argsort(wave)
    return wave[order], dwave[order], depth[order], err[order]


def build_cells():
    """Expand the pair catalog into retrieval cells.

    Returns list of dicts: one per (pair, side) spectrum with its cut data,
    system params, molecule list and lambda range.
    """
    pairs = pd.read_csv(PAIRS_CSV)
    tap = pd.read_csv(TAP_CSV).set_index("pl_name")
    cells = []
    for _, r in pairs.iterrows():
        key = PAIR_KEY[(r.planet, r.axis)]
        if key in SKIP_DONE:
            continue
        st = tap.loc[r.planet]
        data = {}
        for side, f in (("a", r.file_a), ("b", r.file_b)):
            data[side] = load_any(f)
        # cut both sides to the pair's common wavelength range
        lo = max(d[0].min() for d in data.values())
        hi = min(d[0].max() for d in data.values())
        mols = [m.strip().upper() for m in r.molecules.split(",")]
        models = ["flat", mols[0].lower(), f"{mols[0].lower()}_cloud"] + \
                 [m.lower() for m in mols[1:]]
        for side in ("a", "b"):
            wave, dwave, depth, err = data[side]
            m = (wave >= lo - 1e-9) & (wave <= hi + 1e-9)
            variants = ["raw"] + (["infl"] if key in INFLATE else [])
            cells.append(dict(
                pair=key, planet=r.planet, axis=r.axis, side=side,
                spectrum=f"{key}__{side}", file=str(r.file_a if side == "a" else r.file_b),
                wave=wave[m], dwave=dwave[m], depth=depth[m], err=err[m],
                lam_lo=lo / LAM_MARGIN, lam_hi=hi * LAM_MARGIN,
                mp_me=float(st.pl_bmasse), rp_re=float(st.pl_rade),
                rs_sun=float(st.st_rad), teff=float(st.st_teff),
                teq=float(st.pl_eqt), models=models, variants=variants,
            ))
    return cells


# ---------------------------------------------------------------- taurex ---
class TaurexForward:
    """Isothermal trace gas + gray clouds over an H2/He fill."""

    def __init__(self, gas, cell):
        from taurex.cache import OpacityCache
        OpacityCache().set_opacity_path(str(XSEC))
        from taurex.chemistry import TaurexChemistry, ConstantGas
        from taurex.planet import Planet
        from taurex.stellar import BlackbodyStar
        from taurex.temperature import Isothermal
        from taurex.model import TransmissionModel
        from taurex.contributions import (AbsorptionContribution,
                                          RayleighContribution,
                                          SimpleCloudsContribution)
        self.gas = gas
        planet = Planet(planet_mass=cell["mp_me"] * MJUP_ME,
                        planet_radius=cell["rp_re"] * RJUP_RE)
        star = BlackbodyStar(temperature=cell["teff"], radius=cell["rs_sun"])
        chem = TaurexChemistry(fill_gases=["H2", "He"], ratio=HE_H2)
        chem.addGas(ConstantGas(gas, mix_ratio=1e-4))
        self.tm = TransmissionModel(
            planet=planet, star=star, chemistry=chem,
            temperature_profile=Isothermal(T=max(cell["teq"], 300.0)),
            atm_min_pressure=1e-1, atm_max_pressure=1e6, nlayers=50)
        self.tm.add_contribution(AbsorptionContribution())
        self.tm.add_contribution(RayleighContribution())
        self.tm.add_contribution(SimpleCloudsContribution(clouds_pressure=1e7))
        self.tm.build()
        lam = np.exp(np.arange(np.log(cell["lam_lo"]), np.log(cell["lam_hi"]),
                               1 / 15000))
        self.wn = np.sort(10000.0 / lam)
        self.wave_model = 10000.0 / self.wn

    def depths_ppm(self, rp_re, T, logx, logpc=7.0):
        self.tm["planet_radius"] = rp_re * RJUP_RE
        self.tm["T"] = T
        self.tm[self.gas] = 10.0 ** logx
        self.tm["clouds_pressure"] = 10.0 ** logpc
        wn_out, depth, _, _ = self.tm.model(wngrid=self.wn)
        wave = 10000.0 / wn_out
        order = np.argsort(wave)
        self.wave_model = wave[order]
        return depth[order] * 1e6


# ---------------------------------------------------------------- platon ---
class PlatonForward:
    """PLATON 5.4 TransitDepthCalculator, custom trace gas + H2/He."""

    def __init__(self, gas, cell):
        from platon.transit_depth_calculator import TransitDepthCalculator
        from platon.abundance_getter import AbundanceGetter
        from platon.constants import R_sun, R_earth, M_earth
        self.R_sun, self.R_earth, self.M_earth = R_sun, R_earth, M_earth
        self.gas = gas
        self.mp_me, self.rs_sun = cell["mp_me"], cell["rs_sun"]
        self.lam_lo, self.lam_hi = cell["lam_lo"], cell["lam_hi"]
        self.calc = TransitDepthCalculator()
        ab = AbundanceGetter().get(0.0)
        if gas not in ab:
            raise ValueError(f"PLATON: no opacity for {gas}")
        for k in ab:
            ab[k] = np.zeros_like(ab[k])
        self.abund = ab
        self.wave_model = None

    def depths_ppm(self, rp_re, T, logx, logpc=7.0):
        x = 10.0 ** logx
        h2 = (1.0 - x) / (1.0 + HE_H2)
        for k in self.abund:
            self.abund[k][:] = 0.0
        self.abund[self.gas] += x
        self.abund["H2"] += h2
        self.abund["He"] += h2 * HE_H2
        wave, depth = self.calc.compute_depths(
            self.R_sun * self.rs_sun, self.mp_me * self.M_earth,
            rp_re * self.R_earth, T, logZ=None, CO_ratio=None,
            custom_abundances=self.abund,
            cloudtop_pressure=10.0 ** logpc, full_output=False)
        wave_um = np.asarray(wave, dtype=float) * 1e6
        depth = np.asarray(depth, dtype=float)
        m = (wave_um > self.lam_lo) & (wave_um < self.lam_hi)
        self.wave_model = wave_um[m]
        return depth[m] * 1e6


# ------------------------------------------------------- shared machinery ---
# bin_model and lnb_to_sigma are identical to study #1 — import, don't fork
from l9859b_retrieval import bin_model, lnb_to_sigma  # noqa: E402

# dynesty checkpointing (2026-09-03) pickles the whole sampler, prior and
# likelihood included, so these are module-level classes, not closures. The
# heavy forward model (taurex/platon opacity grids) is NOT pickled: it lives
# in this cache and is rebuilt lazily from the cell after a restore (seconds,
# vs the hours a marathon cell costs).
_FWD_CACHE = {}


def get_forward(code, gas, cell):
    key = (code, gas, cell["spectrum"])
    if key not in _FWD_CACHE:
        _FWD_CACHE.clear()   # forwards hold big grids; keep one at a time
        _FWD_CACHE[key] = (TaurexForward(gas, cell) if code == "taurex"
                           else PlatonForward(gas, cell))
    return _FWD_CACHE[key]


class UniformPrior:
    def __init__(self, lo, hi):
        self.lo = np.asarray(lo, dtype=float)
        self.hi = np.asarray(hi, dtype=float)

    def __call__(self, u):
        return self.lo + u * (self.hi - self.lo)


class FlatLogLike:
    def __init__(self, depth, err):
        self.depth, self.err = depth, err

    def __call__(self, theta):
        r = (self.depth - theta[0]) / self.err
        return -0.5 * np.sum(r**2 + np.log(2 * np.pi * self.err**2))


class GasLogLike:
    def __init__(self, code, gas, cell, err):
        self.code, self.gas, self.cell = code, gas, cell
        self.err = err          # variant-adjusted copy (raw or inflated)
        self.fails = 0

    def __call__(self, theta):
        fwd = get_forward(self.code, self.gas, self.cell)
        try:
            model = fwd.depths_ppm(*theta)
            self.fails = 0
        except Exception:
            self.fails += 1
            if self.fails > 50:
                raise
            return -1e300
        binned = bin_model(fwd.wave_model, model, self.cell["wave"],
                           self.cell["dwave"])
        r = (self.cell["depth"] - binned) / self.err
        return -0.5 * np.sum(r**2 + np.log(2 * np.pi * self.err**2))


def run_one(code, cell, model_name, variant, nlive):
    import dynesty

    wave, dwave, depth = cell["wave"], cell["dwave"], cell["depth"]
    err = cell["err"].copy()
    if variant == "infl":
        err *= np.sqrt(INFLATE[cell["pair"]])

    if model_name == "flat":
        lo, hi = depth.min() - 5 * err.max(), depth.max() + 5 * err.max()
        prior = UniformPrior([lo], [hi])
        loglike = FlatLogLike(depth, err)
        ndim = 1
    else:
        cloudy = model_name.endswith("_cloud")
        gas = model_name.removesuffix("_cloud").upper()
        get_forward(code, gas, cell)    # build/reuse now: fail fast, warm cache

        t_lo = max(100.0, T_FACT[0] * cell["teq"])
        t_hi = min(3400.0, T_FACT[1] * cell["teq"])
        rp_lo, rp_hi = (PRIOR_RP_FACT[0] * cell["rp_re"],
                        PRIOR_RP_FACT[1] * cell["rp_re"])
        lows = [rp_lo, t_lo, PRIOR_LOGX[0]]
        highs = [rp_hi, t_hi, PRIOR_LOGX[1]]
        if cloudy:
            lows.append(PRIOR_LOGPC[0])
            highs.append(PRIOR_LOGPC[1])
        prior = UniformPrior(lows, highs)
        loglike = GasLogLike(code, gas, cell, err)
        ndim = 4 if cloudy else 3

    t0 = time.time()
    # sample="rwalk": dynesty's auto choice ("unif" at this ndim) collapses
    # on flat spectra (the 2026-08-19 benchmark lesson: 1.18M calls / 16 h on
    # GJ 1214 b; A/B probe validated rwalk). Switched 2026-08-22 after the
    # taurex matrix spent >5 h CPU on one cell; unif-computed free-model rows
    # were dropped and recomputed for homogeneity (flat rows kept).
    # Marathon guard (2026-08-24, after a taurex cell burned 8.3 h without
    # converging): cap the likelihood calls; if the cap is hit, redo the cell
    # once at nlive/2 (still a complete dlogz=0.1 nested run, just noisier
    # lnZ -- recorded in nlive_used); if that caps too, keep the capped
    # result with capped=1 (lnZ then unreliable -- filter in the report).
    # Checkpointing (2026-09-03): marathon cells (up to ~50 h on taurex) used
    # to be all-or-nothing -- any stop/reboot lost the whole run (STOP #6
    # burned 2 x ~25 h). Non-flat runs now checkpoint every 5 min; an
    # interrupted cell resumes from its .save (rng state is saved too, so the
    # result matches an uninterrupted run). One file per (code, spectrum,
    # model, variant, nlive attempt), deleted once the attempt completes.
    tag = f"[{code}|{cell['spectrum']}|{model_name}|{variant}]"
    maxcall = MAXCALL if model_name != "flat" else 0
    nlive_used, capped = nlive, 0
    chk_used = []
    for attempt_nlive in (nlive, max(nlive // 2, 100)):
        nlive_used = attempt_nlive
        chk = None
        if model_name != "flat":
            chk_dir = OUT / "checkpoints"
            chk_dir.mkdir(parents=True, exist_ok=True)
            chk = chk_dir / (f"chk_{code}_{cell['spectrum']}_{model_name}"
                             f"_{variant}_n{attempt_nlive}.save")
            chk_used.append(chk)
        sampler, resume = None, False
        if chk is not None and chk.exists():
            try:
                sampler = dynesty.NestedSampler.restore(str(chk))
                resume = True
                print(f"{tag} resuming from checkpoint "
                      f"({int(sampler.ncall)} calls done)", flush=True)
            except Exception as e:
                print(f"{tag} checkpoint restore failed ({e!r}), "
                      f"starting fresh", flush=True)
                chk.unlink(missing_ok=True)
                sampler = None
        if sampler is None:
            sampler = dynesty.NestedSampler(loglike, prior, ndim,
                                            nlive=attempt_nlive,
                                            sample="rwalk",
                                            rstate=np.random.default_rng(42))
        # run_nested's maxcall only counts calls made by THIS invocation, so
        # subtract what the resumed sampler already spent: the marathon cap
        # stays cumulative across resumes.
        budget = maxcall or 10 ** 9
        if resume:
            budget = max(budget - int(sampler.ncall), 1)
        run_kw = dict(dlogz=0.1, print_progress=False, maxcall=budget,
                      resume=resume)
        if chk is not None:
            run_kw.update(checkpoint_file=str(chk),
                          checkpoint_every=CHECKPOINT_EVERY)
        sampler.run_nested(**run_kw)
        res = sampler.results
        capped = int(maxcall and int(np.sum(res.ncall)) >= maxcall)
        if not capped:
            break
        print(f"{tag} maxcall hit at nlive={attempt_nlive}, "
              f"{'retrying at nlive//2' if attempt_nlive == nlive else 'keeping capped result'}",
              flush=True)
    # delete checkpoints only once the whole CELL is done: a capped attempt's
    # final state stays on disk, so a stop during the nlive//2 retry does not
    # force redoing the first ~100k-call attempt (a finished run restores and
    # returns its results in seconds)
    for f in chk_used:
        f.unlink(missing_ok=True)
    dt = time.time() - t0

    imax = int(np.argmax(res.logl))
    best = res.samples[imax]
    # posterior mean/std from the importance weights -- the sensitivity
    # metric of stage 3 is Delta(param) in units of these sigmas
    wts = np.exp(res.logwt - res.logz[-1])
    wts /= wts.sum()
    pmean = np.average(res.samples, weights=wts, axis=0)
    pstd = np.sqrt(np.average((res.samples - pmean) ** 2, weights=wts,
                              axis=0))
    row = {
        "spectrum": cell["spectrum"], "pair": cell["pair"],
        "planet": cell["planet"], "axis": cell["axis"], "side": cell["side"],
        "code": code, "model": model_name, "err_variant": variant,
        "n_points": len(wave), "lnZ": float(res.logz[-1]),
        "lnZ_err": float(res.logzerr[-1]),
        "best_params": json.dumps([round(float(x), 4) for x in best]),
        "post_mean": json.dumps([round(float(x), 4) for x in pmean]),
        "post_std": json.dumps([round(float(x), 4) for x in pstd]),
        "max_lnL": float(res.logl[imax]), "ncall": int(np.sum(res.ncall)),
        "runtime_s": round(dt, 1),
        "nlive_used": nlive_used, "capped": capped,
    }
    print(f"[{code}|{row['spectrum']}|{model_name}|{variant}] "
          f"lnZ={row['lnZ']:.2f}+-{row['lnZ_err']:.2f} "
          f"best={row['best_params']} ({dt:.0f} s, {row['ncall']} calls)")
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--code", required=True, choices=["taurex", "platon"])
    ap.add_argument("--pairs", default="")
    ap.add_argument("--models", default="all")
    ap.add_argument("--nlive", type=int, default=400)
    ap.add_argument("--redo", action="store_true")
    ap.add_argument("--max-runs", type=int, default=0)
    ap.add_argument("--list", action="store_true",
                    help="print the planned matrix and exit")
    ap.add_argument("--out-suffix", default="",
                    help="shard mode: write to evidences_<code><suffix>.csv "
                         "(done-skip also honours the main csv); lets "
                         "disjoint --pairs shards run in parallel without "
                         "racing on one file")
    args = ap.parse_args()

    cells = build_cells()
    if args.pairs:
        # --pairs order is binding (2026-09-08): cells run in the order the
        # pairs are listed, so a shard can be told "fastest pairs first"
        # (before: pairs-CSV row order, which fronted the TOI-5205 marathons)
        order = {k: i for i, k in enumerate(args.pairs.split(","))}
        cells = [c for c in cells if c["pair"] in order]
        cells.sort(key=lambda c: order[c["pair"]])

    todo = []
    for c in cells:
        for model in c["models"]:
            if args.models != "all" and model not in args.models.split(","):
                continue
            for variant in c["variants"]:
                todo.append((c, model, variant))

    if args.list:
        n_slow = sum(1 for c, m, v in todo if m != "flat")
        for c, m, v in todo:
            print(f"{c['spectrum']:22s} {m:12s} {v:5s} n={len(c['wave'])}")
        print(f"\n{len(todo)} cells for --code {args.code} "
              f"({n_slow} non-flat, x2 codes = {2 * n_slow} slow runs total)")
        return

    OUT.mkdir(parents=True, exist_ok=True)
    csv = OUT / f"evidences_{args.code}{args.out_suffix}.csv"
    done = pd.read_csv(csv) if csv.exists() else pd.DataFrame(
        columns=["spectrum", "code", "model", "err_variant"])
    skip = done
    main_csv = OUT / f"evidences_{args.code}.csv"
    if args.out_suffix and main_csv.exists():
        skip = pd.concat([done, pd.read_csv(main_csv)], ignore_index=True)

    n_ran = 0
    for c, model, variant in todo:
        exists = ((skip["spectrum"] == c["spectrum"])
                  & (skip["code"] == args.code) & (skip["model"] == model)
                  & (skip["err_variant"] == variant)).any()
        if exists and not args.redo:
            continue
        if args.max_runs and n_ran >= args.max_runs:
            print("max-runs reached, exiting (fresh process resumes)")
            return
        row = run_one(args.code, c, model, variant, args.nlive)
        n_ran += 1
        done = done[~((done["spectrum"] == c["spectrum"])
                      & (done["code"] == args.code)
                      & (done["model"] == model)
                      & (done["err_variant"] == variant))]
        done = pd.concat([done, pd.DataFrame([row])], ignore_index=True)
        done.to_csv(csv, index=False)

    if not args.pairs and args.models == "all":
        # a filtered invocation must never mark the full matrix as done
        (OUT / f"DONE_{args.code}.flag").write_text("all cells complete\n")
        print(f"ALL DONE ({args.code})")
    elif args.pairs:
        # shard completion marker for the shard wrapper (NOT the matrix flag)
        print(f"SHARD DONE ({args.code} pairs={args.pairs})")
    else:
        print("filtered run complete (no DONE flag)")


if __name__ == "__main__":
    main()
