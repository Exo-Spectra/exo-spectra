"""Gap #1 verification, step 2: Bayesian model comparison on every L 98-59 b
spectrum variant with two independent retrieval codes.

Matrix: {pub_firefly, pub_eureka, own_avg, own_avg_infl, own_v1..v4,
         synth_flat} x {taurex, platon} x {flat, so2, co2}

Per (spectrum, code, model): nested sampling (dynesty) -> ln Z.
sigma(model vs flat) via the Bayes-factor -> p-value -> sigma mapping of
Benneke & Seager 2013 (as used by Bello-Arufe et al. 2025).

Models (isothermal, well-mixed):
  flat  1 param  : constant depth
  so2   3 params : Rp [R_earth], T [K], log10 P_surf [Pa]; 99.9% SO2 + N2
  co2   3 params : same, 99.9% CO2 + N2
Fixed: Mp = 0.47 M_earth, Rs = 0.3155 R_sun, Teff = 3415 K
(Bello-Arufe et al. 2025 / Cadieux et al. 2025).

Usage:
  python src/l9859b_retrieval.py --code taurex [--spectra own_avg,pub_eureka]
                                 [--models flat,so2,co2] [--nlive 400]
Results appended to reports/l9859b_verification/evidences_<code>.csv (idempotent:
existing (spectrum, code, model) rows are skipped unless --redo).
"""
import argparse
import json
import os
import sys
import time
import warnings
from pathlib import Path

# cap BLAS pools before numpy import — must coexist with the survey
# downloader and a second retrieval code on a loaded 32 GB system
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed" / "l9859b"
OUT = ROOT / "reports" / "l9859b_verification"
XSEC = ROOT / "data" / "opacities" / "taurex_ready"

# system parameters (Bello-Arufe 2025; stellar radius Cadieux 2025)
RS_SUN = 0.3155
TEFF = 3415.0
MP_ME = 0.47
RP_RE = 0.85
RJUP_RE = 0.0892147   # 1 R_earth in R_jup
MJUP_ME = 1 / 317.83  # 1 M_earth in M_jup

SPECTRA_ALL = ["pub_firefly", "pub_eureka", "own_avg", "own_avg_infl",
               "own_v1", "own_v2", "own_v3", "own_v4", "synth_flat"]
MODELS_ALL = ["flat", "so2", "co2"]

# priors
PRIOR_RP = (0.70, 1.00)      # R_earth
PRIOR_T = (300.0, 1000.0)    # K
PRIOR_LOGP = (1.0, 7.0)      # log10 Pa  (1e-4 .. 100 bar)


def load_data(name):
    df = pd.read_csv(DATA / f"{name}.csv")
    return (df["wave_um"].to_numpy(), df["dwave_um"].to_numpy(),
            df["depth_ppm"].to_numpy(), df["err_ppm"].to_numpy())


_BIN_CACHE = {}


def bin_model(wave_model, depth_model, wave, dwave):
    """Average the model over each data bandpass (same as compare.py).

    Bin edge indices are cached per (model grid, data grid) — the forward
    model grid is static within a retrieval, so this runs searchsorted once.
    """
    key = (len(wave_model), float(wave_model[0]), len(wave), float(wave[0]))
    idx = _BIN_CACHE.get(key)
    if idx is None:
        lo = np.searchsorted(wave_model, wave - dwave / 2, side="left")
        hi = np.searchsorted(wave_model, wave + dwave / 2, side="right")
        idx = (lo, hi)
        _BIN_CACHE[key] = idx
    lo, hi = idx
    csum = np.concatenate(([0.0], np.cumsum(depth_model)))
    n = (hi - lo).astype(float)
    out = np.where(n > 0, (csum[hi] - csum[lo]) / np.maximum(n, 1), np.nan)
    bad = ~np.isfinite(out)
    if bad.any():
        out[bad] = np.interp(wave[bad], wave_model, depth_model)
    return out


def lnb_to_sigma(lnB):
    """Benneke & Seager 2013: Bayes factor -> p-value -> n_sigma."""
    from scipy.optimize import brentq
    from scipy.special import erfcinv
    if lnB <= 0:
        return 0.0
    B = np.exp(min(lnB, 700))
    # B = -1/(e p ln p), valid for p < 1/e
    f = lambda p: -1.0 / (np.e * p * np.log(p)) - B
    try:
        p = brentq(f, 1e-300, 1 / np.e - 1e-12)
    except ValueError:
        return np.inf
    return float(np.sqrt(2) * erfcinv(2 * p))


# ---------------------------------------------------------------- taurex ---
class TaurexForward:
    """One pure-gas isothermal transmission model, parameters mutated in place."""

    def __init__(self, gas):
        from taurex.cache import OpacityCache
        OpacityCache().set_opacity_path(str(XSEC))
        from taurex.chemistry import TaurexChemistry, ConstantGas
        from taurex.planet import Planet
        from taurex.stellar import BlackbodyStar
        from taurex.temperature import Isothermal
        from taurex.model import TransmissionModel
        from taurex.contributions import (AbsorptionContribution,
                                          RayleighContribution)

        self.planet = Planet(planet_mass=MP_ME * MJUP_ME,
                             planet_radius=RP_RE * RJUP_RE)
        star = BlackbodyStar(temperature=TEFF, radius=RS_SUN)
        chem = TaurexChemistry(fill_gases=["N2"], ratio=1.0)
        chem.addGas(ConstantGas(gas, mix_ratio=0.999))
        self.tm = TransmissionModel(
            planet=self.planet, star=star, chemistry=chem,
            temperature_profile=Isothermal(T=600.0),
            atm_min_pressure=1e-1, atm_max_pressure=1e5, nlayers=50)
        self.tm.add_contribution(AbsorptionContribution())
        self.tm.add_contribution(RayleighContribution())
        self.tm.build()
        # native grid covering the data range with margin
        lam = np.exp(np.arange(np.log(2.55), np.log(5.35), 1 / 15000))
        self.wn = np.sort(10000.0 / lam)
        self.wave_model = 10000.0 / self.wn

    def depths_ppm(self, rp_re, T, logp_pa):
        self.tm["planet_radius"] = rp_re * RJUP_RE
        self.tm["T"] = T
        self.tm["atm_max_pressure"] = 10.0 ** logp_pa
        wn_out, depth, _, _ = self.tm.model(wngrid=self.wn)
        # taurex may pad the requested grid -> use the grid it returns;
        # sort ascending in wavelength for the np.interp fallback in bin_model
        wave = 10000.0 / wn_out
        order = np.argsort(wave)
        self.wave_model = wave[order]
        return depth[order] * 1e6


# ---------------------------------------------------------------- platon ---
class PlatonForward:
    """PLATON 5.4 TransitDepthCalculator, single-gas custom abundance grid.

    v5.4 has no gases=/vmrs= API -> build a custom_abundances dict: 99.9%
    of the tested gas + 0.1% N2 on the (T, P) grid of the solar template.
    """

    def __init__(self, gas):
        from platon.transit_depth_calculator import TransitDepthCalculator
        from platon.abundance_getter import AbundanceGetter
        from platon.constants import R_sun, R_earth, M_earth
        self.R_sun, self.R_earth, self.M_earth = R_sun, R_earth, M_earth
        self.calc = TransitDepthCalculator()
        ab = AbundanceGetter().get(0.0)   # solar template -> grid shapes
        if gas not in ab:
            raise ValueError(f"PLATON: no opacity for {gas}; "
                             f"available: {sorted(ab)}")
        for k in ab:
            ab[k] = np.zeros_like(ab[k])
        ab[gas] += 0.999
        ab["N2"] += 0.001
        self.abund = ab
        self.wave_model = None

    def depths_ppm(self, rp_re, T, logp_pa):
        wave, depth = self.calc.compute_depths(
            self.R_sun * RS_SUN, MP_ME * self.M_earth, rp_re * self.R_earth,
            T, logZ=None, CO_ratio=None, custom_abundances=self.abund,
            cloudtop_pressure=10.0 ** logp_pa, full_output=False)
        wave_um = np.asarray(wave, dtype=float) * 1e6
        depth = np.asarray(depth, dtype=float)
        m = (wave_um > 2.55) & (wave_um < 5.35)
        self.wave_model = wave_um[m]
        return depth[m] * 1e6


# -------------------------------------------------------------- retrieval ---
def run_one(code, spec_name, model_name, nlive, forward_cache):
    import dynesty

    wave, dwave, depth, err = load_data(spec_name)

    if model_name == "flat":
        lo, hi = depth.min() - 5 * err.max(), depth.max() + 5 * err.max()

        def prior(u):
            return np.array([lo + u[0] * (hi - lo)])

        def loglike(theta):
            r = (depth - theta[0]) / err
            return -0.5 * np.sum(r**2 + np.log(2 * np.pi * err**2))

        ndim = 1
    else:
        gas = model_name.upper()
        key = (code, gas)
        if key not in forward_cache:
            forward_cache[key] = (TaurexForward(gas) if code == "taurex"
                                  else PlatonForward(gas))
        fwd = forward_cache[key]

        def prior(u):
            return np.array([
                PRIOR_RP[0] + u[0] * (PRIOR_RP[1] - PRIOR_RP[0]),
                PRIOR_T[0] + u[1] * (PRIOR_T[1] - PRIOR_T[0]),
                PRIOR_LOGP[0] + u[2] * (PRIOR_LOGP[1] - PRIOR_LOGP[0]),
            ])

        fails = {"n": 0}

        def loglike(theta):
            try:
                model = fwd.depths_ppm(*theta)
                fails["n"] = 0
            except Exception:
                # tolerate sporadic numerical failures, but a broken forward
                # model (e.g. missing opacity file) must abort loudly instead
                # of silently returning -1e300 for every sample
                fails["n"] += 1
                if fails["n"] > 50:
                    raise
                return -1e300
            binned = bin_model(fwd.wave_model, model, wave, dwave)
            r = (depth - binned) / err
            return -0.5 * np.sum(r**2 + np.log(2 * np.pi * err**2))

        ndim = 3

    t0 = time.time()
    sampler = dynesty.NestedSampler(loglike, prior, ndim, nlive=nlive,
                                    rstate=np.random.default_rng(42))
    sampler.run_nested(dlogz=0.1, print_progress=False)
    res = sampler.results
    dt = time.time() - t0

    imax = int(np.argmax(res.logl))
    best = res.samples[imax]
    row = {
        "spectrum": spec_name, "code": code, "model": model_name,
        "n_points": len(wave), "lnZ": float(res.logz[-1]),
        "lnZ_err": float(res.logzerr[-1]),
        "best_params": json.dumps([round(float(x), 4) for x in best]),
        "max_lnL": float(res.logl[imax]), "ncall": int(np.sum(res.ncall)),
        "runtime_s": round(dt, 1),
    }
    print(f"[{code}|{spec_name}|{model_name}] lnZ={row['lnZ']:.2f}"
          f"+-{row['lnZ_err']:.2f}  best={row['best_params']}"
          f"  ({dt:.0f} s, {row['ncall']} calls)")
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--code", required=True, choices=["taurex", "platon"])
    ap.add_argument("--spectra", default=",".join(SPECTRA_ALL))
    ap.add_argument("--models", default=",".join(MODELS_ALL))
    ap.add_argument("--nlive", type=int, default=400)
    ap.add_argument("--redo", action="store_true")
    ap.add_argument("--max-runs", type=int, default=0,
                    help="exit after N retrievals (fresh process per batch "
                         "works around a slow memory leak in long runs)")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    # per-code file: two codes run concurrently, a shared CSV loses updates
    csv = OUT / f"evidences_{args.code}.csv"
    done = pd.read_csv(csv) if csv.exists() else pd.DataFrame(
        columns=["spectrum", "code", "model"])

    forward_cache = {}
    n_ran = 0
    pending_left = False
    for spec in args.spectra.split(","):
        for model in args.models.split(","):
            exists = ((done["spectrum"] == spec) & (done["code"] == args.code)
                      & (done["model"] == model)).any()
            if exists and not args.redo:
                continue
            if args.max_runs and n_ran >= args.max_runs:
                pending_left = True
                break
            row = run_one(args.code, spec, model, args.nlive, forward_cache)
            n_ran += 1
            done = done[~((done["spectrum"] == spec)
                          & (done["code"] == args.code)
                          & (done["model"] == model))]
            done = pd.concat([done, pd.DataFrame([row])], ignore_index=True)
            done.to_csv(csv, index=False)
        if pending_left:
            break

    if not pending_left:
        (OUT / f"DONE_{args.code}.flag").write_text("all cells complete\n")
        print(f"ALL DONE ({args.code})")

    # summary: sigma vs flat within each (spectrum, code)
    print("\n=== sigma(model vs flat) ===")
    for (spec, code), g in done.groupby(["spectrum", "code"]):
        if code != args.code:
            continue
        flat = g[g["model"] == "flat"]
        if flat.empty:
            continue
        lnz0 = float(flat["lnZ"].iloc[0])
        for _, r in g[g["model"] != "flat"].iterrows():
            lnb = float(r["lnZ"]) - lnz0
            print(f"{spec:14s} {code:7s} {r['model']:4s} "
                  f"lnB={lnb:+7.2f}  sigma={lnb_to_sigma(lnb):.2f}")


if __name__ == "__main__":
    main()
