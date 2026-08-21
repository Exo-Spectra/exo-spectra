"""Transit model + fitting library for the limb-asymmetry study.

Pure numpy quadratic-limb-darkening transit model (no compiled deps):
the occulted flux is integrated radially over the annulus covered by the
planet disk, so the model is exact up to quadrature error (validated to
<10 ppm against pytransit's QuadraticModel, see __main__ self-test).

Kepler orbit with e/omega so the *genuine* ingress/egress duration
asymmetry of an eccentric orbit is modelled and cannot leak into the
fitted radius asymmetry.

Fitting: scipy least_squares with Gaussian priors appended as residuals;
uncertainty for the asymmetry statistic via cyclic-shift residual
bootstrap (preserves red noise) and a binned-RMS beta factor. Optional
emcee escalation. No I/O here — drivers live in limb_asym_run.py.

Times are BMJD_TDB float64 throughout (JWST INT_TIMES convention).
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2 as chi2_dist

TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------- limb dark
def ld_coeffs(q1: float, q2: float) -> tuple[float, float]:
    """Kipping (2013) q1/q2 in [0,1] -> quadratic u1/u2 (always physical)."""
    s = np.sqrt(q1)
    return 2.0 * s * q2, s * (1.0 - 2.0 * q2)


# ------------------------------------------------------------- occultation
def occult(z, k: float, u1: float = 0.0, u2: float = 0.0, n_r: int = 200):
    """Fraction of total stellar flux blocked by an opaque disk of radius k
    (stellar radii) at projected separation(s) z. Quadratic limb darkening
    I(mu) = 1 - u1(1-mu) - u2(1-mu)^2. Radial midpoint quadrature over the
    covered annulus; azimuthal overlap is analytic (arc angle)."""
    z = np.abs(np.asarray(z, float))
    out = np.zeros_like(z)
    hit = z < 1.0 + k
    if not hit.any():
        return out
    zh = z[hit]
    lo = np.maximum(zh - k, 0.0)
    hi = np.minimum(zh + k, 1.0)
    u = (np.arange(n_r) + 0.5) / n_r
    r = lo[:, None] + (hi - lo)[:, None] * u[None, :]
    dr = ((hi - lo) / n_r)[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        cosa = (zh[:, None] ** 2 + r ** 2 - k ** 2) / (2.0 * zh[:, None] * r)
    # z ~ 0: division blows up; clip maps fully-covered rings (r<k) to
    # alpha=pi and uncovered ones to 0, which is the correct limit
    alpha = np.arccos(np.clip(np.nan_to_num(cosa, nan=-1.0), -1.0, 1.0))
    mu = np.sqrt(np.maximum(1.0 - r ** 2, 0.0))
    inten = 1.0 - u1 * (1.0 - mu) - u2 * (1.0 - mu) ** 2
    f_tot = np.pi * (1.0 - u1 / 3.0 - u2 / 6.0)
    out[hit] = np.sum(2.0 * r * alpha * inten * dr, axis=1) / f_tot
    return out


# ------------------------------------------------------------------- orbit
def _kepler_E(M, e: float, n_iter: int = 12):
    """Vectorized Newton solve of Kepler's equation."""
    E = M + e * np.sin(M)
    for _ in range(n_iter):
        E -= (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
    return E


def impact_to_inc(b: float, aRs: float, e: float = 0.0, w_deg: float = 90.0):
    """Impact parameter -> inclination [rad]; b = aRs cos i (1-e^2)/(1+e sin w)."""
    cosi = b / aRs * (1.0 + e * np.sin(np.radians(w_deg))) / (1.0 - e ** 2)
    return np.arccos(np.clip(cosi, 0.0, 1.0))


def proj_sep(t, t0: float, per: float, aRs: float, inc: float,
             e: float = 0.0, w_deg: float = 90.0):
    """Projected star-planet separation [stellar radii] and planet-in-front
    mask. t0 = time of inferior conjunction (mid-transit). inc in radians."""
    t = np.asarray(t, float)
    if e < 1e-6:
        th = TWO_PI * (t - t0) / per
        z = aRs * np.sqrt(np.sin(th) ** 2 + (np.cos(inc) * np.cos(th)) ** 2)
        return z, np.cos(th) > 0
    w = np.radians(w_deg)
    nu0 = np.pi / 2.0 - w
    E0 = 2.0 * np.arctan(np.sqrt((1.0 - e) / (1.0 + e)) * np.tan(nu0 / 2.0))
    M0 = E0 - e * np.sin(E0)
    M = M0 + TWO_PI * (t - t0) / per
    E = _kepler_E(M, e)
    nu = 2.0 * np.arctan2(np.sqrt(1.0 + e) * np.sin(E / 2.0),
                          np.sqrt(1.0 - e) * np.cos(E / 2.0))
    r = aRs * (1.0 - e ** 2) / (1.0 + e * np.cos(nu))
    z = r * np.sqrt(1.0 - (np.sin(w + nu) * np.sin(inc)) ** 2)
    return z, np.sin(w + nu) > 0


def contacts(t0: float, per: float, aRs: float, inc: float, k: float,
             e: float = 0.0, w_deg: float = 90.0, window_d: float = 0.35):
    """Contact times T1..T4 (BMJD) by scanning z(t) on a fine grid around t0.
    nan where a contact does not occur (grazing)."""
    tt = t0 + np.linspace(-window_d, window_d, 40001)
    z, front = proj_sep(tt, t0, per, aRs, inc, e, w_deg)
    z = np.where(front, z, np.inf)

    def crossings(level):
        s = np.sign(z - level)
        idx = np.flatnonzero(np.diff(s) != 0)
        times = []
        for i in idx:
            f = (level - z[i]) / (z[i + 1] - z[i])
            times.append(tt[i] + f * (tt[i + 1] - tt[i]))
        return times

    o = crossings(1.0 + k)
    i = crossings(1.0 - k)
    t1 = o[0] if len(o) >= 2 else np.nan
    t4 = o[-1] if len(o) >= 2 else np.nan
    t2 = i[0] if len(i) >= 2 else np.nan
    t3 = i[-1] if len(i) >= 2 else np.nan
    return t1, t2, t3, t4


# ------------------------------------------------------------------- ramps
RAMP_NPAR = {"linear": 1, "quad": 2, "quad_exp": 4}


def ramp_flux(name: str, coeffs, x, xe):
    """Multiplicative systematics model. x = (t - t_ref)/span (conditioned),
    xe = t - t_min [d] (for the settling exponential)."""
    if name == "linear":
        return 1.0 + coeffs[0] * x
    if name == "quad":
        return 1.0 + coeffs[0] * x + coeffs[1] * x ** 2
    if name == "quad_exp":
        tau = max(coeffs[3], 1e-4)
        return ((1.0 + coeffs[0] * x + coeffs[1] * x ** 2)
                * (1.0 + coeffs[2] * np.exp(-xe / tau)))
    raise ValueError(name)


# ------------------------------------------------------------- model config
@dataclass
class TransitConfig:
    per: float                # orbital period [d], fixed
    e: float = 0.0            # eccentricity, fixed
    w_deg: float = 90.0       # argument of periastron [deg], fixed
    ramp: str = "quad"
    asym: bool = False        # k -> (k_in, k_eg) split at fitted t0
    t_ref: float = 0.0        # conditioning: x = (t - t_ref)/span
    span: float = 1.0
    t_min: float = 0.0
    n_r: int = 200

    def names(self):
        base = (["t0", "k_in", "k_eg"] if self.asym else ["t0", "k"])
        return base + ["aRs", "b", "q1", "q2"] + \
            [f"c{i}" for i in range(RAMP_NPAR[self.ramp])]


def model_flux(t, theta, cfg: TransitConfig):
    """Full model: limb-darkened transit x multiplicative ramp."""
    i = 0
    t0 = theta[0]
    if cfg.asym:
        k_in, k_eg = theta[1], theta[2]
        i = 3
    else:
        k_in = k_eg = theta[1]
        i = 2
    aRs, b, q1, q2 = theta[i:i + 4]
    coeffs = theta[i + 4:]
    u1, u2 = ld_coeffs(q1, q2)
    inc = impact_to_inc(b, aRs, cfg.e, cfg.w_deg)
    z, front = proj_sep(t, t0, cfg.per, aRs, inc, cfg.e, cfg.w_deg)
    tr = np.ones_like(z)
    for k, side in ((k_in, t <= t0), (k_eg, t > t0)):
        m = front & side & (z < 1.0 + k)
        if m.any():
            tr[m] = 1.0 - occult(z[m], k, u1, u2, cfg.n_r)
    x = (t - cfg.t_ref) / cfg.span
    xe = t - cfg.t_min
    return tr * ramp_flux(cfg.ramp, coeffs, x, xe)


# ---------------------------------------------------------------- fitting
@dataclass
class FitResult:
    theta: np.ndarray
    names: list
    errors: np.ndarray
    sigma: float              # per-point white scatter [rel. flux]
    chi2: float
    dof: int
    bic: float
    residuals: np.ndarray = field(repr=False)
    model: np.ndarray = field(repr=False)
    cfg: TransitConfig = None

    def __getitem__(self, name):
        return self.theta[self.names.index(name)]

    def err(self, name):
        return self.errors[self.names.index(name)]


def _residual_fn(t, f, sigma, cfg, priors):
    def fn(theta):
        res = (f - model_flux(t, theta, cfg)) / sigma
        extra = [(theta[i] - mu) / sd for i, mu, sd in priors]
        return np.concatenate([res, extra]) if extra else res
    return fn


def _bounds(cfg: TransitConfig, t0_guess, dur_d):
    lo, hi = [], []
    for n in cfg.names():
        if n == "t0":
            lo.append(t0_guess - 0.5 * dur_d)
            hi.append(t0_guess + 0.5 * dur_d)
        elif n in ("k", "k_in", "k_eg"):
            lo.append(1e-4)
            hi.append(0.5)
        elif n == "aRs":
            lo.append(1.0)
            hi.append(500.0)
        elif n == "b":
            lo.append(0.0)
            hi.append(1.3)
        elif n in ("q1", "q2"):
            lo.append(0.0)
            hi.append(1.0)
        elif n == "c3":          # quad_exp timescale [d]
            lo.append(1e-4)
            hi.append(0.5)
        else:                    # ramp polynomial coefficients
            lo.append(-0.5)
            hi.append(0.5)
    return np.array(lo), np.array(hi)


def robust_sigma(f):
    """Per-point white scatter from median absolute successive difference."""
    d = np.diff(f[np.isfinite(f)])
    return 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2.0)


def fit_transit(t, f, init: dict, cfg: TransitConfig, priors: dict,
                dur_d: float, n_restarts: int = 3, seed: int = 0):
    """Least-squares fit with Gaussian priors (appended residuals).

    init: {name: value} starting point (missing ramp coeffs start at 0).
    priors: {name: (mu, sd)} Gaussian priors, e.g. aRs/b from the archive.
    Returns FitResult with post-hoc rescaled errors.
    """
    ok = np.isfinite(f) & np.isfinite(t)
    t, f = t[ok], f[ok]
    names = cfg.names()
    theta0 = np.array([init.get(n, 0.0) for n in names])
    if cfg.ramp == "quad_exp" and init.get("c3", 0.0) == 0.0:
        theta0[names.index("c3")] = 0.01
    sigma = robust_sigma(f)
    prior_list = [(names.index(n), mu, sd) for n, (mu, sd) in priors.items()
                  if n in names]
    fn = _residual_fn(t, f, sigma, cfg, prior_list)
    lo, hi = _bounds(cfg, init["t0"], dur_d)
    theta0 = np.clip(theta0, lo + 1e-12, hi - 1e-12)

    rng = np.random.default_rng(seed)
    best = None
    for j in range(n_restarts):
        start = theta0.copy()
        if j:
            start[0] += rng.normal(0.0, 0.05 * dur_d)
            start = np.clip(start, lo + 1e-12, hi - 1e-12)
        sol = least_squares(fn, start, bounds=(lo, hi), x_scale="jac",
                            method="trf")
        if best is None or sol.cost < best.cost:
            best = sol

    model = model_flux(t, best.x, cfg)
    res = f - model
    ndata = len(f)
    npar = len(names)
    dof = max(ndata - npar, 1)
    chi2 = float(np.sum((res / sigma) ** 2))
    # post-hoc error scaling: covariance from J^T J times reduced chi2
    J = best.jac[:ndata]         # drop prior rows for the data covariance
    try:
        cov = np.linalg.inv(J.T @ J) * chi2 / dof
        errors = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    except np.linalg.LinAlgError:
        errors = np.full(npar, np.nan)
    bic = chi2 + npar * np.log(ndata)
    return FitResult(best.x, names, errors, sigma, chi2, dof, bic,
                     res, model, cfg), t, f


def fit_with_ramp_selection(t, f, init, priors, dur_d, per, e, w_deg,
                            ramps=("linear", "quad", "quad_exp"), n_r=200):
    """Symmetric fit per ramp model; returns (best FitResult, {ramp: bic})."""
    t_ref, span, t_min = float(np.median(t)), float(t.max() - t.min()), \
        float(t.min())
    results = {}
    for ramp in ramps:
        cfg = TransitConfig(per=per, e=e, w_deg=w_deg, ramp=ramp,
                            t_ref=t_ref, span=span, t_min=t_min, n_r=n_r)
        try:
            fit, tt, ff = fit_transit(t, f, init, cfg, priors, dur_d)
            results[ramp] = (fit, tt, ff)
        except Exception:
            continue
    if not results:
        raise RuntimeError("all ramp fits failed")
    bics = {r: v[0].bic for r, v in results.items()}
    best = min(bics, key=bics.get)
    return results[best], bics


def fit_asymmetric(t, f, sym: FitResult, priors, dur_d):
    """Refit with independent ingress/egress radii; everything else shared.
    Starts from the symmetric solution."""
    cfg = TransitConfig(per=sym.cfg.per, e=sym.cfg.e, w_deg=sym.cfg.w_deg,
                        ramp=sym.cfg.ramp, asym=True, t_ref=sym.cfg.t_ref,
                        span=sym.cfg.span, t_min=sym.cfg.t_min,
                        n_r=sym.cfg.n_r)
    init = {n: sym[n] for n in sym.names}
    init["k_in"] = init["k_eg"] = init.pop("k")
    return fit_transit(t, f, init, cfg, priors, dur_d, n_restarts=1)


def delta_depth_ppm(asym: FitResult) -> float:
    """Asymmetry statistic: egress - ingress depth [ppm]."""
    return (asym["k_eg"] ** 2 - asym["k_in"] ** 2) * 1e6


# ------------------------------------------------------- uncertainty tools
def beta_factor(t, res, bin_minutes=(5, 10, 15, 20, 30)):
    """Red-noise inflation: max over bin sizes of binned RMS vs the white
    expectation (Winn 2008 style). Time-ordered binning on a uniform grid."""
    res = np.asarray(res, float)
    ok = np.isfinite(res)
    t, res = t[ok], res[ok]
    s1 = np.std(res)
    if s1 == 0 or len(res) < 50:
        return 1.0
    betas = []
    for bm in bin_minutes:
        w = bm / (24.0 * 60.0)
        idx = np.floor((t - t[0]) / w).astype(int)
        nb = idx.max() + 1
        if nb < 8:
            continue
        sums = np.bincount(idx, weights=res, minlength=nb)
        cnts = np.bincount(idx, minlength=nb)
        good = cnts >= max(2, int(0.5 * np.median(cnts[cnts > 0])))
        if good.sum() < 8:
            continue
        means = sums[good] / cnts[good]
        n = cnts[good].mean()
        expected = s1 / np.sqrt(n) * np.sqrt(len(means) / max(len(means) - 1.0, 1.0))
        betas.append(np.std(means) / expected)
    return max(1.0, max(betas)) if betas else 1.0


def bootstrap_delta_depth(t, f, asym: FitResult, priors, dur_d,
                          n_boot: int = 300, seed: int = 1):
    """Cyclic-shift residual bootstrap (preserves red noise): shift the
    residuals of the asymmetric fit by a random offset, add back to the
    model, refit from the MAP start. Returns (dd_err_ppm, dd_samples)."""
    ok = np.isfinite(f) & np.isfinite(t)
    tt, _ = t[ok], f[ok]
    model, res = asym.model, asym.residuals
    rng = np.random.default_rng(seed)
    init = {n: asym[n] for n in asym.names}
    samples = []
    for _ in range(n_boot):
        shift = rng.integers(1, len(res) - 1)
        fb = model + np.roll(res, shift)
        try:
            fit, _, _ = fit_transit(tt, fb, init, asym.cfg, priors, dur_d,
                                    n_restarts=1)
            samples.append(delta_depth_ppm(fit))
        except Exception:
            continue
    samples = np.asarray(samples)
    if len(samples) < 20:
        return np.nan, samples
    lo, hi = np.percentile(samples, [15.865, 84.135])
    return 0.5 * (hi - lo), samples


def fold_residual_test(t, res, sigma, t0, t1, t4):
    """Method B: fold symmetric-fit residuals about t0 within [T1, T4] and
    compare ingress side vs mirrored egress side. Returns dict with the
    mean-difference z-score, chi2 p-value and a sign-test p-value."""
    ok = np.isfinite(res)
    t, res = t[ok], res[ok]
    inn = (t >= t1) & (t <= t4)
    ing = inn & (t < t0)
    egr = inn & (t > t0)
    if ing.sum() < 10 or egr.sum() < 10:
        return {"mb_n": 0, "mb_z": np.nan, "mb_p_chi2": np.nan,
                "mb_p_sign": np.nan}
    dt_ing = t0 - t[ing]
    order = np.argsort(t[egr] - t0)
    dt_egr = (t[egr] - t0)[order]
    r_egr = res[egr][order]
    lo, hi = max(dt_ing.min(), dt_egr.min()), min(dt_ing.max(), dt_egr.max())
    use = (dt_ing >= lo) & (dt_ing <= hi)
    if use.sum() < 10:
        return {"mb_n": 0, "mb_z": np.nan, "mb_p_chi2": np.nan,
                "mb_p_sign": np.nan}
    mirror = np.interp(dt_ing[use], dt_egr, r_egr)
    d = res[ing][use] - mirror
    n = len(d)
    sd = sigma * np.sqrt(2.0)          # both sides carry white noise
    z = d.mean() / (sd / np.sqrt(n))
    chi2 = float(np.sum((d / sd) ** 2))
    p_chi2 = float(chi2_dist.sf(chi2, n))
    npos = int((d > 0).sum())
    from scipy.stats import binomtest
    p_sign = float(binomtest(npos, n, 0.5).pvalue)
    return {"mb_n": n, "mb_z": float(z), "mb_p_chi2": p_chi2,
            "mb_p_sign": p_sign}


# ------------------------------------------------------------------- emcee
def mcmc_asym(t, f, asym: FitResult, priors, dur_d, n_walkers: int = 40,
              n_burn: int = 1500, n_steps: int = 3000, seed: int = 2):
    """Escalation path: emcee posterior for the asymmetric model with a
    free white-noise jitter. Returns (dd_med, dd_err, flat_dd_samples)."""
    import emcee
    ok = np.isfinite(f) & np.isfinite(t)
    t, f = t[ok], f[ok]
    cfg = asym.cfg
    names = asym.names
    prior_list = [(names.index(n), mu, sd) for n, (mu, sd) in priors.items()
                  if n in names]
    lo, hi = _bounds(cfg, asym["t0"], dur_d)
    lns0 = np.log(asym.sigma)

    def log_prob(theta):
        p, lns = theta[:-1], theta[-1]
        if np.any(p < lo) or np.any(p > hi) or not (lns0 - 2 < lns < lns0 + 2):
            return -np.inf
        s2 = np.exp(2.0 * lns)
        res = f - model_flux(t, p, cfg)
        ll = -0.5 * np.sum(res ** 2 / s2 + np.log(TWO_PI * s2))
        lp = sum(-0.5 * ((p[i] - mu) / sd) ** 2 for i, mu, sd in prior_list)
        return ll + lp

    rng = np.random.default_rng(seed)
    ndim = len(names) + 1
    center = np.append(asym.theta, lns0)
    scale = np.append(np.where(np.isfinite(asym.errors) & (asym.errors > 0),
                               asym.errors, 1e-4), 0.05)
    p0 = center + scale * 0.5 * rng.standard_normal((n_walkers, ndim))
    p0[:, :-1] = np.clip(p0[:, :-1], lo + 1e-12, hi - 1e-12)
    sampler = emcee.EnsembleSampler(n_walkers, ndim, log_prob)
    state = sampler.run_mcmc(p0, n_burn, progress=False)
    sampler.reset()
    sampler.run_mcmc(state, n_steps, progress=False)
    chain = sampler.get_chain(flat=True)
    i_in, i_eg = names.index("k_in"), names.index("k_eg")
    dd = (chain[:, i_eg] ** 2 - chain[:, i_in] ** 2) * 1e6
    lo_p, med, hi_p = np.percentile(dd, [15.865, 50.0, 84.135])
    return med, 0.5 * (hi_p - lo_p), dd


# -------------------------------------------------------------- self-test
if __name__ == "__main__":
    import sys
    import time
    from pathlib import Path

    print("== limb_asym_model self-test ==")

    # 1. quadrature convergence: n_r=200 vs n_r=20000 reference
    k, u1, u2 = 0.05, 0.35, 0.18
    zz = np.linspace(0.0, 1.0 + k, 4001)
    f200 = occult(zz, k, u1, u2, n_r=200)
    fref = occult(zz, k, u1, u2, n_r=20000)
    err_ppm = np.max(np.abs(f200 - fref)) * 1e6
    print(f"quadrature n_r=200 vs 20000: max |diff| = {err_ppm:.3f} ppm")
    assert err_ppm < 2.0, "quadrature error too large"

    # 2. against pytransit QuadraticModel (celerite stub on path)
    sys.path.insert(0, str(Path(__file__).parent / "_stubs"))
    try:
        from pytransit import QuadraticModel
        tm = QuadraticModel()
        tgrid = np.linspace(-0.06, 0.06, 2000)
        worst = 0.0
        for kk, b, uu1, uu2 in [(0.02, 0.1, 0.2, 0.1), (0.05, 0.3, 0.4, 0.2),
                                (0.074, 0.23, 0.5, 0.15), (0.03, 0.9, 0.3, 0.3)]:
            aRs, per = 25.0, 10.0
            inc = impact_to_inc(b, aRs)
            tm.set_data(tgrid)
            f_pt = tm.evaluate(k=kk, ldc=[uu1, uu2], t0=0.0, p=per, a=aRs,
                               i=float(inc))
            z, front = proj_sep(tgrid, 0.0, per, aRs, inc)
            f_own = 1.0 - occult(z, kk, uu1, uu2)
            worst = max(worst, np.max(np.abs(f_own - np.asarray(f_pt))) * 1e6)
        print(f"vs pytransit over 4 param sets: max |diff| = {worst:.3f} ppm")
        assert worst < 10.0, "disagreement with pytransit > 10 ppm"
    except ImportError as exc:
        print(f"pytransit unavailable ({exc}) — skipping cross-check")

    # 3. end-to-end: inject + recover a symmetric transit, then null dd
    rng = np.random.default_rng(42)
    per, aRs, b, kk = 24.737, 95.3, 0.23, 0.0739     # LHS 1140 b-like
    inc = impact_to_inc(b, aRs)
    t = 60873.0 + np.linspace(0.0, 0.28, 5000)       # 6.7 h visit
    t0_true = 60873.14
    z, front = proj_sep(t, t0_true, per, aRs, inc)
    u1t, u2t = ld_coeffs(0.5, 0.35)
    flux = 1.0 - occult(z, kk, u1t, u2t)
    flux *= 1.0 + 0.001 * ((t - t.mean()) / 0.28) - 0.0008 * ((t - t.mean()) / 0.28) ** 2
    noise = 250e-6
    f_obs = flux + rng.normal(0.0, noise, t.size)

    dur_d = 2.15 / 24.0
    init = {"t0": t0_true + 0.002, "k": 0.06, "aRs": aRs, "b": b,
            "q1": 0.4, "q2": 0.3}
    priors = {"aRs": (aRs, 0.1 * aRs), "b": (b, 0.1)}
    tic = time.time()
    (sym, tt, ff), bics = fit_with_ramp_selection(
        t, f_obs, init, priors, dur_d, per, 0.0, 90.0)
    dt_fit = time.time() - tic
    depth_fit = sym["k"] ** 2 * 1e6
    print(f"symmetric fit ({dt_fit:.1f} s, ramp={sym.cfg.ramp}, "
          f"BICs={ {r: round(v, 1) for r, v in bics.items()} }): "
          f"k={sym['k']:.5f} (true {kk}), depth={depth_fit:.0f} ppm "
          f"(true {kk**2*1e6:.0f}), t0 err={(sym['t0']-t0_true)*86400:.1f} s")
    assert abs(sym["k"] - kk) < 0.002, "k recovery off"
    assert abs(sym["t0"] - t0_true) * 86400 < 30, "t0 recovery off"

    tic = time.time()
    asym, ta, fa = fit_asymmetric(tt, ff, sym, priors, dur_d)
    dd = delta_depth_ppm(asym)
    dt_asym = time.time() - tic
    print(f"asymmetric fit ({dt_asym:.1f} s): dd = {dd:+.1f} ppm "
          f"(true 0), k_in={asym['k_in']:.5f} k_eg={asym['k_eg']:.5f}")

    tic = time.time()
    dd_err, samples = bootstrap_delta_depth(ta, fa, asym, priors, dur_d,
                                            n_boot=60)
    print(f"bootstrap x60 ({time.time()-tic:.1f} s): dd_err = {dd_err:.1f} ppm"
          f" -> dd/err = {dd/dd_err:+.2f} sigma (expect |.|<3)")
    assert abs(dd) < 4 * dd_err, "null injection came out significant"

    t1, t2, t3, t4 = contacts(sym["t0"], per, aRs,
                              impact_to_inc(sym["b"], sym["aRs"]), sym["k"])
    mb = fold_residual_test(ta, asym.residuals, asym.sigma, sym["t0"], t1, t4)
    print(f"contacts: T14 = {(t4-t1)*24:.3f} h (archive 2.15 h), "
          f"method B: z={mb['mb_z']:+.2f}, p_chi2={mb['mb_p_chi2']:.3f}")
    print("== all self-tests passed ==")
