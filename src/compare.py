"""Pair-comparison statistic for two spectra of the same planet & type.

Model: within the overlapping wavelength range, spectrum B equals spectrum A
plus a constant vertical offset c (free parameter absorbing instrument /
normalization systematics). We test whether the residuals after fitting c
are consistent with the quoted uncertainties.

    chi2 = sum_i w_i * (B_i - A_i - c)^2,  w_i = 1 / (sigma_A_i^2 + sigma_B_i^2)
    dof  = N - 1  (one fitted parameter)

Matching wavelength grids: the sparser spectrum defines the grid. For each of
its points we average the other spectrum's points falling inside the bandpass
[wave - dwave/2, wave + dwave/2] (inverse-variance weighted); if the bandpass
is empty/undefined we fall back to linear interpolation (errors interpolated
too — an approximation, flagged in n_interp).
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats


@dataclass
class PairResult:
    n: int              # matched points used
    n_interp: int       # points obtained by interpolation fallback
    offset: float       # fitted vertical offset c (B - A)
    offset_err: float
    chi2: float
    dof: int
    p_value: float      # P(chi2 >= observed | consistent) — small = discrepant
    chi2_no_offset: float
    p_no_offset: float
    # offset + slope model: B - A = c + m * (wave - mean(wave))
    slope: float | None = None
    slope_err: float | None = None
    chi2_slope: float | None = None
    dof_slope: int | None = None
    p_slope: float | None = None   # small = discrepant even allowing a linear trend


def _sym_sigma(df: pd.DataFrame) -> np.ndarray:
    """Symmetrized 1-sigma uncertainty; falls back to the one finite side."""
    hi, lo = df["err_hi"].to_numpy(float), df["err_lo"].to_numpy(float)
    sig = np.nanmean(np.vstack([hi, lo]), axis=0)
    return sig


def _match(sparse: pd.DataFrame, dense: pd.DataFrame):
    """Project `dense` onto `sparse`'s wavelength grid. Returns (a, sa, b, sb, n_interp)."""
    w_d = dense["wave"].to_numpy(float)
    v_d = dense["value"].to_numpy(float)
    s_d = _sym_sigma(dense)

    vals, sigs, interp_flags = [], [], []
    for _, row in sparse.iterrows():
        w, dw = row["wave"], row["dwave"]
        if np.isfinite(dw) and dw > 0:
            m = (w_d >= w - dw / 2) & (w_d <= w + dw / 2)
        else:
            m = np.zeros_like(w_d, bool)
        if m.sum() >= 1:
            wgt = 1.0 / s_d[m] ** 2
            vals.append(np.sum(v_d[m] * wgt) / np.sum(wgt))
            sigs.append(np.sqrt(1.0 / np.sum(wgt)))
            interp_flags.append(False)
        else:
            if not (w_d.min() <= w <= w_d.max()):
                vals.append(np.nan); sigs.append(np.nan); interp_flags.append(False)
                continue
            vals.append(np.interp(w, w_d, v_d))
            sigs.append(np.interp(w, w_d, s_d))
            interp_flags.append(True)

    a = sparse["value"].to_numpy(float)
    sa = _sym_sigma(sparse)
    wav = sparse["wave"].to_numpy(float)
    b, sb = np.asarray(vals), np.asarray(sigs)
    ok = np.isfinite(a) & np.isfinite(sa) & np.isfinite(b) & np.isfinite(sb) & (sa > 0) & (sb > 0)
    return wav[ok], a[ok], sa[ok], b[ok], sb[ok], int(np.sum(np.asarray(interp_flags)[ok]))


def compare_pair(spec_a: pd.DataFrame, spec_b: pd.DataFrame, min_points: int = 3) -> PairResult | None:
    """Compare two normalized spectra (usable points only). None if too few matched points."""
    # overlap region
    lo = max(spec_a["wave"].min(), spec_b["wave"].min())
    hi = min(spec_a["wave"].max(), spec_b["wave"].max())
    pad = 1e-9
    a_in = spec_a[(spec_a["wave"] >= lo - pad) & (spec_a["wave"] <= hi + pad)]
    b_in = spec_b[(spec_b["wave"] >= lo - pad) & (spec_b["wave"] <= hi + pad)]
    if len(a_in) == 0 or len(b_in) == 0:
        return None

    # sparser spectrum defines the grid
    if len(a_in) <= len(b_in):
        wav, a, sa, b, sb, n_interp = _match(a_in, spec_b)
        sign = 1.0   # offset means B - A
    else:
        wav, a, sa, b, sb, n_interp = _match(b_in, spec_a)
        sign = -1.0  # roles swapped -> flip offset back to B - A

    n = len(a)
    if n < min_points:
        return None

    diff = b - a
    var = sa**2 + sb**2
    w = 1.0 / var
    c = np.sum(diff * w) / np.sum(w)
    c_err = np.sqrt(1.0 / np.sum(w))
    chi2 = float(np.sum((diff - c) ** 2 * w))
    dof = n - 1
    chi2_no = float(np.sum(diff**2 * w))
    res = PairResult(
        n=n, n_interp=n_interp,
        offset=float(sign * c), offset_err=float(c_err),
        chi2=chi2, dof=dof,
        p_value=float(stats.chi2.sf(chi2, dof)),
        chi2_no_offset=chi2_no,
        p_no_offset=float(stats.chi2.sf(chi2_no, n)),
    )

    # offset + slope model (needs >=2 distinct wavelengths and n >= min_points+1)
    x = wav - np.average(wav, weights=w)
    if n >= min_points + 1 and np.ptp(x) > 0:
        # weighted least squares for diff = c2 + m*x (x is weighted-centered,
        # so the two parameters are orthogonal under these weights)
        m = np.sum(diff * x * w) / np.sum(x**2 * w)
        m_err = np.sqrt(1.0 / np.sum(x**2 * w))
        c2 = np.sum((diff - m * x) * w) / np.sum(w)
        chi2_s = float(np.sum((diff - c2 - m * x) ** 2 * w))
        res.slope = float(sign * m)
        res.slope_err = float(m_err)
        res.chi2_slope = chi2_s
        res.dof_slope = n - 2
        res.p_slope = float(stats.chi2.sf(chi2_s, n - 2))
    return res
