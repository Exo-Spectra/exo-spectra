"""Phase 5 (Tier B): shape-anomaly detection within homogeneous instrument
cohorts ("oddball hunting", cf. Matchev+ 2022, arXiv:2601.02324 — both on
synthetic spectra; here applied to the real archive for the first time).

Per cohort (instrument regex + wavelength window):
  1. resample every member onto a fixed common grid (inverse-variance-weighted
     bin means); membership requires coverage of >= MIN_COVER of the bins,
     remaining gaps are imputed with the cohort median (per bin, flagged);
  2. normalize each spectrum to a pure *shape*: subtract the weighted mean,
     divide by the std of the binned values (removes (Rp/Rs)^2 scale and
     scale-height amplitude — model-independent);
  3. PCA via SVD on the mean-centered shape matrix; k = components explaining
     >= PCA_VAR of variance (capped at PCA_KMAX);
  4. three outlier scores per spectrum: reconstruction residual outside the
     k-PC subspace, robust Mahalanobis distance in PC space (iteratively
     trimmed center/covariance), and mean kNN distance in PC space.
     Final rank = median of the three ranks.

`amp_snr` (shape amplitude / median bin error) is carried along so that
low-SNR spectra — whose "shape" is mostly noise — can be recognized when
interpreting the ranking.

Outputs: data/processed/phase5_cohort_scores.csv,
         reports/phase5_pca_<cohort>.png, reports/phase5_oddball_<cohort>_<rank>.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).parent))
from compare import _sym_sigma
from spectra_io import load_spectrum, usable_points

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed"
REPORTS = ROOT / "reports"

MIN_MEMBERS = 15
MIN_COVER = 0.8
PCA_VAR = 0.90
PCA_KMAX = 8
KNN_K = 5
TRIM_FRAC = 0.10   # robust Mahalanobis: fraction trimmed per iteration
TOP_PLOTS = 5

# name, spec_type, instrument regex (case-insensitive), exclude regex, lo, hi, n_bins
COHORTS = [
    ("T_WFC3_G141",   "Transmission", r"WFC3|Wide Field",  None,      1.10, 1.66, 18),
    ("T_STIS",        "Transmission", r"Imaging Spectrograph", None,  0.30, 1.02, 12),
    ("T_G395",        "Transmission", r"NIRSpec",          r"PRISM",  2.87, 5.10, 20),
    ("T_PRISM",       "Transmission", r"PRISM",            None,      0.60, 5.20, 20),
    ("T_NIRISS",      "Transmission", r"NIRISS",           None,      0.65, 2.80, 20),
    ("T_NIRCam",      "Transmission", r"NIRCam",           None,      2.45, 3.95, 16),
    ("T_MIRI_LRS",    "Transmission", r"MIRI",             None,      5.00, 12.0, 14),
    ("E_WFC3_G141",   "Eclipse",      r"WFC3|Wide Field",  None,      1.10, 1.66, 18),
    ("E_NIRCam",      "Eclipse",      r"NIRCam",           None,      2.45, 3.95, 16),
]


def resample(spec: pd.DataFrame, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Inverse-variance-weighted mean per bin. Returns (values, sigmas), NaN where empty."""
    sig = _sym_sigma(spec)
    ok = np.isfinite(sig) & (sig > 0)
    wave = spec["wave"].to_numpy(float)[ok]
    val = spec["value"].to_numpy(float)[ok]
    sig = sig[ok]
    nb = len(edges) - 1
    v = np.full(nb, np.nan)
    s = np.full(nb, np.nan)
    idx = np.digitize(wave, edges) - 1
    for b in range(nb):
        m = idx == b
        if m.any():
            w = 1.0 / sig[m] ** 2
            v[b] = np.sum(val[m] * w) / np.sum(w)
            s[b] = np.sqrt(1.0 / np.sum(w))
    return v, s


def robust_mahalanobis(pcs: np.ndarray) -> np.ndarray:
    """Mahalanobis distance with center/covariance from iteratively trimmed data."""
    keep = np.ones(len(pcs), bool)
    d = np.zeros(len(pcs))
    for _ in range(3):
        mu = pcs[keep].mean(axis=0)
        cov = np.cov(pcs[keep].T) + 1e-9 * np.eye(pcs.shape[1])
        inv = np.linalg.inv(cov)
        diff = pcs - mu
        d = np.sqrt(np.einsum("ij,jk,ik->i", diff, inv, diff))
        cut = np.quantile(d, 1 - TRIM_FRAC)
        keep = d <= cut
    return d


def analyze_cohort(name, spec_type, lo, hi, n_bins, members, specs) -> pd.DataFrame | None:
    edges = np.linspace(lo, hi, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    rows, mat, sigs = [], [], []
    for sid, meta in members.iterrows():
        v, s = resample(specs[sid], edges)
        cover = np.isfinite(v).mean()
        if cover < MIN_COVER:
            continue
        rows.append({"spec_id": sid, "pl_name": meta.pl_name, "authors": meta.authors,
                     "instrument": meta.instrument, "bibcode": meta.bibcode,
                     "coverage": cover})
        mat.append(v)
        sigs.append(s)
    if len(rows) < MIN_MEMBERS:
        print(f"  {name}: only {len(rows)} members with coverage >= {MIN_COVER} — SKIP")
        return None

    X = np.array(mat)
    S = np.array(sigs)
    info = pd.DataFrame(rows)

    # shape normalization: weighted mean out, unit std
    w = np.where(np.isfinite(S), 1.0 / S**2, 0.0)
    vfill = np.where(np.isfinite(X), X, 0.0)
    mean = np.sum(vfill * w, axis=1) / np.sum(w, axis=1)
    Z = X - mean[:, None]
    std = np.nanstd(Z, axis=1)
    std[std == 0] = 1.0
    Z = Z / std[:, None]
    info["amp_snr"] = std / np.nanmedian(S, axis=1)

    # impute missing bins with the cohort median shape
    med_shape = np.nanmedian(Z, axis=0)
    miss = ~np.isfinite(Z)
    Z[miss] = np.broadcast_to(med_shape, Z.shape)[miss]

    # PCA via SVD on column-centered shapes
    col_mean = Z.mean(axis=0)
    Zc = Z - col_mean
    U, sv, Vt = np.linalg.svd(Zc, full_matrices=False)
    var = sv**2 / np.sum(sv**2)
    k = min(int(np.searchsorted(np.cumsum(var), PCA_VAR)) + 1, PCA_KMAX, len(sv))
    pcs = U[:, :k] * sv[:k]

    # scores
    recon = Zc @ Vt[:k].T @ Vt[:k]
    info["score_recon"] = np.sum((Zc - recon) ** 2, axis=1)
    info["score_mahal"] = robust_mahalanobis(pcs)
    tree = cKDTree(pcs)
    dist, _ = tree.query(pcs, k=min(KNN_K + 1, len(pcs)))
    info["score_knn"] = dist[:, 1:].mean(axis=1)
    ranks = np.vstack([info[c].rank(ascending=False) for c in
                       ("score_recon", "score_mahal", "score_knn")])
    info["oddball_rank"] = np.median(ranks, axis=0)
    info["cohort"] = name
    info["pca_k"] = k
    info = info.sort_values("oddball_rank").reset_index(drop=True)
    print(f"  {name}: {len(info)} spectra, {info.pl_name.nunique()} planets, "
          f"k={k} PCs ({np.cumsum(var)[k-1]:.0%} var); top oddball: "
          f"{info.pl_name.iloc[0]} ({info.authors.iloc[0]})")

    # --- plots ----------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(pcs[:, 0], pcs[:, 1] if k > 1 else np.zeros(len(pcs)), s=18, alpha=0.6)
    for i in range(min(TOP_PLOTS, len(info))):
        j = info.index[i]
        row = info.loc[j]
        # info was sorted; map back to matrix row via spec_id order
        mrow = np.nonzero(np.array([r["spec_id"] for r in rows]) == row.spec_id)[0][0]
        ax.scatter(pcs[mrow, 0], pcs[mrow, 1] if k > 1 else 0, s=60, color="C3")
        ax.annotate(row.pl_name, (pcs[mrow, 0], pcs[mrow, 1] if k > 1 else 0),
                    fontsize=7, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    ax.set_title(f"{name}: shape space ({len(info)} spectra), top {TOP_PLOTS} oddballs")
    fig.tight_layout()
    fig.savefig(REPORTS / f"phase5_pca_{name}.png", dpi=150)
    plt.close(fig)

    band_lo = np.nanquantile(Z, 0.25, axis=0)
    band_hi = np.nanquantile(Z, 0.75, axis=0)
    for i in range(min(TOP_PLOTS, len(info))):
        row = info.iloc[i]
        mrow = np.nonzero(np.array([r["spec_id"] for r in rows]) == row.spec_id)[0][0]
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.fill_between(centers, band_lo, band_hi, alpha=0.25, label="cohort IQR")
        ax.plot(centers, med_shape, lw=1.2, color="k", label="cohort median shape")
        ax.plot(centers, Z[mrow], "o-", color="C3", ms=4,
                label=f"{row.pl_name} — {row.authors}")
        ax.set_xlabel("wavelength [microns]")
        ax.set_ylabel("normalized shape (mean 0, std 1)")
        ax.set_title(f"{name} oddball #{i+1}: {row.pl_name} "
                     f"(recon {row.score_recon:.1f}, mahal {row.score_mahal:.1f}, "
                     f"amp_snr {row.amp_snr:.1f})")
        ax.legend(fontsize=8)
        fig.tight_layout()
        slug = row.pl_name.replace(" ", "_")
        fig.savefig(REPORTS / f"phase5_oddball_{name}_{i+1:02d}_{slug}.png", dpi=150)
        plt.close(fig)
    return info


def main() -> None:
    summary = pd.read_csv(OUT / "spectra_summary.csv").set_index("spec_id")
    REPORTS.mkdir(exist_ok=True)
    results = []
    for name, st, pat, excl, lo, hi, nb in COHORTS:
        m = ((summary.spec_type == st)
             & summary.instrument.str.contains(pat, case=False, na=False, regex=True)
             & (summary.n_usable >= 5)
             & (summary.wave_min <= lo + 0.3 * (hi - lo))
             & (summary.wave_max >= hi - 0.3 * (hi - lo)))
        if excl:
            m &= ~summary.instrument.str.contains(excl, case=False, na=False, regex=True)
        members = summary[m]
        specs = {}
        for sid, r in members.iterrows():
            specs[sid] = usable_points(load_spectrum(ROOT / "data" / "spectra" / r.file, r.spec_type))
        res = analyze_cohort(name, st, lo, hi, nb, members, specs)
        if res is not None:
            results.append(res)
    if results:
        all_scores = pd.concat(results, ignore_index=True)
        all_scores.to_csv(OUT / "phase5_cohort_scores.csv", index=False)
        print(f"\n{len(all_scores)} spectra scored in {len(results)} cohorts "
              f"-> data/processed/phase5_cohort_scores.csv")


if __name__ == "__main__":
    main()
