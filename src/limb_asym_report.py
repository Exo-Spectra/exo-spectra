"""Limb-asymmetry study — catalog statistics and summary report.

Inputs: reports/limb_asymmetry/limb_asym_catalog.csv (+ upper_limits.csv
if the injection stage has run). Produces:
  * BH-FDR 1% over the Method-A p-values
  * per-planet repeatability: inverse-variance weighted mean dd across
    visit-detectors + consistency chi2 (real limb asymmetry repeats
    across visits AND detectors; spots don't)
  * per-visit NRS1 vs NRS2 agreement
  * summary figure + reports/limb_asymmetry/summary.md

Usage: python limb_asym_report.py
Run after: limb_asym_run.py (and ideally limb_asym_inject.py)
"""
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_dist, norm

sys.path.insert(0, str(Path(__file__).parent))
from survey_analyze import ROOT                              # noqa: E402

OUT = ROOT / "reports" / "limb_asymmetry"


def bh_fdr(p, alpha=0.01):
    """Benjamini-Hochberg significance mask (same style as survey_analyze)."""
    p = np.asarray(p, float)
    m = np.isfinite(p).sum()
    thresh = 0.0
    for i, pv in enumerate(np.sort(p[np.isfinite(p)]), start=1):
        if pv <= alpha * i / m:
            thresh = pv
    return (p <= thresh) & np.isfinite(p), thresh


def weighted_mean(dd, err):
    w = 1.0 / err ** 2
    mean = np.sum(w * dd) / np.sum(w)
    mean_err = np.sqrt(1.0 / np.sum(w))
    chi2 = float(np.sum(w * (dd - mean) ** 2))
    dof = len(dd) - 1
    p = float(chi2_dist.sf(chi2, dof)) if dof > 0 else np.nan
    return mean, mean_err, chi2, dof, p


def main():
    cat = pd.read_csv(OUT / "limb_asym_catalog.csv")
    ok = cat[(cat.partial_flag == 0) & cat.dd_ppm.notna()
             & (cat.dd_err_ppm > 0)].copy()
    print(f"{len(cat)} catalog rows, {len(ok)} with an asymmetry fit, "
          f"{int(cat.partial_flag.sum())} partial transits")

    # calibrated errors from injection-recovery, when the stage has run:
    # sigma_cal (zero-injection scatter) supersedes the beta-inflated error
    # (injection shows beta double counts red noise the cyclic-shift
    # bootstrap already carries; sigma_cal/boot median ~0.9)
    ul_path = OUT / "upper_limits.csv"
    ul = None
    if ul_path.exists():
        ul = pd.read_csv(ul_path)
        key = ["target", "planet", "visit", "det"]
        ok = ok.merge(ul[key + ["sigma_cal_ppm", "sigma_cal_resid_ppm",
                                "ul95_abs_dd_ppm"]], on=key, how="left")
    ok["dd_err_cal_ppm"] = ok.get("sigma_cal_ppm", pd.Series(np.nan,
                                  index=ok.index)).fillna(ok.dd_err_ppm)
    ok["dd_sigma_cal"] = ok.dd_ppm / ok.dd_err_cal_ppm
    ok["p_value_cal"] = 2.0 * norm.sf(ok.dd_sigma_cal.abs())

    ok["significant"], thresh = bh_fdr(ok.p_value_cal)
    nsig = int(ok.significant.sum())
    print(f"BH-FDR 1% (calibrated errors): {nsig}/{len(ok)} significant "
          f"(p <= {thresh:.2e})")
    sig_beta, thr_beta = bh_fdr(ok.p_value)
    print(f"BH-FDR 1% (conservative beta-inflated): "
          f"{int(sig_beta.sum())}/{len(ok)} (p <= {thr_beta:.2e})")

    # ---------------------------------------------- per-planet repeatability
    lines = ["# Limb asymmetry — summary", "",
             f"Catalog: {len(ok)} visit-detector fits, "
             f"{int(cat.partial_flag.sum())} partials excluded, "
             f"{int(ok.spot_flag.sum())} spot-flagged.",
             f"Errors: injection-calibrated sigma_cal"
             f"{' (available)' if ul is not None else ' NOT available — '
                'falling back to beta-inflated'}.",
             f"Method A @ BH-FDR 1% (calibrated): {nsig} significant; "
             f"conservative beta-inflated: {int(sig_beta.sum())}.", "",
             "## Per-planet weighted mean Delta depth (egress - ingress)",
             "", "| planet | n fits | dd_w [ppm] | err [ppm] | consistency "
             "chi2/dof | p_cons | spot-flagged excluded |",
             "|---|---|---|---|---|---|---|"]
    rep_rows = []
    for planet, g in ok.groupby("planet"):
        gclean = g[g.spot_flag == 0]
        use = gclean if len(gclean) >= 2 else g
        mean, err, chi2, dof, p = weighted_mean(
            use.dd_ppm.to_numpy(), use.dd_err_cal_ppm.to_numpy())
        rep_rows.append({"planet": planet, "n": len(use), "dd_w_ppm": mean,
                         "dd_w_err_ppm": err, "chi2": chi2, "dof": dof,
                         "p_consistency": p,
                         "dd_w_sigma": mean / err if err > 0 else np.nan})
        lines.append(f"| {planet} | {len(use)} | {mean:+.0f} | {err:.0f} | "
                     f"{chi2:.1f}/{dof} | {p if np.isnan(p) else round(p,3)} |"
                     f" {len(g) - len(gclean)} |")
        print(f"{planet}: dd_w = {mean:+.0f} +- {err:.0f} ppm "
              f"({mean/err:+.1f} sigma), consistency chi2/dof {chi2:.1f}/{dof}"
              f" (p={p:.3f})" if dof > 0 else
              f"{planet}: dd_w = {mean:+.0f} +- {err:.0f} ppm (single fit)")
    rep = pd.DataFrame(rep_rows)
    rep.to_csv(OUT / "per_planet_summary.csv", index=False)

    # ------------------------------------------------ NRS1 vs NRS2 agreement
    lines += ["", "## NRS1 vs NRS2 per visit", "",
              "| planet | visit | dd_nrs1 | dd_nrs2 | diff/err |",
              "|---|---|---|---|---|"]
    for (planet, visit), g in ok.groupby(["planet", "visit"]):
        if set(g.det) != {"nrs1", "nrs2"}:
            continue
        a = g[g.det == "nrs1"].iloc[0]
        b = g[g.det == "nrs2"].iloc[0]
        z = (a.dd_ppm - b.dd_ppm) / np.sqrt(a.dd_err_cal_ppm ** 2
                                            + b.dd_err_cal_ppm ** 2)
        lines.append(f"| {planet} | {visit} | {a.dd_ppm:+.0f} | "
                     f"{b.dd_ppm:+.0f} | {z:+.1f} |")

    # ------------------------------------------------------ emcee escalation
    if "dd_mcmc_ppm" in ok.columns and ok.dd_mcmc_ppm.notna().any():
        lines += ["", "## emcee escalation (rows > 2 sigma pre-calibration)",
                  "", "| planet | visit | det | dd LSQ [ppm] | dd emcee [ppm]"
                  " | emcee sigma |", "|---|---|---|---|---|---|"]
        for _, r in ok[ok.dd_mcmc_ppm.notna()].iterrows():
            zmc = r.dd_mcmc_ppm / r.dd_mcmc_err_ppm
            lines.append(f"| {r.planet} | {r.visit} | {r.det} | "
                         f"{r.dd_ppm:+.0f} +- {r.dd_err_cal_ppm:.0f} | "
                         f"{r.dd_mcmc_ppm:+.0f} +- {r.dd_mcmc_err_ppm:.0f} | "
                         f"{zmc:+.1f} |")

    # ---------------------------------------------------------- upper limits
    if ul is not None:
        lines += ["", "## Sensitivity (injection-recovery)", "",
                  f"median sigma_cal = {ul.sigma_cal_ppm.median():.0f} ppm; "
                  f"median UL95 |dd| = {ul.ul95_abs_dd_ppm.median():.0f} ppm; "
                  f"calibration sigma_cal/boot = "
                  f"{ul.cal_vs_boot.median():.2f}, "
                  f"sigma_cal/beta-inflated = {ul.cal_vs_beta.median():.2f}; "
                  f"pooled-residual variant sigma_cal_resid/sigma_cal = "
                  f"{(ul.sigma_cal_resid_ppm / ul.sigma_cal_ppm).median():.2f}"
                  f"; median bias slope = {ul.bias_slope.median():.2f}",
                  ""]

    # -------------------------------------------------------- summary figure
    fig, ax = plt.subplots(figsize=(12, 6))
    planets = sorted(ok.planet.unique())
    xticks, xlabels = [], []
    x = 0
    for planet in planets:
        g = ok[ok.planet == planet].sort_values(["visit", "det"])
        xs = np.arange(len(g)) + x
        for det, col in (("nrs1", "C0"), ("nrs2", "C3")):
            sel = (g.det == det).to_numpy()
            if sel.any():
                ax.errorbar(xs[sel], g.dd_ppm[sel],
                            yerr=g.dd_err_cal_ppm[sel],
                            fmt="o", ms=5, color=col, capsize=2,
                            label=det if x == 0 else None)
        spot = (g.spot_flag == 1).to_numpy()
        if spot.any():
            ax.plot(xs[spot], g.dd_ppm[spot], "x", ms=10, color="k",
                    label="spot flag" if x == 0 else None)
        xticks.append(x + 0.5 * (len(g) - 1))
        xlabels.append(planet)
        x += len(g) + 1.5
    ax.axhline(0, color="0.7", lw=0.8)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=30, ha="right")
    ax.set_ylabel("Delta depth (egress - ingress) [ppm]")
    ax.set_title("Limb asymmetry across all visits — Method A")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "summary_dd.png", dpi=150)
    plt.close(fig)

    ok.to_csv(OUT / "limb_asym_catalog_fdr.csv", index=False)
    (OUT / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"-> {OUT / 'summary.md'}, summary_dd.png, per_planet_summary.csv")


if __name__ == "__main__":
    main()
