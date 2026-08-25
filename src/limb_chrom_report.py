"""Chromatic limb-asymmetry study (gap #8), stage 3b: statistics + report.

Works on the (possibly partial) calibration table. All significances use
the injection-calibrated sigma_cal (study-05 rule).

Per calibrated bin: z = dd / sigma_cal; BH-FDR at 1% across all bins.
Per (visit, det, planet): chromatic-structure chi2 of dd(lambda) against
its own weighted mean (does the asymmetry DEPEND on wavelength within one
transit?).
Per planet: inverse-variance stacked dd per bin across visits/detectors
(bins aligned by detector+index), molecular-band vs continuum contrast
(CH4 3.20-3.45, CO2 4.20-4.45, CO 4.50-4.75 um), and an upper-limit table.

-> reports/limb_chrom/summary.md + limb_chrom_fdr.csv + stack_planet_bin.csv
   + band_contrast.csv + dd_lambda__<planet>.png

Usage: python limb_chrom_report.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402
import pandas as pd               # noqa: E402
from scipy import stats           # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from survey_analyze import ROOT   # noqa: E402

CHROM = ROOT / "reports" / "limb_chrom"
BANDS = {"CH4": (3.20, 3.45), "CO2": (4.20, 4.45), "CO": (4.50, 4.75)}
FDR_Q = 0.01


def bh_fdr(p, q=FDR_Q):
    p = np.asarray(p)
    n = len(p)
    order = np.argsort(p)
    passed = np.zeros(n, bool)
    thresh = q * (np.arange(1, n + 1)) / n
    ok = p[order] <= thresh
    if ok.any():
        passed[order[:np.max(np.where(ok)[0]) + 1]] = True
    return passed


def wmean(x, e):
    w = 1.0 / np.asarray(e, float) ** 2
    m = np.sum(w * x) / np.sum(w)
    return m, np.sqrt(1.0 / np.sum(w))


def band_of(lo, hi):
    c = 0.5 * (lo + hi)
    for name, (blo, bhi) in BANDS.items():
        if blo <= c < bhi:
            return name
    return "cont"


def main():
    cat = pd.read_csv(CHROM / "limb_chrom_catalog.csv")
    cal = pd.read_csv(CHROM / "limb_chrom_calibration.csv")
    d = pd.merge(cat, cal[["visit", "det", "planet", "bin", "sigma_cal_ppm",
                           "sigma_resid_ppm", "bias_slope"]],
                 on=["visit", "det", "planet", "bin"])
    d = d[np.isfinite(d.sigma_cal_ppm) & (d.sigma_cal_ppm > 0)].copy()
    d["z_cal"] = d.dd_ppm / d.sigma_cal_ppm
    d["p_cal"] = 2 * stats.norm.sf(np.abs(d.z_cal))
    d["fdr_pass"] = bh_fdr(d.p_cal.to_numpy())
    d["band"] = [band_of(a, b) for a, b in zip(d.wave_lo_um, d.wave_hi_um)]
    d["ul95_ppm"] = np.abs(d.dd_ppm) + 1.645 * d.sigma_cal_ppm
    d.to_csv(CHROM / "limb_chrom_fdr.csv", index=False)

    # chromatic structure per fit: dd(lambda) vs its own weighted mean
    struct = []
    for (v, det, pl), g in d.groupby(["visit", "det", "planet"]):
        if len(g) < 4:
            continue
        m, _ = wmean(g.dd_ppm, g.sigma_cal_ppm)
        chi2 = float(np.sum(((g.dd_ppm - m) / g.sigma_cal_ppm) ** 2))
        dof = len(g) - 1
        struct.append({"visit": v, "det": det, "planet": pl, "n_bins": len(g),
                       "dd_mean_ppm": m, "chi2": chi2, "dof": dof,
                       "p_struct": float(stats.chi2.sf(chi2, dof))})
    struct = pd.DataFrame(struct)
    struct["fdr_pass"] = bh_fdr(struct.p_struct.to_numpy())
    struct.to_csv(CHROM / "chromatic_structure.csv", index=False)

    # per-planet stack: weighted mean dd per (det, bin) across visits
    stacks = []
    for (pl, det, b), g in d.groupby(["planet", "det", "bin"]):
        m, e = wmean(g.dd_ppm, g.sigma_cal_ppm)
        stacks.append({"planet": pl, "det": det, "bin": b,
                       "wave_lo_um": g.wave_lo_um.iloc[0],
                       "wave_hi_um": g.wave_hi_um.iloc[0],
                       "band": g.band.iloc[0], "n_fits": len(g),
                       "dd_ppm": m, "err_ppm": e, "z": m / e})
    stacks = pd.DataFrame(stacks)
    stacks["p"] = 2 * stats.norm.sf(np.abs(stacks.z))
    stacks["fdr_pass"] = bh_fdr(stacks.p.to_numpy())
    stacks.to_csv(CHROM / "stack_planet_bin.csv", index=False)

    # band-vs-continuum contrast per planet
    rows = []
    for pl, g in d.groupby("planet"):
        cont = g[g.band == "cont"]
        if not len(cont):
            continue
        mc, ec = wmean(cont.dd_ppm, cont.sigma_cal_ppm)
        for band in BANDS:
            gb = g[g.band == band]
            if not len(gb):
                continue
            mb, eb = wmean(gb.dd_ppm, gb.sigma_cal_ppm)
            dz = (mb - mc) / np.hypot(eb, ec)
            rows.append({"planet": pl, "band": band, "n_bins": len(gb),
                         "dd_band_ppm": mb, "err_band_ppm": eb,
                         "dd_cont_ppm": mc, "err_cont_ppm": ec,
                         "contrast_ppm": mb - mc,
                         "contrast_err_ppm": float(np.hypot(eb, ec)),
                         "contrast_z": float(dz),
                         "p": float(2 * stats.norm.sf(abs(dz)))})
    bands = pd.DataFrame(rows)
    if len(bands):
        bands["fdr_pass"] = bh_fdr(bands.p.to_numpy())
    bands.to_csv(CHROM / "band_contrast.csv", index=False)

    # dd(lambda) plot per planet
    for pl, g in stacks.groupby("planet"):
        fig, ax = plt.subplots(figsize=(9, 4.5))
        wc = 0.5 * (g.wave_lo_um + g.wave_hi_um)
        ax.errorbar(wc, g.dd_ppm, yerr=g.err_ppm, fmt="o", ms=4, capsize=2)
        ax.axhline(0, color="0.6", lw=0.8)
        for name, (blo, bhi) in BANDS.items():
            ax.axvspan(blo, bhi, alpha=0.10, color="C2")
            ax.text(0.5 * (blo + bhi), ax.get_ylim()[1] * 0.9, name,
                    ha="center", fontsize=8, color="C2")
        ax.set_xlabel("wavelength [um]")
        ax.set_ylabel("stacked dd = egress - ingress [ppm]")
        ax.set_title(f"{pl} - chromatic limb-asymmetry contrast "
                     f"(divide-white)")
        fig.tight_layout()
        fig.savefig(CHROM / f"dd_lambda__{pl.replace(' ', '_')}.png", dpi=140)
        plt.close(fig)

    with open(CHROM / "summary.md", "w", encoding="utf-8") as fh:
        fh.write("# Chromatic limb asymmetry - summary\n\n")
        fh.write(f"Calibrated bins: {len(d)} (of {len(cat)} fitted); "
                 f"median sigma_cal = {d.sigma_cal_ppm.median():.0f} ppm; "
                 f"median UL95 = {d.ul95_ppm.median():.0f} ppm\n\n")
        fh.write(f"Per-bin FDR {FDR_Q:.0%}: "
                 f"{int(d.fdr_pass.sum())} significant of {len(d)}\n\n")
        fh.write(f"Chromatic structure (dd(lambda) vs const per fit) "
                 f"@FDR {FDR_Q:.0%}: {int(struct.fdr_pass.sum())} of "
                 f"{len(struct)}\n\n")
        fh.write(f"Planet-stacked bins @FDR {FDR_Q:.0%}: "
                 f"{int(stacks.fdr_pass.sum())} of {len(stacks)}\n\n")
        if len(bands):
            fh.write(f"Band-vs-continuum contrasts @FDR {FDR_Q:.0%}: "
                     f"{int(bands.fdr_pass.sum())} of {len(bands)}\n\n")
            fh.write("## Band contrasts\n\n")
            fh.write(bands.round(1).to_string(index=False))
            fh.write("\n\n")
        fh.write("## Most significant stacked bins\n\n")
        fh.write(stacks.nlargest(10, "z", keep="all")[
            ["planet", "det", "bin", "wave_lo_um", "wave_hi_um", "band",
             "n_fits", "dd_ppm", "err_ppm", "z"]].round(2).to_string(
                 index=False))
        fh.write("\n")
    print(open(CHROM / "summary.md", encoding="utf-8").read())


if __name__ == "__main__":
    main()
