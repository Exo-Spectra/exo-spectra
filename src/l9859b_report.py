"""Study 4, step 3: build the results matrix and figures from
reports/l9859b_verification/evidences_{taurex,platon}.csv.

Outputs (reports/l9859b_verification/):
  sigma_matrix.csv    sigma(model vs flat) per (spectrum, code, model)
  sigma_matrix.png    heatmap of the SO2/CO2-vs-flat sigmas
  fit_<spectrum>.png  data + best-fit flat/SO2/CO2 models (TauREx curves)
Run after both retrieval matrices are complete.
"""
import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from l9859b_retrieval import (SPECTRA_ALL, TaurexForward, bin_model,
                              lnb_to_sigma, load_data)  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "l9859b_verification"

SPECTRA = SPECTRA_ALL
LABELS = {
    "pub_firefly": "published FIREFLy (80 pts)",
    "pub_eureka": "published Eureka! (218 pts)",
    "own_avg": "own mean of 4 visits",
    "own_avg_infl": "own mean, errors ×√χ²",
    "own_v1": "own visit 1", "own_v2": "own visit 2",
    "own_v3": "own visit 3", "own_v4": "own visit 4",
    "synth_flat": "flat synthetic (control)",
}


def sigma_table(ev):
    rows = []
    for (spec, code), g in ev.groupby(["spectrum", "code"]):
        flat = g[g["model"] == "flat"]
        if flat.empty:
            continue
        lnz0 = float(flat["lnZ"].iloc[0])
        for _, r in g[g["model"] != "flat"].iterrows():
            lnb = float(r["lnZ"]) - lnz0
            rows.append({
                "spectrum": spec, "code": code, "model": r["model"],
                "lnB_vs_flat": round(lnb, 2),
                "sigma": round(lnb_to_sigma(lnb), 2),
                "best_params": r["best_params"],
            })
    return pd.DataFrame(rows)


def heatmap(tab):
    specs = [s for s in SPECTRA if s in set(tab["spectrum"])]
    codes = ["taurex", "platon"]
    models = ["so2", "co2"]
    cols = [(c, m) for m in models for c in codes]
    grid = np.full((len(specs), len(cols)), np.nan)
    for i, s in enumerate(specs):
        for j, (c, m) in enumerate(cols):
            r = tab[(tab.spectrum == s) & (tab.code == c) & (tab.model == m)]
            if len(r):
                grid[i, j] = r["sigma"].iloc[0]
    fig, ax = plt.subplots(figsize=(8, 0.6 * len(specs) + 2))
    im = ax.imshow(grid, cmap="RdYlGn_r", vmin=0, vmax=4, aspect="auto")
    ax.set_xticks(range(len(cols)),
                  [f"{m.upper()}\n{c}" for c, m in cols])
    ax.set_yticks(range(len(specs)), [LABELS.get(s, s) for s in specs])
    for i in range(len(specs)):
        for j in range(len(cols)):
            if np.isfinite(grid[i, j]):
                ax.text(j, i, f"{grid[i, j]:.1f}", ha="center", va="center",
                        fontsize=9)
    ax.set_title("L 98-59 b: model preference vs flat line [sigma]\n"
                 "(Bello-Arufe 2025: SO2 at 2.2-3.6 sigma depending on code)")
    fig.colorbar(im, label="sigma")
    fig.tight_layout()
    fig.savefig(OUT / "sigma_matrix.png", dpi=150)
    plt.close(fig)


def fit_plots(ev):
    fwd = {}
    for spec in SPECTRA:
        g = ev[(ev.spectrum == spec) & (ev.code == "taurex")]
        if g.empty:
            continue
        wave, dwave, depth, err = load_data(spec)
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.errorbar(wave, depth, yerr=err, xerr=dwave / 2, fmt="o", ms=3,
                    lw=0.8, color="k", alpha=0.65, label="data")
        colors = {"flat": "gray", "so2": "tab:red", "co2": "tab:blue"}
        for _, r in g.iterrows():
            p = json.loads(r["best_params"])
            if r["model"] == "flat":
                ax.axhline(p[0], color=colors["flat"], ls="--", lw=1,
                           label="flat (best)")
                continue
            gas = r["model"].upper()
            if gas not in fwd:
                fwd[gas] = TaurexForward(gas)
            model = fwd[gas].depths_ppm(*p)
            binned = bin_model(fwd[gas].wave_model, model, wave, dwave)
            ax.plot(wave, binned, color=colors[r["model"]], lw=1.5,
                    label=f"{gas} (best, taurex)")
        ax.set_xlabel("wavelength [µm]")
        ax.set_ylabel("transit depth [ppm]")
        ax.set_title(f"L 98-59 b — {LABELS.get(spec, spec)}")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT / f"fit_{spec}.png", dpi=150)
        plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ev = pd.concat([pd.read_csv(p) for p in OUT.glob("evidences_*.csv")],
                   ignore_index=True)
    tab = sigma_table(ev)
    tab.to_csv(OUT / "sigma_matrix.csv", index=False)
    pt = tab.pivot_table(index="spectrum", columns=["model", "code"],
                         values="sigma")
    pt = pt.reindex([s for s in SPECTRA if s in pt.index])
    print(pt.to_string())
    heatmap(tab)
    fit_plots(ev)
    print(f"\nwritten: sigma_matrix.csv/png + fit_*.png in {OUT}")


if __name__ == "__main__":
    main()
