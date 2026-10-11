"""Gap #7 (retrieval robustness), stage 3: sensitivity tables and figures
from reports/retrieval_robustness/evidences_*.csv.

The L 98-59 b reduction pair is imported from study #1
(reports/l9859b_verification/evidences_*.csv, pub_firefly vs pub_eureka) --
those runs predate the post_mean/post_std columns, so the pair contributes
Delta(sigma) only, no parameter shifts. Its two sides were also retrieved on
their full (uncut) wavelength ranges -- caveat for the write-up.

Outputs (reports/retrieval_robustness/):
  sigma_all.csv          sigma(model vs flat) per (pair, side, code, model,
                         err_variant), with posterior mean/std
  sensitivity_pairs.csv  the headline table: per (pair, code, model, variant)
                         sigma_a, sigma_b, Delta(sigma) and parameter shifts
                         in units of the combined posterior sigma
  codes_agreement.csv    TauREx vs PLATON per (pair, side, model, variant)
  delta_sigma.png        heatmap of |Delta(sigma)| across the matrix
  fit__<pair>.png        both sides' data + best-fit primary-gas model
Run on a partial matrix at any time (missing cells are simply absent).
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
from gap7_retrieval import OUT, TaurexForward, build_cells  # noqa: E402
from l9859b_retrieval import bin_model, lnb_to_sigma  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]

SIGMA_CAP = 6.99  # solver brentq (xtol~2e-12) i tak saturuje przy p~1e-12


def lnb_to_sigma_capped(lnb):
    """lnb_to_sigma z jawnym capem: przy lnB>~683 brentq pada i bazowa
    funkcja zwraca inf — dla raportu to to samo 'poza skala' co cap 6.99."""
    s = lnb_to_sigma(lnb)
    if not np.isfinite(s):
        return SIGMA_CAP if lnb > 0 else 0.0
    return min(s, SIGMA_CAP)
L9859B = ROOT / "reports" / "l9859b_verification"

PARAM_NAMES = {3: ["rp_re", "T", "logx"], 4: ["rp_re", "T", "logx", "logpc"]}


def load_evidence():
    # tylko main CSV per kod; glob("evidences_*.csv") lapalby tez pliki
    # shardow (evidences_<code>_s0.csv...), dublujac kazda komorke
    paths = [p for p in sorted(OUT.glob("evidences_*.csv"))
             if "_s" not in p.stem.removeprefix("evidences_")]
    frames = [pd.read_csv(p) for p in paths]
    ev = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if not ev.empty:
        ev = ev.drop_duplicates(
            subset=["spectrum", "code", "model", "err_variant"], keep="last")
    l1 = [pd.read_csv(p) for p in sorted(L9859B.glob("evidences_*.csv"))]
    if l1:
        l1 = pd.concat(l1, ignore_index=True)
        l1 = l1[l1.spectrum.isin(["pub_firefly", "pub_eureka"])].copy()
        l1["side"] = l1.spectrum.map({"pub_firefly": "a", "pub_eureka": "b"})
        l1["pair"], l1["planet"] = "l9859b", "L 98-59 b"
        l1["axis"], l1["err_variant"] = "reduction", "raw"
        l1["spectrum"] = "l9859b__" + l1["side"]
        ev = pd.concat([ev, l1], ignore_index=True)
    for col in ("post_mean", "post_std"):
        if col not in ev.columns:
            ev[col] = np.nan
    return ev


def sigma_table(ev):
    rows = []
    for (pair, side, code, var), g in ev.groupby(
            ["pair", "side", "code", "err_variant"]):
        flat = g[g.model == "flat"]
        if flat.empty:
            continue
        lnz0 = float(flat.lnZ.iloc[0])
        for _, r in g[g.model != "flat"].iterrows():
            lnb = float(r.lnZ) - lnz0
            rows.append({
                "pair": pair, "planet": r.planet, "axis": r.axis,
                "side": side, "code": code, "model": r.model,
                "err_variant": var, "lnB_vs_flat": round(lnb, 2),
                "sigma": round(lnb_to_sigma_capped(lnb), 2),
                "best_params": r.best_params,
                "post_mean": r.post_mean, "post_std": r.post_std,
            })
    return pd.DataFrame(rows)


def pair_table(tab):
    """The headline metric: how much does the conclusion move when you swap
    the reduction / epoch / pipeline / instrument (side a -> side b)?"""
    a = tab[tab.side == "a"].set_index(["pair", "code", "model",
                                        "err_variant"])
    b = tab[tab.side == "b"].set_index(["pair", "code", "model",
                                        "err_variant"])
    rows = []
    for key in a.index.intersection(b.index):
        ra, rb = a.loc[key], b.loc[key]
        row = {
            "pair": key[0], "planet": ra.planet, "axis": ra.axis,
            "code": key[1], "model": key[2], "err_variant": key[3],
            "sigma_a": ra.sigma, "sigma_b": rb.sigma,
            "dsigma": round(rb.sigma - ra.sigma, 2),
        }
        try:  # parameter shifts in units of the combined posterior sigma
            ma, sa = json.loads(ra.post_mean), json.loads(ra.post_std)
            mb, sb = json.loads(rb.post_mean), json.loads(rb.post_std)
            names = PARAM_NAMES[len(ma)]
            for i, name in enumerate(names):
                denom = np.hypot(sa[i], sb[i])
                row[f"z_{name}"] = (round(abs(ma[i] - mb[i]) / denom, 2)
                                    if denom > 0 else np.nan)
        except (TypeError, ValueError):
            pass  # study #1 import has no posterior columns
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["pair", "model", "err_variant",
                                           "code"])


def codes_table(tab):
    t = tab[tab.code == "taurex"].set_index(["pair", "side", "model",
                                             "err_variant"])
    p = tab[tab.code == "platon"].set_index(["pair", "side", "model",
                                             "err_variant"])
    rows = []
    for key in t.index.intersection(p.index):
        rows.append({
            "pair": key[0], "side": key[1], "model": key[2],
            "err_variant": key[3], "sigma_taurex": t.loc[key].sigma,
            "sigma_platon": p.loc[key].sigma,
            "dsigma_codes": round(p.loc[key].sigma - t.loc[key].sigma, 2),
            "dlnB_codes": round(p.loc[key].lnB_vs_flat
                                - t.loc[key].lnB_vs_flat, 2),
        })
    return pd.DataFrame(rows)


def heatmap(pt):
    if pt.empty:
        return
    pt = pt.copy()
    pt["row"] = pt.pair + np.where(pt.err_variant == "infl", " (infl)", "")
    pt["col"] = pt.model + "\n" + pt.code
    rows = sorted(pt.row.unique())
    cols = sorted(pt.col.unique())
    grid = np.full((len(rows), len(cols)), np.nan)
    for _, r in pt.iterrows():
        grid[rows.index(r.row), cols.index(r.col)] = abs(r.dsigma)
    fig, ax = plt.subplots(figsize=(1.1 * len(cols) + 3,
                                    0.55 * len(rows) + 2.5))
    im = ax.imshow(grid, cmap="RdYlGn_r", vmin=0, vmax=3, aspect="auto")
    ax.set_xticks(range(len(cols)), cols, fontsize=8)
    ax.set_yticks(range(len(rows)), rows, fontsize=9)
    for _, r in pt.iterrows():
        ax.text(cols.index(r.col), rows.index(r.row),
                f"{r.sigma_a:.1f}→{r.sigma_b:.1f}",
                ha="center", va="center", fontsize=7)
    ax.set_title("Gap #7: |Delta sigma| detekcji przy zamianie strony pary "
                 "(a→b)\nkolor = |Delta sigma|, tekst = sigma_a→"
                 "sigma_b")
    fig.colorbar(im, label="|Delta sigma|")
    fig.tight_layout()
    fig.savefig(OUT / "delta_sigma.png", dpi=150)
    plt.close(fig)


def fit_plots(ev):
    cells = {(c["pair"], c["side"]): c for c in build_cells()}
    for pair in sorted({p for p, _ in cells}):
        fig, ax = plt.subplots(figsize=(9, 5))
        plotted = False
        for side, color in (("a", "tab:blue"), ("b", "tab:red")):
            c = cells[(pair, side)]
            gas = c["models"][1].upper()   # primary gas (index 0 = flat)
            ax.errorbar(c["wave"], c["depth"], yerr=c["err"],
                        xerr=c["dwave"] / 2, fmt="o", ms=3, lw=0.8,
                        color=color, alpha=0.55,
                        label=f"strona {side} (dane)")
            g = ev[(ev.spectrum == c["spectrum"]) & (ev.code == "taurex")
                   & (ev.model == gas.lower())
                   & (ev.err_variant == "raw")]
            if g.empty:
                continue
            fwd = TaurexForward(gas, c)
            model = fwd.depths_ppm(*json.loads(g.best_params.iloc[0]))
            binned = bin_model(fwd.wave_model, model, c["wave"], c["dwave"])
            ax.plot(c["wave"], binned, color=color, lw=1.5,
                    label=f"strona {side}: best {gas} (taurex)")
            plotted = True
        if not plotted:
            plt.close(fig)
            continue
        c = cells[(pair, "a")]
        ax.set_xlabel("długość fali [µm]")
        ax.set_ylabel("głębokość tranzytu [ppm]")
        ax.set_title(f"{c['planet']} — para {pair} (oś: {c['axis']})")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT / f"fit__{pair}.png", dpi=150)
        plt.close(fig)


def main():
    ev = load_evidence()
    if ev.empty:
        print("brak evidences_*.csv")
        return
    tab = sigma_table(ev)
    tab.to_csv(OUT / "sigma_all.csv", index=False)
    pt = pair_table(tab)
    pt.to_csv(OUT / "sensitivity_pairs.csv", index=False)
    ct = codes_table(tab)
    ct.to_csv(OUT / "codes_agreement.csv", index=False)

    n_total = len(ev)
    print(f"wierszy evidence: {n_total} | komorek sigma: {len(tab)} "
          f"| par policzonych obustronnie: {len(pt)}\n")
    if not pt.empty:
        show = [c for c in ("pair", "code", "model", "err_variant", "sigma_a",
                            "sigma_b", "dsigma", "z_rp_re", "z_T", "z_logx",
                            "z_logpc") if c in pt.columns]
        print(pt[show].to_string(index=False))
    if not ct.empty:
        bad = ct[ct.dlnB_codes.abs() > 1]
        print(f"\nzgodnosc kodow: {len(ct) - len(bad)}/{len(ct)} komorek "
              f"z |dlnB| <= 1" + ("" if bad.empty else
              f"; ROZJAZDY:\n{bad.to_string(index=False)}"))
    heatmap(pt)
    fit_plots(ev)
    print(f"\nzapisano: sigma_all / sensitivity_pairs / codes_agreement .csv"
          f" + delta_sigma.png + fit__*.png w {OUT}")


if __name__ == "__main__":
    main()
