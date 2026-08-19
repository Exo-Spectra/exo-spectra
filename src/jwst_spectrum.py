"""Own transmission spectrum of GJ 1132 b from JWST rateints:
per-wavelength-bin photometry -> transit depth per bin -> comparison with
the published May et al. 2023 spectrum from the NASA Exoplanet Archive.

Wavelength solution: taken from the official x1dints product (per-column
wavelengths); photometry, background, detrending and depths are our own.

Usage: python jwst_spectrum.py [visit]   (default: 1)
Transit windows per visit come from the matched-filter detection in
jwst_whitelight.py (identical extraction code for every visit).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).parent))
from spectra_io import load_spectrum, local_name, usable_points
from compare import compare_pair

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "reports"
APER_HALF = 6
BG_MARGIN = 10
NBINS = 22          # per detector -> R ~ 100-ish

# per-visit: transit window [hours since start] from jwst_whitelight.py
# (matched box filter), and the archive `note` of the published counterpart
VISITS = {
    "1": dict(t_in=(1.65, 2.25), t_out_gap=(1.45, 2.40),
              pub_note=r"Eureka R100 .Visit 1."),
    "2": dict(t_in=(1.61, 2.21), t_out_gap=(1.43, 2.38),
              pub_note=r"Eureka R100 .Visit 2."),
}


import re


def column_wavelengths(x1d_path: Path, col_profile: np.ndarray) -> np.ndarray:
    """Per-column wavelength [um] from the x1dints table.

    The extract_2d step cuts the slit out of the full detector; the cutout's
    x-extent is not in any FITS keyword, but the pipeline log embedded in the
    ASDF extension states it verbatim ("Subarray x-extents are: A B").
    Wavelengths are identical across integrations (verified).
    """
    nx = len(col_profile)
    with fits.open(x1d_path) as hdul:
        tab = hdul["EXTRACT1D"].data
        wl = np.asarray(tab["WAVELENGTH"], float)[0]
        txt = hdul["ASDF"].data.tobytes().decode("utf-8", errors="ignore")
    m = re.search(r"Subarray\s*\\?n?\s*x-extents are:\s*\\?n?\s*(\d+)\s+(\d+)", txt)
    if not m:
        raise RuntimeError(f"x-extents not found in ASDF log of {x1d_path.name}")
    x0, x1 = int(m.group(1)), int(m.group(2))
    if x1 - x0 != len(wl):
        raise RuntimeError(f"extent {x0}..{x1} != {len(wl)} wavelength points")
    print(f"  slit cutout columns {x0}..{x1} (from pipeline log)")
    wave = np.full(nx, np.nan)
    wave[x0:x1] = wl
    return wave


def bin_depths(rate_path: Path, x1d_path: Path, t_in, t_out_gap):
    with fits.open(rate_path) as hdul:
        sci = hdul["SCI"].data.astype(float)
        dq = hdul["DQ"].data
        t = hdul["INT_TIMES"].data["int_mid_BJD_TDB"]
    sci[dq != 0] = np.nan
    nints, ny, nx = sci.shape
    h = (t - t[0]) * 24

    prof = np.nanmedian(sci, axis=(0, 2))
    ic = int(np.nanargmax(prof))
    rows = np.arange(ny)
    bg_rows = np.abs(rows - ic) >= BG_MARGIN
    ap_rows = np.abs(rows - ic) <= APER_HALF
    template = np.nanmedian(sci, axis=0)

    col_profile = np.nansum(template[ap_rows], axis=0)
    wave = column_wavelengths(x1d_path, col_profile)
    good_cols = np.isfinite(wave) & (np.nan_to_num(col_profile) > 0)
    wmin, wmax = np.nanmin(wave[good_cols]), np.nanmax(wave[good_cols])
    edges = np.linspace(wmin, wmax, NBINS + 1)

    # per-integration, per-column aperture flux (background-subtracted)
    colflux = np.empty((nints, nx))
    for i in range(nints):
        img = sci[i].copy()
        bad = ~np.isfinite(img)
        img[bad] = template[bad]
        bg = np.nanmedian(img[bg_rows], axis=0)
        img = img - bg[None, :]
        colflux[i] = np.nansum(img[ap_rows], axis=0)

    inn = (h > t_in[0]) & (h < t_in[1])
    out = (h < t_out_gap[0]) | (h > t_out_gap[1])

    res = []
    for k in range(NBINS):
        sel = good_cols & (wave >= edges[k]) & (wave < edges[k + 1])
        if sel.sum() < 5:
            continue
        lc = colflux[:, sel].sum(axis=1)
        # sigma-clip
        med, sd = np.nanmedian(lc), np.nanstd(lc)
        lc[np.abs(lc - med) > 5 * sd] = np.nan
        # linear detrend on out-of-transit baseline
        okout = out & np.isfinite(lc)
        coef = np.polyfit(h[okout], lc[okout], 1)
        lc = lc / np.polyval(coef, h)
        d = 1 - np.nanmedian(lc[inn]) / np.nanmedian(lc[out])
        sd_out = np.nanstd(lc[out])
        err = sd_out * np.sqrt(1 / np.sum(inn & np.isfinite(lc)) + 1 / np.sum(okout))
        res.append({"wave": 0.5 * (edges[k] + edges[k + 1]),
                    "dwave": edges[k + 1] - edges[k],
                    "depth_pct": d * 100, "err_pct": err * 100})
    return pd.DataFrame(res)


def main() -> None:
    visit = sys.argv[1] if len(sys.argv) > 1 else "1"
    cfg = VISITS[visit]
    if cfg["t_in"] is None:
        raise SystemExit(f"visit {visit}: transit window not set — run "
                         f"jwst_whitelight.py {visit} first and fill VISITS")
    raw = ROOT / "data" / "jwst_raw" / f"gj1132b_visit{visit}"
    # science TSO rateints (skip small target-acquisition images), pair with x1dints
    rates = [f for f in sorted(raw.glob("jw*_nrs?_rateints.fits"))
             if f.stat().st_size > 100e6]
    pairs = [(r, r.with_name(r.name.replace("_rateints", "_x1dints"))) for r in rates]
    print(f"visit {visit}: {len(pairs)} science rateints files")
    ours = pd.concat([bin_depths(r, x, cfg["t_in"], cfg["t_out_gap"])
                      for r, x in pairs], ignore_index=True)
    ours = ours.sort_values("wave").reset_index(drop=True)
    # quality cut: drop low-flux red-end bins where the box extraction breaks down
    ncut = (ours.err_pct >= 0.02).sum()
    if ncut:
        print(f"quality cut: dropping {ncut} bins with err >= 200 ppm")
    ours = ours[ours.err_pct < 0.02].reset_index(drop=True)
    ours.to_csv(REPORTS / f"gj1132b_visit{visit}_own_spectrum.csv", index=False)

    # published May et al. 2023 Eureka R100 counterpart from the archive
    idx = pd.read_csv(ROOT / "data" / "raw" / "spectra_index.csv")
    row = idx[(idx.pl_name == "GJ 1132 b") & (idx.authors.str.contains("May"))
              & (idx.note.str.contains(cfg["pub_note"], na=False))].iloc[0]
    pub = usable_points(load_spectrum(ROOT / "data" / "spectra" / local_name(row.spec_path), "Transmission"))

    # comparison statistic (our spectrum in the same normalized frame)
    mine = pd.DataFrame({
        "wave": ours.wave, "dwave": ours.dwave, "value": ours.depth_pct,
        "err_hi": ours.err_pct, "err_lo": ours.err_pct, "lim": 0, "obs_date": np.nan,
    })
    cmp_res = compare_pair(mine, pub)

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.errorbar(pub["wave"].to_numpy(), pub["value"].to_numpy(),
                yerr=[pub["err_lo"].to_numpy(), pub["err_hi"].to_numpy()],
                fmt="s", ms=4, alpha=0.8, color="C1",
                label=f"May et al. 2023 (Eureka!, Visit {visit}) — published")
    ax.errorbar(mine["wave"].to_numpy(), mine["value"].to_numpy(),
                yerr=mine["err_hi"].to_numpy(),
                fmt="o", ms=5, color="C0", label="this work — own extraction from rateints")
    ax.set_xlabel("wavelength [microns]")
    ax.set_ylabel("transit depth [%]")
    title = f"GJ 1132 b visit {visit} — own transmission spectrum vs published"
    if cmp_res:
        title += f"  (chi2_red={cmp_res.chi2/cmp_res.dof:.2f}, p={cmp_res.p_value:.3f}, offset={cmp_res.offset:+.4f}%)"
    ax.set_title(title, fontsize=11)
    ax.legend()
    fig.tight_layout()
    fig.savefig(REPORTS / f"gj1132b_visit{visit}_spectrum_comparison.png", dpi=150)

    print(ours.to_string(index=False))
    if cmp_res:
        print(f"\nvs published: n={cmp_res.n}, offset={cmp_res.offset:+.4f}%, "
              f"chi2_red={cmp_res.chi2/cmp_res.dof:.2f}, p={cmp_res.p_value:.4f}")
    print(f"-> reports/gj1132b_visit{visit}_spectrum_comparison.png")


if __name__ == "__main__":
    main()
