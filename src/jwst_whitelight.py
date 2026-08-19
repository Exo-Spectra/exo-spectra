"""Extract a white-light transit curve of GJ 1132 b from JWST NIRSpec
RATEINTS files (per-integration slope images) using our own extraction:

  per integration:
    1. column-wise background: median of rows far from the spectral trace
    2. aperture sum of rows around the trace -> flux

Then normalize per detector and plot flux vs time. If the reduction is sane,
the ~0.26%-deep transit of GJ 1132 b must be visible.

Usage: python jwst_whitelight.py [visit]   (default: 1)
Also reports the transit window found by a matched box filter (same box
geometry as visit 1: in-transit half-width 0.3 h, baseline gap 0.475 h).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

VISIT = sys.argv[1] if len(sys.argv) > 1 else "1"
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "jwst_raw" / f"gj1132b_visit{VISIT}"
REPORTS = ROOT / "reports"
IN_HALF = 0.3     # h, in-transit half-width of the box (visit-1 calibrated)
GAP_HALF = 0.475  # h, |t-c| beyond this counts as baseline
APER_HALF = 6     # rows around trace center
BG_MARGIN = 10    # rows beyond this |offset| used as background


def white_light(path: Path):
    with fits.open(path) as hdul:
        sci = hdul["SCI"].data.astype(float)      # (nints, ny, nx)
        dq = hdul["DQ"].data
        t = hdul["INT_TIMES"].data["int_mid_BJD_TDB"]
    sci[dq != 0] = np.nan

    # trace = row of max median flux (stable for BOTS)
    prof = np.nanmedian(sci, axis=(0, 2))
    ic = int(np.nanargmax(prof))
    ny = sci.shape[1]
    rows = np.arange(ny)
    bg_rows = np.abs(rows - ic) >= BG_MARGIN
    ap_rows = np.abs(rows - ic) <= APER_HALF

    # median image over time -> template used to in-fill DQ-masked pixels,
    # so that flagged pixels don't read as missing flux
    template = np.nanmedian(sci, axis=0)

    flux = np.empty(sci.shape[0])
    for i in range(sci.shape[0]):
        img = sci[i].copy()
        bad = ~np.isfinite(img)
        img[bad] = template[bad]
        bg = np.nanmedian(img[bg_rows], axis=0)   # per-column background
        img = img - bg[None, :]
        flux[i] = np.nansum(img[ap_rows])

    # 5-sigma clip in time (cosmic-ray hit integrations)
    med, sd = np.nanmedian(flux), np.nanstd(flux)
    for _ in range(3):
        ok = np.abs(flux - med) < 5 * sd
        med, sd = np.nanmedian(flux[ok]), np.nanstd(flux[ok])
    flux[~ok] = np.nan
    return t, flux


def find_transit(h: np.ndarray, rel: np.ndarray):
    """Matched box filter: transit center maximizing the in/out depth."""
    best = (None, -np.inf)
    for c in np.arange(GAP_HALF, h.max() - GAP_HALF, 0.01):
        inn = np.abs(h - c) <= IN_HALF
        out = np.abs(h - c) >= GAP_HALF
        if inn.sum() < 20 or out.sum() < 20:
            continue
        d = np.nanmedian(rel[out]) - np.nanmedian(rel[inn])
        if d > best[1]:
            best = (c, d)
    return best


def main() -> None:
    files = sorted(RAW.glob("jw*04102*rateints.fits"))  # science TSO only
    print(f"{len(files)} rateints files")
    fig, axes = plt.subplots(len(files), 1, figsize=(11, 4 * len(files)), sharex=True)
    if len(files) == 1:
        axes = [axes]
    for ax, f in zip(axes, files):
        det = "NRS1" if "nrs1" in f.name else "NRS2"
        t, flux = white_light(f)
        norm = np.nanmedian(flux)
        rel = flux / norm
        # light binning for display
        nb = 64
        edges = np.linspace(t.min(), t.max(), nb + 1)
        ib = np.digitize(t, edges) - 1
        tb = np.array([t[ib == k].mean() for k in range(nb)])
        fb = np.array([np.nanmedian(rel[ib == k]) for k in range(nb)])
        ax.plot((t - t[0]) * 24, rel, ".", ms=2, alpha=0.3, color="gray")
        ax.plot((tb - t[0]) * 24, fb, "o-", ms=4, color="C0")
        ax.set_ylabel(f"{det} relative flux")
        ax.set_title(f"{f.name}  (n_int={len(t)})")
        rms = np.nanstd(rel) * 1e6
        h = (t - t[0]) * 24
        c, depth = find_transit(h, rel)
        print(f"{det}: {len(t)} integrations, per-int scatter {rms:.0f} ppm; "
              f"transit center {c:.2f} h, depth {depth*1e6:.0f} ppm, "
              f"T_IN=({c-IN_HALF:.2f}, {c+IN_HALF:.2f}), "
              f"T_OUT_GAP=({c-GAP_HALF:.2f}, {c+GAP_HALF:.2f})")
        ax.axvspan(c - IN_HALF, c + IN_HALF, alpha=0.1, color="C3")
    axes[-1].set_xlabel("hours since start (BJD_TDB)")
    fig.suptitle(f"GJ 1132 b visit {VISIT} (NIRSpec G395H) — own extraction from rateints")
    fig.tight_layout()
    out = REPORTS / f"gj1132b_visit{VISIT}_whitelight.png"
    fig.savefig(out, dpi=150)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
