"""Chromatic limb-asymmetry study (gap #8), stage 1: per-wavelength-bin
light curves for every full-transit (visit, detector) of the white-light
study (study-05).

Reuses the survey extraction (segment_colflux) and bins the per-column
aperture flux into NBINS_CHROM uniform wavelength bins per detector
(coarser than the 22-bin survey spectra: each bin must support a light-curve
fit on its own; band-level grouping happens at the statistics stage).
Per bin: 3x5-sigma clip, median normalization, and the same per-segment
out-of-transit renormalization as the white-light study, with the OOT mask
built from the study-05 catalog (t0_bmjd, t14_h, GAP_HALF padding) over ALL
planets transiting in the visit.

Outputs (reports/limb_chrom/):
  lightcurves/<target>__<visit>_<det>__bins.csv  -- t_bmjd + bin00..NN flux
  bins_index.csv  -- (target, visit, det, bin, wave_lo/hi um, n_cols, n_ok)

Run after: study-05 (needs reports/limb_asymmetry/limb_asym_catalog.csv and
the JWST rateints of the survey + GJ 1132 pilot on disk).
Usage: python limb_chrom_extract.py [target ...] [--visit jwNNNNNNNNNNN]
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from survey_analyze import (ROOT, GAP_HALF, segment_colflux)   # noqa: E402
from jwst_spectrum import column_wavelengths                   # noqa: E402
from limb_asym_run import collect_visits, renorm_segments      # noqa: E402

import numpy as np    # noqa: E402
import pandas as pd   # noqa: E402

OUT = ROOT / "reports" / "limb_chrom"
CAT = ROOT / "reports" / "limb_asymmetry" / "limb_asym_catalog.csv"
NBINS_CHROM = 8


def extract_colflux(seg_files):
    """Concatenate segments -> (t, colflux[nints, nx], wave, seg bounds)."""
    ts, cfs, bounds, ap, bg, profile, n = [], [], [], None, None, None, 0
    for f in sorted(seg_files):
        t, cf, prof, ap, bg = segment_colflux(f, ap, bg)
        ts.append(t)
        cfs.append(cf)
        n += len(t)
        bounds.append(n)
        if profile is None:
            profile = prof
    x1d = sorted(seg_files)[0]
    x1d = x1d.with_name(x1d.name.replace("_rateints", "_x1dints"))
    wave = column_wavelengths(x1d, profile)
    return (np.concatenate(ts), np.vstack(cfs), np.asarray(wave),
            bounds[:-1])


def bin_curve(colflux, cols):
    """Normalized, clipped light curve of one wavelength bin."""
    flux = np.nansum(colflux[:, cols], axis=1)
    med, sd = np.nanmedian(flux), np.nanstd(flux)
    ok = np.isfinite(flux)
    for _ in range(3):
        ok = np.abs(flux - med) < 5 * sd
        med, sd = np.nanmedian(flux[ok]), np.nanstd(flux[ok])
    flux = flux.astype(float)
    flux[~ok] = np.nan
    return flux / np.nanmedian(flux)


def main() -> None:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    only_visit = None
    for i, a in enumerate(sys.argv[1:]):
        if a == "--visit":
            only_visit = sys.argv[1:][i + 1]
    only = set(a for a in args if a != only_visit)

    cat = pd.read_csv(CAT)
    cat = cat[cat.partial_flag == 0]
    (OUT / "lightcurves").mkdir(parents=True, exist_ok=True)

    groups = collect_visits(only, only_visit)
    index_rows = []
    for (target, visit, det), segs in sorted(groups.items()):
        rows = cat[(cat.visit == visit) & (cat.det == det)]
        if not len(rows):
            continue                       # partial or not in the catalog
        t, colflux, wave, bounds = extract_colflux(segs)

        # OOT mask over every catalogued transit in this (visit, det)
        oot = np.ones_like(t, bool)
        for _, r in rows.iterrows():
            oot &= np.abs(t - r.t0_bmjd) > GAP_HALF * (r.t14_h / 24.0)

        good = np.isfinite(wave)
        lo, hi = np.nanmin(wave[good]), np.nanmax(wave[good])
        edges = np.linspace(lo, hi, NBINS_CHROM + 1)
        data = {"t_bmjd": t}
        for k in range(NBINS_CHROM):
            cols = np.where(good & (wave >= edges[k])
                            & (wave < edges[k + 1] + (k == NBINS_CHROM - 1)
                               * 1e-9))[0]
            rel = bin_curve(colflux, cols)
            rel = renorm_segments(t, rel, bounds, oot)
            data[f"bin{k:02d}"] = rel
            index_rows.append({
                "target": target, "visit": visit, "det": det, "bin": k,
                "wave_lo_um": edges[k], "wave_hi_um": edges[k + 1],
                "n_cols": len(cols),
                "n_ok": int(np.isfinite(rel).sum()),
                "oot_scatter_ppm": float(np.nanstd(rel[oot]) * 1e6),
            })
        stem = f"{target}__{visit}_{det}__bins.csv"
        pd.DataFrame(data).to_csv(OUT / "lightcurves" / stem, index=False)
        print(f"{target} {visit} {det}: {NBINS_CHROM} bins "
              f"({lo:.2f}-{hi:.2f} um), {len(t)} ints -> {stem}", flush=True)
        del colflux

    idx = pd.DataFrame(index_rows)
    out_idx = OUT / "bins_index.csv"
    if out_idx.exists() and (only or only_visit):   # partial run: merge
        old = pd.read_csv(out_idx)
        key = ["target", "visit", "det", "bin"]
        idx = pd.concat([old[~old.set_index(key).index
                             .isin(idx.set_index(key).index)], idx])
    idx.sort_values(["target", "visit", "det", "bin"]).to_csv(
        out_idx, index=False)
    print(f"index: {len(idx)} rows -> {out_idx}", flush=True)


if __name__ == "__main__":
    main()
