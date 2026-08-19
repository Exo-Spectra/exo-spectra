"""Mini-survey stage 1: uniform re-extraction of every public G395H/M visit
of the top-5 multi-visit hosts with the SAME method as the GJ 1132 b pilot,
then pairwise visit-to-visit comparison per planet.

Per visit and detector:
  1. per-integration, per-column aperture photometry from rateints segments
     (column background, DQ in-fill from median template) — segments
     concatenated in time
  2. wavelengths from the x1dints pipeline log (see jwst_spectrum)
  3. the visit is assigned to a planet via archive ephemerides
     (survey_ephemerides.csv); visits with no predicted transit in-window
     (e.g. eclipse observations) are skipped and logged
  4. transit center refined by a matched box filter around the prediction;
     in-transit = |t-c| < 0.375*dur, baseline = |t-c| > 0.6*dur, minus the
     windows of any other planet transiting in the same visit
  5. NBINS per-detector wavelength bins -> depths + errors (same estimator
     as the pilot), quality cut err < max(200 ppm, 3x median)

Then: all visit pairs per (target, planet) through compare_pair
(offset and offset+slope models), BH-FDR at 1% over the pair set.

Run after: survey_download.py, survey_ephem.py
"""
import os
import re
import sys
from pathlib import Path

# cap BLAS thread pools BEFORE numpy import — each thread preallocates
# buffers and the analysis must coexist with a loaded 32 GB system
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).parent))
from jwst_spectrum import column_wavelengths, APER_HALF, BG_MARGIN, NBINS  # noqa: E402
from compare import compare_pair  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "jwst_raw" / "survey"
OUT = ROOT / "reports" / "survey"
IN_HALF = 0.375   # in-transit half-window, in units of transit duration
GAP_HALF = 0.6    # baseline starts beyond this, in units of duration
ERR_CUT_PPM = 200.0


def segment_colflux(path: Path, ap_rows=None, bg_rows=None):
    """Background-subtracted aperture flux per (integration, column)."""
    with fits.open(path) as hdul:
        sci = hdul["SCI"].data.astype(np.float32)  # float32: halve the footprint
        dq = hdul["DQ"].data
        t = np.asarray(hdul["INT_TIMES"].data["int_mid_BJD_TDB"], float)
    sci[dq != 0] = np.nan
    del dq
    nints, ny, nx = sci.shape
    if ap_rows is None:
        prof = np.nanmedian(sci, axis=(0, 2))
        ic = int(np.nanargmax(prof))
        rows = np.arange(ny)
        bg_rows = np.abs(rows - ic) >= BG_MARGIN
        ap_rows = np.abs(rows - ic) <= APER_HALF
    template = np.nanmedian(sci, axis=0)
    colflux = np.empty((nints, nx))
    for i in range(nints):
        img = sci[i].copy()
        bad = ~np.isfinite(img)
        img[bad] = template[bad]
        bg = np.nanmedian(img[bg_rows], axis=0)
        img = img - bg[None, :]
        colflux[i] = np.nansum(img[ap_rows], axis=0)
    col_profile = np.nansum(template[ap_rows], axis=0)
    return t, colflux, col_profile, ap_rows, bg_rows


def extract_visit(seg_files: list[Path]):
    """Concatenate segments of one (visit, detector) exposure."""
    ts, cfs, ap, bg, profile = [], [], None, None, None
    for f in sorted(seg_files):
        t, cf, prof, ap, bg = segment_colflux(f, ap, bg)
        ts.append(t)
        cfs.append(cf)
        if profile is None:
            profile = prof
    x1d = sorted(seg_files)[0]
    x1d = x1d.with_name(x1d.name.replace("_rateints", "_x1dints"))
    wave = column_wavelengths(x1d, profile)
    return np.concatenate(ts), np.vstack(cfs), wave


def white_light(colflux: np.ndarray) -> np.ndarray:
    flux = np.nansum(colflux, axis=1)
    med, sd = np.nanmedian(flux), np.nanstd(flux)
    for _ in range(3):
        ok = np.abs(flux - med) < 5 * sd
        med, sd = np.nanmedian(flux[ok]), np.nanstd(flux[ok])
    flux[~ok] = np.nan
    return flux / np.nanmedian(flux)


def predicted_transits(planets: pd.DataFrame, t0: float, t1: float):
    """(pl_name, predicted center BJD, duration [d]) inside the visit."""
    cands = []
    for _, p in planets.iterrows():
        dur_d = p.pl_trandur / 24.0
        # archive tranmid is full JD; JWST INT_TIMES BJD_TDB is JD - 2400000.5
        tmid = p.pl_tranmid - 2400000.5
        n = np.round((0.5 * (t0 + t1) - tmid) / p.pl_orbper)
        for k in (n - 1, n, n + 1):
            c = tmid + k * p.pl_orbper
            if t0 + 0.2 * dur_d < c < t1 - 0.2 * dur_d:
                cands.append((p.pl_name, c, dur_d))
    return cands


def refine_center(t, rel, c0, dur_d):
    best, best_d = c0, -np.inf
    for c in np.arange(c0 - 0.05, c0 + 0.05, 0.0005):
        inn = np.abs(t - c) <= IN_HALF * dur_d
        out = np.abs(t - c) >= GAP_HALF * dur_d
        if inn.sum() < 10 or out.sum() < 20:
            continue
        d = np.nanmedian(rel[out]) - np.nanmedian(rel[inn])
        if d > best_d:
            best, best_d = c, d
    return best, best_d


def masks(t, c, dur_d, others):
    inn = np.abs(t - c) <= IN_HALF * dur_d
    out = np.abs(t - c) >= GAP_HALF * dur_d
    for oc, odur in others:
        out &= np.abs(t - oc) >= GAP_HALF * odur
    return inn, out


def bin_depths(t, colflux, wave, inn, out):
    good = np.isfinite(wave) & (np.nan_to_num(np.nanmedian(colflux, axis=0)) > 0)
    if good.sum() < NBINS:
        return pd.DataFrame()
    edges = np.linspace(np.nanmin(wave[good]), np.nanmax(wave[good]), NBINS + 1)
    h = (t - t[0]) * 24
    res = []
    for k in range(NBINS):
        sel = good & (wave >= edges[k]) & (wave < edges[k + 1])
        if sel.sum() < 5:
            continue
        lc = colflux[:, sel].sum(axis=1)
        med, sd = np.nanmedian(lc), np.nanstd(lc)
        lc[np.abs(lc - med) > 5 * sd] = np.nan
        okout = out & np.isfinite(lc)
        if okout.sum() < 20 or (inn & np.isfinite(lc)).sum() < 10:
            continue
        coef = np.polyfit(h[okout], lc[okout], 1)
        lcn = lc / np.polyval(coef, h)
        d = 1 - np.nanmedian(lcn[inn]) / np.nanmedian(lcn[out])
        sd_out = np.nanstd(lcn[out])
        err = sd_out * np.sqrt(1 / np.sum(inn & np.isfinite(lcn)) + 1 / okout.sum())
        res.append({"wave": 0.5 * (edges[k] + edges[k + 1]),
                    "dwave": edges[k + 1] - edges[k],
                    "depth_pct": d * 100, "err_pct": err * 100})
    return pd.DataFrame(res)


def main() -> None:
    (OUT / "spectra").mkdir(parents=True, exist_ok=True)
    eph = pd.read_csv(ROOT / "data" / "processed" / "survey_ephemerides.csv")
    meta_rows = []

    only = set(sys.argv[1:])  # optional target filter (default: all)
    for tdir in sorted(BASE.iterdir()):
        if not tdir.is_dir():
            continue
        target = tdir.name
        if only and target not in only:
            print(f"\n=== {target}: skipped (not in target filter) ===")
            continue
        planets = eph[eph.mast_target == target]
        rates = [f for f in tdir.glob("jw*_nrs?_rateints.fits")
                 if f.stat().st_size > 100e6]
        groups: dict = {}
        for f in rates:
            m = re.match(r"(jw\d{11})_\d+_\d+-seg\d+_(nrs\d)_rateints", f.name)
            if m:
                groups.setdefault((m.group(1), m.group(2)), []).append(f)
        visits = sorted({v for v, _ in groups})
        print(f"\n=== {target}: {len(visits)} science visits ===", flush=True)

        for visit in visits:
            dets = sorted(d for v, d in groups if v == visit)
            parts, info = [], None
            for det in dets:
                t, colflux, wave = extract_visit(groups[(visit, det)])
                rel = white_light(colflux)
                if info is None:
                    cands = predicted_transits(planets, t.min(), t.max())
                    if not cands:
                        print(f"{visit}: no predicted transit in window "
                              f"(eclipse/phase obs?) — SKIPPED", flush=True)
                        break
                    if len(cands) > 1:
                        print(f"{visit}: {len(cands)} transits in window "
                              f"({[c[0] for c in cands]}) — analyzing each "
                              f"with the others' windows masked", flush=True)
                    info = []
                    for name, c0, dur_d in cands:
                        c, depth = refine_center(t, rel, c0, dur_d)
                        info.append((name, c, dur_d, depth))
                        print(f"{visit} {det} {name}: center {c:.4f} "
                              f"(pred {c0:.4f}), dur {dur_d*24:.2f} h, "
                              f"white depth {depth*1e6:.0f} ppm", flush=True)
                for name, c, dur_d, _ in info:
                    others = [(oc, od) for on, oc, od, _ in info if on != name]
                    inn, out = masks(t, c, dur_d, others)
                    df = bin_depths(t, colflux, wave, inn, out)
                    if len(df):
                        df["planet"], df["det"] = name, det
                        parts.append(df)
            if not parts:
                continue
            allb = pd.concat(parts, ignore_index=True)
            for name, g in allb.groupby("planet"):
                g = g.sort_values("wave").reset_index(drop=True)
                cut = max(ERR_CUT_PPM / 1e4, 3 * g.err_pct.median())
                ncut = int((g.err_pct >= cut).sum())
                if ncut:
                    print(f"{visit} {name}: quality cut drops {ncut} bins "
                          f"(err >= {cut*1e4:.0f} ppm)", flush=True)
                g = g[g.err_pct < cut]
                fn = OUT / "spectra" / (f"{target}__{visit}__"
                                        f"{name.replace(' ', '_')}.csv")
                g.to_csv(fn, index=False)
                meta_rows.append({"target": target, "visit": visit,
                                  "planet": name, "n_bins": len(g),
                                  "t0_bjd": t.min(), "file": fn.name})

    meta = pd.DataFrame(meta_rows)
    meta.to_csv(OUT / "survey_visits.csv", index=False)
    print(f"\n{len(meta)} visit-spectra extracted")

    # pairwise comparisons per (target, planet)
    rows = []
    for (target, planet), g in meta.groupby(["target", "planet"]):
        g = g.sort_values("t0_bjd").reset_index(drop=True)
        specs = {}
        for _, r in g.iterrows():
            d = pd.read_csv(OUT / "spectra" / r.file)
            specs[r.visit] = pd.DataFrame({
                "wave": d.wave, "dwave": d.dwave, "value": d.depth_pct,
                "err_hi": d.err_pct, "err_lo": d.err_pct, "lim": 0,
                "obs_date": np.nan})
        vs = list(specs)
        for i in range(len(vs)):
            for j in range(i + 1, len(vs)):
                r = compare_pair(specs[vs[i]], specs[vs[j]])
                if r is None:
                    continue
                has_slope = r.slope is not None
                rows.append({
                    "target": target, "planet": planet,
                    "visit_a": vs[i], "visit_b": vs[j], "n": r.n,
                    "offset_ppm": r.offset * 1e4,
                    "chi2_red": r.chi2 / r.dof, "p_value": r.p_value,
                    "slope_ppm_um": r.slope * 1e4 if has_slope else np.nan,
                    "slope_sigma": (abs(r.slope / r.slope_err)
                                    if has_slope else np.nan),
                    "chi2_red_slope": (r.chi2_slope / r.dof_slope
                                       if has_slope else np.nan),
                    "p_slope": r.p_slope if has_slope else np.nan,
                })
        # overlay plot per planet
        fig, ax = plt.subplots(figsize=(11, 6))
        for k, (v, s) in enumerate(specs.items()):
            ax.errorbar(s.wave, s.value, yerr=s.err_hi, fmt="o", ms=4,
                        alpha=0.8, color=f"C{k}", label=v)
        ax.set_xlabel("wavelength [microns]")
        ax.set_ylabel("transit depth [%]")
        ax.set_title(f"{planet} — {len(specs)} visits, own uniform extraction")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT / f"{target}__{planet.replace(' ', '_')}.png", dpi=150)
        plt.close(fig)

    if not rows:
        print("\nno comparable visit pairs — nothing to test")
        return
    pairs = pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)
    # Benjamini-Hochberg FDR 1% on the offset-model p-values
    m = len(pairs)
    thresh = 0.0
    for i, p in enumerate(np.sort(pairs.p_value.to_numpy()), start=1):
        if p <= 0.01 * i / m:
            thresh = p
    pairs["significant"] = pairs.p_value <= thresh
    pairs.to_csv(OUT / "survey_pairs.csv", index=False)
    print(f"\n{m} visit pairs, {int(pairs.significant.sum())} significant @FDR1%")
    print(pairs.to_string(index=False, max_colwidth=18))


if __name__ == "__main__":
    main()
