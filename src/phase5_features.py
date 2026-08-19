"""Phase 5 (Tier A + C): model-independent per-spectrum structure statistics
and a point-anomaly catalog over the whole archive.

Tier A — for every Transmission/Eclipse spectrum with >= MIN_POINTS usable
points, test the data against three model-free nulls:
    flat    value = const                      (chi2_flat, dof n-1)
    slope   value = c + m * (wave - <wave>_w)  (chi2_slope, dof n-2)
    smooth  value = weighted polynomial        (chi2_smooth, dof n-(d+1))
plus residual diagnostics: lag-1 autocorrelation (coherent bands vs white
noise) and the largest single-point deviation from a local running median.
BH-FDR over p_flat within each spec_type flags "spectrum has significant
structure" — expected physics for real atmospheres, the baseline for Tier B/C.

Tier C — every point deviating > Z_POINT sigma from its own local median is
catalogued; each anomaly is cross-checked for repeatability (do other spectra
of the same planet deviate the same way at that wavelength?) and aggregated
per instrument into wavelength "hotspots" (recurring bin across planets =
suspected instrument systematic, not astrophysics).

Outputs: data/processed/phase5_features.csv,
         data/processed/phase5_point_anomalies.csv,
         data/processed/phase5_instrument_hotspots.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from compare import _sym_sigma
from spectra_io import load_spectrum, usable_points
from phase3_full_archive import bh_fdr, FDR_LEVEL

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed"

MIN_POINTS = 5      # usable points required for the feature set
Z_POINT = 4.0       # local-median deviation threshold for the point catalog
Z_CONFIRM = 2.0     # same-sign deviation in another spectrum counts as confirming
MEDIAN_WINDOW = 5   # rolling-median window for the local continuum
HOTSPOT_BIN = 0.01  # instrument hotspot bins: 1% in log10(wavelength)
SPEC_TYPES = ("Transmission", "Eclipse")


def local_median_z(wave: np.ndarray, val: np.ndarray, sig: np.ndarray) -> np.ndarray:
    """Deviation of each point from a running median of its neighborhood,
    in units of the point's own quoted sigma. Robust to broad molecular
    bands (the median tracks them), sensitive to single-bin spikes."""
    s = pd.Series(val)
    med = s.rolling(MEDIAN_WINDOW, center=True, min_periods=3).median().to_numpy()
    return (val - med) / sig


def spectrum_features(spec: pd.DataFrame) -> dict | None:
    """Tier A statistics for one spectrum (usable points, symmetrized errors)."""
    sig = _sym_sigma(spec)
    ok = np.isfinite(sig) & (sig > 0)
    wave = spec["wave"].to_numpy(float)[ok]
    val = spec["value"].to_numpy(float)[ok]
    sig = sig[ok]
    n = len(val)
    if n < MIN_POINTS:
        return None

    w = 1.0 / sig**2
    # flat: weighted mean
    mean = np.sum(val * w) / np.sum(w)
    chi2_flat = float(np.sum((val - mean) ** 2 * w))
    # slope: c + m*x with weighted-centered x (orthogonal params, cf. compare.py)
    x = wave - np.average(wave, weights=w)
    if np.ptp(x) > 0:
        m = np.sum((val - mean) * x * w) / np.sum(x**2 * w)
        c2 = np.sum((val - m * x) * w) / np.sum(w)
        chi2_slope = float(np.sum((val - c2 - m * x) ** 2 * w))
        slope_snr = float(abs(m) * np.sqrt(np.sum(x**2 * w)))
    else:
        chi2_slope, slope_snr = chi2_flat, 0.0
    # smooth: weighted polynomial, degree scaled to n and capped by the number
    # of distinct wavelengths (repeat photometry can have all points at one wave);
    # fit in a centered/scaled variable for numerical conditioning
    deg = 3 if n >= 20 else (2 if n >= 10 else 1)
    deg = min(deg, len(np.unique(wave)) - 1)
    if deg >= 1:
        u = (wave - wave.mean()) / wave.std()
        coef = np.polyfit(u, val, deg, w=1.0 / sig)
        resid_smooth = (val - np.polyval(coef, u)) / sig
        chi2_smooth = float(np.sum(resid_smooth**2))
        dof_smooth = n - (deg + 1)
    else:  # single distinct wavelength: "smooth" degenerates to the flat model
        chi2_smooth, dof_smooth = chi2_flat, n - 1

    # residual diagnostics on the slope-model residuals
    z = (val - c2 - m * x) / sig if np.ptp(x) > 0 else (val - mean) / sig
    r1 = float(np.sum(z[:-1] * z[1:]) / np.sum(z**2)) if n >= 3 else np.nan
    z_local = local_median_z(wave, val, sig)

    return {
        "n_used": n,
        "chi2_red_flat": chi2_flat / (n - 1),
        "p_flat": float(stats.chi2.sf(chi2_flat, n - 1)),
        "chi2_red_slope": chi2_slope / (n - 2) if n > 2 else np.nan,
        "slope_snr": slope_snr,
        "chi2_red_smooth": chi2_smooth / dof_smooth if dof_smooth > 0 else np.nan,
        "poly_deg": deg,
        "acf_lag1": r1,
        "acf_lag1_snr": r1 * np.sqrt(n) if np.isfinite(r1) else np.nan,
        "max_abs_z_local": float(np.nanmax(np.abs(z_local))),
        "n_pts_gt3sig": int(np.nansum(np.abs(z_local) > 3)),
        "struct_amp": float(np.sqrt(max(np.sum((val - mean) ** 2 * w) / np.sum(w)
                                        - np.mean(sig**2), 0.0)) )   # excess scatter beyond errors, data units
    }


def point_anomalies(sid: int, meta, spec: pd.DataFrame) -> list[dict]:
    """Tier C: points > Z_POINT sigma from the local median."""
    sig = _sym_sigma(spec)
    ok = np.isfinite(sig) & (sig > 0)
    wave = spec["wave"].to_numpy(float)[ok]
    val = spec["value"].to_numpy(float)[ok]
    sig = sig[ok]
    if len(val) < MIN_POINTS:
        return []
    z = local_median_z(wave, val, sig)
    rows = []
    for i in np.nonzero(np.abs(z) > Z_POINT)[0]:
        rows.append({
            "spec_id": sid, "pl_name": meta.pl_name, "spec_type": meta.spec_type,
            "instrument": meta.instrument, "authors": meta.authors, "bibcode": meta.bibcode,
            "wave": wave[i], "value": val[i], "sigma": sig[i], "z_local": float(z[i]),
        })
    return rows


def repeatability(anoms: pd.DataFrame, specs: dict, summary: pd.DataFrame) -> pd.DataFrame:
    """For each anomaly: do other spectra of the same planet & type covering
    that wavelength deviate the same way (|z|>Z_CONFIRM, same sign)?"""
    n_other, n_conf, n_contra = [], [], []
    for a in anoms.itertuples():
        others = summary[(summary.pl_name == a.pl_name)
                         & (summary.spec_type == a.spec_type)
                         & (summary.index != a.spec_id)
                         & (summary.wave_min <= a.wave) & (summary.wave_max >= a.wave)]
        no = nc = nx = 0
        for sid in others.index:
            spec = specs.get(sid)
            if spec is None or len(spec) < MIN_POINTS:
                continue
            sig = _sym_sigma(spec)
            ok = np.isfinite(sig) & (sig > 0)
            wave = spec["wave"].to_numpy(float)[ok]
            val = spec["value"].to_numpy(float)[ok]
            sg = sig[ok]
            # points of the other spectrum near the anomalous wavelength (2% window)
            near = np.abs(wave - a.wave) <= 0.02 * a.wave
            if not near.any():
                continue
            z = local_median_z(wave, val, sg)[near]
            z = z[np.isfinite(z)]
            if len(z) == 0:
                continue
            no += 1
            zpick = z[np.argmax(np.abs(z))]
            if abs(zpick) > Z_CONFIRM and np.sign(zpick) == np.sign(a.z_local):
                nc += 1
            elif abs(zpick) > Z_CONFIRM:
                nx += 1
        n_other.append(no); n_conf.append(nc); n_contra.append(nx)
    anoms = anoms.copy()
    anoms["n_other_specs"] = n_other
    anoms["n_confirming"] = n_conf
    anoms["n_contradicting"] = n_contra
    return anoms


def hotspots(anoms: pd.DataFrame) -> pd.DataFrame:
    """Instrument x wavelength bins that collect anomalies from >= 2 planets:
    recurring bin across different targets = suspected instrument systematic."""
    a = anoms.copy()
    a["logw_bin"] = np.floor(np.log10(a.wave) / HOTSPOT_BIN).astype(int)
    grp = (a.groupby(["instrument", "logw_bin"])
             .agg(n_anomalies=("spec_id", "size"),
                  n_planets=("pl_name", "nunique"),
                  n_spectra=("spec_id", "nunique"),
                  wave_lo=("wave", "min"), wave_hi=("wave", "max"),
                  mean_z=("z_local", "mean"))
             .reset_index())
    hs = grp[(grp.n_anomalies >= 3) & (grp.n_planets >= 2)]
    return hs.sort_values("n_anomalies", ascending=False).drop(columns="logw_bin")


def main() -> None:
    summary = pd.read_csv(OUT / "spectra_summary.csv").set_index("spec_id")
    todo = summary[summary.spec_type.isin(SPEC_TYPES) & (summary.n_usable >= MIN_POINTS)]
    print(f"spectra with >= {MIN_POINTS} usable points: {len(todo)} "
          f"(Transmission {int((todo.spec_type == 'Transmission').sum())}, "
          f"Eclipse {int((todo.spec_type == 'Eclipse').sum())})")

    specs, rows, anom_rows, failed = {}, [], [], []
    for sid, r in todo.iterrows():
        try:
            spec = usable_points(load_spectrum(ROOT / "data" / "spectra" / r.file, r.spec_type))
        except Exception as e:  # noqa: BLE001 — log and continue, phase-1 parity
            failed.append((sid, str(e)))
            continue
        specs[sid] = spec
        feat = spectrum_features(spec)
        if feat is None:
            continue
        rows.append({"spec_id": sid, "pl_name": r.pl_name, "spec_type": r.spec_type,
                     "instrument": r.instrument, "authors": r.authors,
                     "bibcode": r.bibcode, **feat})
        anom_rows.extend(point_anomalies(sid, r, spec))
    if failed:
        print(f"WARNING: {len(failed)} spectra failed to parse: {failed[:5]}")

    feats = pd.DataFrame(rows)
    # BH-FDR within each spec_type: "has significant structure beyond flat"
    feats["structured"] = False
    for st in SPEC_TYPES:
        m = feats.spec_type == st
        if m.any():
            feats.loc[m, "structured"] = bh_fdr(feats.loc[m, "p_flat"].to_numpy(), FDR_LEVEL)
    feats = feats.sort_values("p_flat").reset_index(drop=True)
    feats.to_csv(OUT / "phase5_features.csv", index=False)

    anoms = pd.DataFrame(anom_rows)
    if len(anoms):
        anoms = repeatability(anoms, specs, summary)
        anoms = anoms.sort_values("z_local", key=np.abs, ascending=False).reset_index(drop=True)
    anoms.to_csv(OUT / "phase5_point_anomalies.csv", index=False)
    hs = hotspots(anoms) if len(anoms) else pd.DataFrame()
    hs.to_csv(OUT / "phase5_instrument_hotspots.csv", index=False)

    # ---- console summary + sanity checks -------------------------------
    for st in SPEC_TYPES:
        f = feats[feats.spec_type == st]
        print(f"\n{st}: {len(f)} spectra, structured @FDR{FDR_LEVEL}: "
              f"{int(f.structured.sum())} ({f.structured.mean():.0%}), "
              f"median chi2_red_flat {f.chi2_red_flat.median():.2f}")
    print(f"\npoint anomalies |z|>{Z_POINT}: {len(anoms)} in "
          f"{anoms.spec_id.nunique() if len(anoms) else 0} spectra; "
          f"confirmed by another spectrum: {int((anoms.n_confirming > 0).sum()) if len(anoms) else 0}")
    print(f"instrument hotspots (>=3 anomalies, >=2 planets): {len(hs)}")

    def check(name, planet, col, substr, expect_high):
        f = feats[(feats.pl_name == planet)
                  & feats[col].str.contains(substr, case=False, na=False)
                  & (feats.spec_type == "Transmission")]
        if len(f) == 0:
            print(f"  {name}: no spectra found — SKIP")
            return
        med = f.chi2_red_flat.median()
        frac = f.structured.mean()
        verdict = "OK" if (frac >= 0.5) == expect_high else "CHECK"
        print(f"  {name}: {len(f)} spectra, median chi2_red_flat {med:.2f}, "
              f"structured {frac:.0%} -> {verdict}")

    print("\nsanity checks (Transmission):")
    check("WASP-39 b / JWST NIRSpec (deep features expected)", "WASP-39 b", "instrument", "NIRSpec", True)
    check("HD 189733 b / WFC3 (water band expected)", "HD 189733 b", "instrument", "Wide Field", True)
    check("GJ 1214 b / Kreidberg 2014 (canonically flat)", "GJ 1214 b", "authors", "Kreidberg", False)


if __name__ == "__main__":
    main()
