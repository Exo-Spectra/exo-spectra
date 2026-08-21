"""SIDE QUEST (separate from the variability survey — writes only to
reports/planet_hunt/): hunt for unknown planets in the survey white-light
curves.

Two searches per visit (NRS1 detector only — white light needs no more):
  1. unexplained dips: mask known-planet transit windows, detrend the
     baseline (quadratic), slide boxes of several durations and flag
     depth detections with SNR >= 5 outside known windows
  2. O-C timing: matched-filter transit centers of KNOWN planets vs the
     archive ephemeris prediction — systematic offsets hint at unseen
     perturbers (and stale ephemerides; both worth a table)

Usage: python survey_planet_hunt.py [targets...]   (default: all in survey/)
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from survey_analyze import (  # noqa: E402 — reuse the survey extraction
    BASE, ROOT, IN_HALF, GAP_HALF, extract_visit, white_light,
    predicted_transits, refine_center,
)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OUT = ROOT / "reports" / "planet_hunt"
SNR_MIN = 5.0
BOX_HOURS = (0.5, 1.0, 2.0)


def known_mask(t: np.ndarray, known) -> np.ndarray:
    """True where t is OUTSIDE every known planet's (gap-padded) window."""
    out = np.ones_like(t, bool)
    for _, c, dur_d, _ in known:
        out &= np.abs(t - c) >= GAP_HALF * dur_d
    return out


def dip_scan(t, rel, base_ok):
    """Box scan on the detrended baseline; returns candidate events."""
    h = (t - t[0]) * 24
    ok = base_ok & np.isfinite(rel)
    if ok.sum() < 100:
        return [], rel
    coef = np.polyfit(h[ok], rel[ok], 2)
    flat = rel / np.polyval(coef, h)
    sd = np.nanstd(flat[ok])
    events = []
    for dur_h in BOX_HOURS:
        half = dur_h / 48.0  # half-duration in days
        for c in np.arange(t.min() + half, t.max() - half, half / 2):
            inn = (np.abs(t - c) <= half) & base_ok
            outw = (np.abs(t - c) >= 2.4 * half) & base_ok
            if inn.sum() < 15 or outw.sum() < 60:
                continue
            depth = np.nanmedian(flat[outw]) - np.nanmedian(flat[inn])
            err = sd * np.sqrt(1 / inn.sum() + 1 / outw.sum())
            if depth / err >= SNR_MIN:
                events.append({"center_bmjd": c, "dur_h": dur_h,
                               "depth_ppm": depth * 1e6, "snr": depth / err})
    # merge overlapping detections: keep the strongest per cluster
    events.sort(key=lambda e: -e["snr"])
    kept = []
    for e in events:
        if all(abs(e["center_bmjd"] - k["center_bmjd"]) > e["dur_h"] / 24
               for k in kept):
            kept.append(e)
    return kept, flat


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    eph = pd.read_csv(ROOT / "data" / "processed" / "survey_ephemerides.csv")
    only = set(sys.argv[1:])
    cands, ocs = [], []

    for tdir in sorted(BASE.iterdir()):
        if not tdir.is_dir() or (only and tdir.name not in only):
            continue
        target = tdir.name
        planets = eph[eph.mast_target == target]
        groups: dict = {}
        for f in tdir.glob("jw*_nrs1_rateints.fits"):
            if f.stat().st_size < 100e6:
                continue
            m = re.match(r"(jw\d{11})_", f.name)
            groups.setdefault(m.group(1), []).append(f)
        for visit, files in sorted(groups.items()):
            t, colflux, _ = extract_visit(files)
            rel = white_light(colflux)
            del colflux
            known = []
            for name, c0, dur_d in predicted_transits(planets, t.min(), t.max()):
                c, depth = refine_center(t, rel, c0, dur_d)
                known.append((name, c, dur_d, c0))
                ocs.append({"target": target, "visit": visit, "planet": name,
                            "pred_bmjd": c0, "obs_bmjd": c,
                            "o_minus_c_min": (c - c0) * 24 * 60,
                            "white_depth_ppm": depth * 1e6})
            base_ok = known_mask(t, known)
            events, flat = dip_scan(t, rel, base_ok)
            for e in events:
                e.update({"target": target, "visit": visit})
                cands.append(e)
            print(f"{target} {visit}: {len(known)} known transits, "
                  f"{len(events)} unexplained dip(s)", flush=True)

            # diagnostic plot: binned curve + known windows + candidates
            h = (t - t[0]) * 24
            nb = 200
            edges = np.linspace(h.min(), h.max(), nb + 1)
            ib = np.digitize(h, edges) - 1
            hb = np.array([h[ib == k].mean() if (ib == k).any() else np.nan
                           for k in range(nb)])
            fb = np.array([np.nanmedian(flat[ib == k]) if (ib == k).any()
                           else np.nan for k in range(nb)])
            fig, ax = plt.subplots(figsize=(11, 4))
            ax.plot(h, flat, ".", ms=1, alpha=0.2, color="gray")
            ax.plot(hb, fb, "o-", ms=3, color="C0")
            for name, c, dur_d, _ in known:
                ax.axvspan((c - IN_HALF * dur_d - t[0]) * 24,
                           (c + IN_HALF * dur_d - t[0]) * 24,
                           alpha=0.15, color="C2", label=name)
            for e in events:
                ax.axvline((e["center_bmjd"] - t[0]) * 24, color="C3", ls="--")
            ax.set_xlabel("hours since start")
            ax.set_ylabel("relative flux (detrended)")
            ax.set_title(f"{target} {visit}")
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(OUT / f"{target}__{visit}_lightcurve.png", dpi=130)
            plt.close(fig)

    pd.DataFrame(ocs).to_csv(OUT / "known_transit_timing_oc.csv", index=False)
    cdf = pd.DataFrame(cands)
    cdf.to_csv(OUT / "dip_candidates.csv", index=False)
    print(f"\n{len(ocs)} known-transit timings, {len(cdf)} dip candidates")
    if len(cdf):
        print(cdf.to_string(index=False))
    odf = pd.DataFrame(ocs)
    if len(odf):
        print("\nO-C [min] per planet (mean +/- std):")
        print(odf.groupby("planet").o_minus_c_min.agg(["mean", "std", "count"])
              .to_string())


if __name__ == "__main__":
    main()
