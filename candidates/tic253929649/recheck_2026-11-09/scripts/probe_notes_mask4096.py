"""Kontrola kandydatow z not RNAAS #1-#4 pod katem bitu QLP 4096 (9.11.2026).

Dla kazdego kandydata i kazdego sektora QLP liczy glebokosc lokalna tranzytu przy
efemerydzie z noty, w trzech wariantach maski QLP:
  bits0  = bez maski QLP (stan sprzed 8.11, tak powstaly noty),
  old    = 2^30 | 2048 (stan 8.11),
  new    = 2^30 | 2048 | 4096 (stan od 9.11).
Dodatkowo liczy, ile kadencji in-transit i w calym sektorze ma bit 4096 oraz
jak gleboki jest najnizszy pojedynczy punkt in-transit (test "jedna kadencja").
Srodek tranzytu w sektorze: T_obs z timing_per_series.csv, jesli seria QLP tam jest,
w przeciwnym razie efemeryda liniowa z noty.
Wyjscie: reports/tess_hunt/rnaas_mask4096.csv
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np, pandas as pd, lightkurve as lk
sys.path.insert(0, str(Path(__file__).parent))
import tess_hunt_deep as thd

ROOT = Path(__file__).resolve().parents[1]
D = str(ROOT / "data" / "tess_hunt" / "lc")
# TIC: (P, t0, T14 [d]) z not / README kandydatow
C = {321290256: (2.331884, 1545.2201, 1.90 / 24),
     194121063: (4.732307, 1714.9131, 2.7 / 24),
     253929649: (4.814602, 1329.651, 1.1 / 24),
     142825877: (7.521089, 1684.36146, 2.4 / 24)}
MASKS = {"bits0": 0, "old": (1 << 30) | 2048, "new": (1 << 30) | 2048 | 4096}


def tobs_table(tic):
    f = ROOT / "reports" / "tess_hunt" / f"deep_{tic}" / "timing_per_series.csv"
    out = {}
    if f.exists():
        d = pd.read_csv(f)
        tc = "T_obs" if "T_obs" in d.columns else "T_obs_btjd"
        for _, r in d.iterrows():
            if str(r["series"]).startswith("QLP-s"):
                out[int(str(r["series"]).split("s")[-1])] = float(r[tc])
    return out


def depth(t, f, P, tc, dur):
    dt = ((t - tc) / P + 0.5) % 1 * P - 0.5 * P
    inn = np.abs(dt) < 0.4 * dur
    out = (np.abs(dt) > dur) & (np.abs(dt) < 3 * dur)
    if inn.sum() < 3 or out.sum() < 10:
        return np.nan, np.nan, int(inn.sum()), np.nan
    sd = 1.4826 * np.median(np.abs(f[out] - np.median(f[out])))
    d = np.median(f[out]) - np.median(f[inn])
    e = 1.25 * sd * np.sqrt(1 / inn.sum() + 1 / out.sum())
    return d * 1e6, e * 1e6, int(inn.sum()), (np.median(f[out]) - f[inn].min()) / sd


rows = []
for tic, (P, t0, dur) in C.items():
    tobs = tobs_table(tic)
    sr = lk.search_lightcurve(f"TIC {tic}", mission="TESS", author="QLP")
    for row in sr:
        sec = int(str(row.mission[0]).split()[-1])
        try:
            lc = row.download(download_dir=D)
        except Exception as ex:
            print("FAIL", tic, sec, ex); continue
        q = np.asarray(lc.quality.value).astype(np.int64)
        tt = lc.time.value
        tc = tobs.get(sec, t0)
        dtr = ((tt - tc) / P + 0.5) % 1 * P - 0.5 * P
        inn_raw = np.abs(dtr) < 0.5 * dur
        r = dict(tic=tic, sector=sec, tc_source="T_obs" if sec in tobs else "linear",
                 n_cad=len(tt), n_4096=int(((q & 4096) > 0).sum()),
                 n_in_raw=int(inn_raw.sum()), n_in_4096=int(((q[inn_raw] & 4096) > 0).sum()),
                 n_in_bit30=int(((q[inn_raw] & (1 << 30)) > 0).sum()))
        for name, bits in MASKS.items():
            os.environ["QLP_MASK_BITS"] = str(bits)
            t, f, _ = thd.prep(lc)
            d, e, n, deepest = depth(t, f, P, tc, dur) if len(t) else (np.nan,) * 4
            r[f"depth_{name}"] = d; r[f"err_{name}"] = e; r[f"nin_{name}"] = n
            r[f"deepest_pt_sigma_{name}"] = deepest
        rows.append(r)
        print(r, flush=True)
df = pd.DataFrame(rows)
df.to_csv(ROOT / "reports" / "tess_hunt" / "rnaas_mask4096.csv", index=False)
print("\n=== srednie wazone po sektorach ===")
for tic, g in df.groupby("tic"):
    s = f"TIC {tic}:"
    for name in MASKS:
        ok = np.isfinite(g[f"depth_{name}"]) & (g[f"err_{name}"] > 0)
        w = 1 / g.loc[ok, f"err_{name}"] ** 2
        s += f" {name} {np.sum(w * g.loc[ok, f'depth_{name}']) / w.sum():.0f}+-{1 / np.sqrt(w.sum()):.0f} ppm ({ok.sum()} sekt.);"
    s += f" kadencje in-transit z 4096: {g.n_in_4096.sum()}/{g.n_in_raw.sum()}"
    print(s)
