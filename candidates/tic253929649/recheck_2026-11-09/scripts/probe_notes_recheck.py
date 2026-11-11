"""Ponowna kontrola kandydatow z not RNAAS #1-#4 po bledzie z bitem QLP 4096 (9.11.2026).

Pytanie: czy sygnal tranzytu przy efemerydzie Z NOTY istnieje, gdy
  (a) pliki QLP maja wyciete kadencje z bitami 2^30 | 2048 | 4096 (obecna maska potoku),
  (b) z plikow SPOC / TESS-SPOC wytniemy te same CHWILE, ktore QLP oznaczyl jako zle
      (okna rozproszonego swiatla dotycza tych samych pikseli, niezaleznie od redukcji)?
Miary dla kazdej serii (pipeline, sektor), dane zbinowane do 10 min jak w notach:
  depth_lin  glebokosc przy efemerydzie liniowej z noty (okno 0,6 T14),
  z_best     najwiekszy z = glebokosc/blad przy przesunieciu srodka w zakresie +-SPAN,
  frac_flag  udzial surowych kadencji in-transit oznaczonych przez QLP.
Test zerowy: te same miary przy 300 losowych fazach (to samo przesuniecie fazy dla
wszystkich serii) -> empiryczny rozklad sredniej wazonej depth_lin i sumy z_best.
Wyjscie: reports/tess_hunt/rnaas_recheck_series.csv, rnaas_recheck_summary.csv
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np, pandas as pd, lightkurve as lk
sys.path.insert(0, str(Path(__file__).parent))
from tess_hunt_deep import prep
from qlp_quality import DEFAULT_BITS

ROOT = Path(__file__).resolve().parents[1]
DATA = str(ROOT / "data" / "tess_hunt" / "qlp")
C = {321290256: (2.331884, 1545.2201, 1.90 / 24, 1.0),
     194121063: (4.732307, 1714.9131, 2.7 / 24, 1.0),
     253929649: (4.814602, 1329.651, 1.1 / 24, 2.0),
     142825877: (7.521089, 1684.36146, 2.4 / 24, 1.5)}   # P, T0, T14, SPAN [h]
NRAND = 300
rng = np.random.default_rng(20261109)


def stat(t, f, P, T0, DUR, span_h):
    ph = ((t - T0 + 0.5 * P) % P) - 0.5 * P
    def one(dt):
        intr = np.abs(ph - dt) < 0.6 * DUR
        outm = (np.abs(ph - dt) > 1.2 * DUR) & (np.abs(ph - dt) < 6 * DUR)
        if intr.sum() < 5 or outm.sum() < 50:
            return np.nan, np.nan, 0
        sig = 1.4826 * np.median(np.abs(f[outm] - np.median(f[outm])))
        d = np.median(f[outm]) - np.mean(f[intr])
        e = sig * np.sqrt(1 / intr.sum() + 1 / outm.sum())
        return d * 1e6, e * 1e6, int(intr.sum())
    d0, e0, n0 = one(0.0)
    zb, dtb = np.nan, np.nan
    for dt in np.arange(-span_h / 24, span_h / 24 + 1e-9, 5 / 1440):
        d, e, n = one(dt)
        if np.isfinite(d) and e > 0 and (not np.isfinite(zb) or d / e > zb):
            zb, dtb = d / e, dt
    return d0, e0, n0, zb, dtb


def wmean(d, e):
    ok = np.isfinite(d) & np.isfinite(e) & (e > 0)
    if ok.sum() == 0:
        return np.nan, np.nan
    w = 1 / e[ok] ** 2
    return np.sum(w * d[ok]) / w.sum(), 1 / np.sqrt(w.sum())


ser_rows, sum_rows = [], []
for tic, (P, T0, DUR, SPAN) in C.items():
    sr = lk.search_lightcurve(f"TIC {tic}", mission="TESS", author=("SPOC", "TESS-SPOC", "QLP"))
    raw, seen = [], set()
    for row in sr:
        key = (str(row.author[0]), str(row.mission[0]), int(row.exptime.value[0]))
        if key in seen:
            continue
        seen.add(key)
        try:
            raw.append((key, row.download(download_dir=DATA)))
        except Exception as ex:
            print("FAIL", tic, key, ex)
    # chwile oznaczone przez QLP, per sektor
    flagged = {}
    for (au, mis, exp), lc in raw:
        if au == "QLP":
            q = np.asarray(lc.quality.value).astype(np.int64)
            tt = np.asarray(lc.time.value)
            flagged[mis] = (tt[(q & DEFAULT_BITS) != 0], exp)
    for variant in ("mask_qlp_only", "mask_all_pipelines"):
        series = []
        for (au, mis, exp), lc in raw:
            frac = np.nan
            tt = np.asarray(lc.time.value)
            ph = ((tt - T0 + 0.5 * P) % P) - 0.5 * P
            intr_raw = np.abs(ph) < 0.6 * DUR
            if mis in flagged and len(flagged[mis][0]):
                ft, fexp = flagged[mis]
                idx = np.clip(np.searchsorted(ft, tt), 1, len(ft) - 1)
                near = np.minimum(np.abs(tt - ft[idx - 1]), np.abs(tt - ft[idx])) < 0.6 * fexp / 86400
                frac = float(near[intr_raw].mean()) if intr_raw.sum() else np.nan
                if variant == "mask_all_pipelines" and au != "QLP":
                    lc = lc[~near]
            elif mis in flagged:
                frac = 0.0
            try:
                t, f, col = prep(lc, bin_min=10)
            except Exception as ex:
                print("prep FAIL", tic, au, mis, ex); continue
            if len(t) < 300:
                continue
            series.append((au, mis, t, f, frac))
        real = []
        for au, mis, t, f, frac in series:
            d0, e0, n0, zb, dtb = stat(t, f, P, T0, DUR, SPAN)
            real.append((d0, e0, zb))
            ser_rows.append(dict(tic=tic, variant=variant, author=au, sector=int(mis.split()[-1]),
                                 depth_lin_ppm=d0, err_ppm=e0, n_in=n0, z_lin=d0 / e0 if e0 else np.nan,
                                 z_best=zb, dt_best_min=dtb * 1440 if np.isfinite(dtb) else np.nan,
                                 frac_in_transit_flagged_by_qlp=frac))
        real = np.array(real, dtype=float)
        null_w, null_s = [], []
        for _ in range(NRAND):
            sh = rng.uniform(0.1, 0.9) * P
            rr = np.array([stat(t, f, P, T0 + sh, DUR, SPAN)[i] for (_, _, t, f, _) in series
                           for i in (0, 1, 3)], dtype=float).reshape(-1, 3)
            null_w.append(wmean(rr[:, 0], rr[:, 1])[0])
            null_s.append(np.nansum(rr[:, 2]))
        null_w, null_s = np.array(null_w), np.array(null_s)
        for grp, sel in (("ALL", lambda a: True), ("QLP", lambda a: a == "QLP"),
                         ("SPOC+TESS-SPOC", lambda a: a != "QLP")):
            m = np.array([sel(s[0]) for s in series])
            dw, ew = wmean(real[m, 0], real[m, 1])
            r = dict(tic=tic, variant=variant, group=grp, n_series=int(m.sum()),
                     depth_w_ppm=dw, err_w_ppm=ew, z_formal=dw / ew if ew else np.nan)
            if grp == "ALL":
                r.update(null_depth_std_ppm=np.nanstd(null_w),
                         z_empirical_depth=(dw - np.nanmean(null_w)) / np.nanstd(null_w),
                         p_null_depth=float((null_w >= dw).mean()),
                         sum_z_best=float(np.nansum(real[:, 2])),
                         z_empirical_sumz=(np.nansum(real[:, 2]) - null_s.mean()) / null_s.std(),
                         p_null_sumz=float((null_s >= np.nansum(real[:, 2])).mean()))
            sum_rows.append(r)
            print(r, flush=True)
pd.DataFrame(ser_rows).to_csv(ROOT / "reports" / "tess_hunt" / "rnaas_recheck_series.csv", index=False)
pd.DataFrame(sum_rows).to_csv(ROOT / "reports" / "tess_hunt" / "rnaas_recheck_summary.csv", index=False)
