"""Noty RNAAS #1 (TIC 253929649) i #2 (TIC 142825877): czy opublikowana statystyka
"timing per seria" (max SNR glebokosci w oknie +-2 h, serie z SNR>=2) jest istotna wzgledem
losowych faz, i jak zalezy od wersji funkcji prep() (9.11.2026).

Warianty prep():
  A_old      clip symetryczny 4 sigma, bez maski QLP  (stan z dnia pisania not, X.2026)
  B_asym     clip tylko w gore (sigma_lower=50), bez maski QLP (stan 7.11)
  C_current  clip tylko w gore + maska QLP 2^30|2048|4096 (stan 9.11)
  D_sym_mask clip symetryczny 4 sigma + maska QLP (kontrola: sama maska)
Statystyka (jak w src/tess_hunt_timing.py::local_dt): dla kazdej serii max po dt z
SNR = (mediana poza - mediana w oknie 0,5 T14) / (sigma / sqrt(n_in)), siatka co 2 min.
Zerowy rozklad: ta sama statystyka dla 400 losowych przesuniec fazy (0,1-0,9 P),
to samo przesuniecie dla wszystkich serii.
Wyjscie: reports/tess_hunt/rnaas_timing_variants.csv
"""
import sys, warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np, pandas as pd, lightkurve as lk
sys.path.insert(0, str(Path(__file__).parent))
from qlp_quality import drop_qlp_flagged, DEFAULT_BITS

ROOT = Path(__file__).resolve().parents[1]
DATA = str(ROOT / "data" / "tess_hunt" / "qlp")
C = {253929649: (4.8146, 1329.63239, 1.1 / 24), 142825877: (7.521089, 1684.36146, 2.0 / 24),
     321290256: (2.331884, 1545.2201, 1.90 / 24), 194121063: (4.732307, 1714.9131, 2.7 / 24)}
VARIANTS = {"A_old": (4, 0), "B_asym": (50, 0), "C_current": (50, DEFAULT_BITS), "D_sym_mask": (4, DEFAULT_BITS)}
NR = 400
rng = np.random.default_rng(1109)


def prep_v(lc, bin_min, sig_lower, bits):
    lc = drop_qlp_flagged(lc, bits=bits) if bits else lc
    cols = lc.colnames
    col = ("pdcsap_flux" if "pdcsap_flux" in cols else "kspsap_flux" if "kspsap_flux" in cols
           else "det_flux" if "det_flux" in cols else "sap_flux")
    use = lc.select_flux(col).remove_nans()
    use = use.remove_outliers(sigma_lower=sig_lower, sigma_upper=4).normalize()
    if bin_min:
        use = use.bin(time_bin_size=bin_min / 60.0 / 24.0).remove_nans()
    use = use.flatten(window_length=min(721, 2 * (len(use) // 6) + 1))
    t = use.time.value; f = np.asarray(use.flux, dtype=float)
    ok = np.isfinite(t) & np.isfinite(f) & (f > 0.5) & (f < 1.5)
    return t[ok], f[ok]


GRID = np.arange(-0.085, 0.0851, 1.0 / 24 / 30)


def best_snr(t, f, sig, P, T0, DUR):
    best = (np.nan, np.nan, np.nan)
    for dt in GRID:
        ph = ((t - (T0 + dt) + 0.5 * P) % P) - 0.5 * P
        intr = np.abs(ph) < 0.5 * DUR
        outm = (np.abs(ph) > 0.75 * DUR) & (np.abs(ph) < 5 * DUR)
        if intr.sum() < 4 or outm.sum() < 30:
            continue
        d = np.median(f[outm]) - np.median(f[intr]); e = sig / np.sqrt(intr.sum())
        if not np.isfinite(best[0]) or d / e > best[0]:
            best = (d / e, dt, d)
    return best


def stats(series, P, T0, DUR):
    b = np.array([best_snr(t, f, s, P, T0, DUR) for (_, t, f, s) in series], dtype=float)
    snr = b[:, 0]
    det = snr >= 2
    return np.nansum(snr), int(det.sum()), np.nansum(np.where(det, snr, 0) ** 2) ** 0.5, b


rows = []
for tic, (P, T0, DUR) in C.items():
    sr = lk.search_lightcurve(f"TIC {tic}", author=("SPOC", "TESS-SPOC", "QLP"))
    raw = []
    for item in sr:
        try:
            raw.append((str(np.atleast_1d(item.author)[0]), float(np.atleast_1d(item.exptime.value)[0]),
                        item.download(download_dir=DATA)))
        except Exception:
            pass
    for vn, (sl, bits) in VARIANTS.items():
        series = []
        for au, exp, lc in raw:
            try:
                t, f = prep_v(lc, 10 if exp < 500 else None, sl, bits)
            except Exception:
                continue
            if len(t) < 100:
                continue
            sig = max(1.4826 * np.median(np.abs(f - np.median(f))), 1e-6)
            series.append((f"{au}-s{lc.meta.get('SECTOR')}", t, f, sig))
        for grp, sel in (("ALL", lambda n: True), ("QLP", lambda n: n.startswith("QLP")),
                         ("SPOC2min+TESS-SPOC", lambda n: not n.startswith("QLP"))):
            ss = [s for s in series if sel(s[0])]
            if not ss:
                continue
            S, ndet, q, b = stats(ss, P, T0, DUR)
            nul = np.array([stats(ss, P, T0 + rng.uniform(0.1, 0.9) * P, DUR)[:3] for _ in range(NR)], dtype=float)
            r = dict(tic=tic, variant=vn, group=grp, n_series=len(ss), sum_snr=round(S, 2),
                     null_mean=round(nul[:, 0].mean(), 2), null_std=round(nul[:, 0].std(), 2),
                     z_emp=round((S - nul[:, 0].mean()) / nul[:, 0].std(), 2),
                     p_emp=float((nul[:, 0] >= S).mean()), n_det_snr2=ndet,
                     null_median_n_det=float(np.median(nul[:, 1])),
                     p_n_det=float((nul[:, 1] >= ndet).mean()),
                     per_series=" ".join(f"{n}:{x[0]:.1f}@{x[1]*1440:+.0f}" for (n, *_), x in zip(ss, b)))
            rows.append(r); print(r, flush=True)
pd.DataFrame(rows).to_csv(ROOT / "reports" / "tess_hunt" / "rnaas_timing_variants.csv", index=False)
