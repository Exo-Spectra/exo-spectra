"""Blind box search 0.9-15 d on dips and on the inverted light curve (bumps) - the same search on
sign-flipped data is the empirical null for the highest peak. astropy BoxLeastSquares is used as a library."""
import sys, numpy as np, pandas as pd
from astropy.timeseries import BoxLeastSquares
from core import *
PRIORITY = {"SPOC": 0, "TESS-SPOC": 1, "QLP": 2}

def prep_all(raw, cfg):
    out = []
    for s in raw:
        pr = prepare(s, **cfg)
        if pr is None: continue
        pr.update(author=s["author"], sector=s["sector"], name=label(s)); out.append(pr)
    return out

def top_peaks(per, snr, k=8, sep=0.01):
    idx = np.argsort(snr)[::-1]; chosen = []
    for i in idx:
        if all(abs(per[i] - per[j]) / per[j] > sep for j in chosen):
            chosen.append(i)
        if len(chosen) >= k: break
    return chosen

SETS = {
    "one-per-sector SPOC>TESS-SPOC>QLP, strict mask, running median": (dict(mask="strict", flux="auto", detrend="runmed", window=0.75), "ops"),
    "QLP only, strict mask, own detrending of SAP": (dict(mask="strict", flux="auto", detrend="runmed", window=0.75), "qlp"),
    "QLP only, lightkurve default mask, QLP-detrended flux (discovery-like setting)": (dict(mask="lkdef", flux="det", detrend="runmed", window=0.75), "qlp"),
}
rows = []
for tic, e in TARGETS.items():
    P, T0, T14 = e["P"], e["T0"], e["T14"]
    raw = load_raw(tic)
    for sname, (cfg, sel) in SETS.items():
        ser = prep_all(raw, cfg)
        if sel == "ops":
            ser = [min([x for x in ser if x["sector"] == sec], key=lambda x: PRIORITY[x["author"]]) for sec in sorted({x["sector"] for x in ser})]
        else:
            ser = [x for x in ser if x["author"] == "QLP"]
        t = np.concatenate([x["t"] for x in ser]); r = np.concatenate([x["r"] for x in ser])
        dy = np.concatenate([np.full(x["t"].size, x["sigma"]) for x in ser])
        o = np.argsort(t); t, r, dy = t[o], r[o], dy[o]
        base = t[-1] - t[0]
        durs = np.array([0.75, 1.0, 1.5]) * T14
        c = durs.min() / (3.0 * base)
        periods = np.exp(np.arange(np.log(0.9), np.log(15.0), c))
        print(f"\n##### TIC {tic} | {sname}\n  series: {[x['name'] for x in ser]}  n={t.size} baseline={base:.0f} d  n_periods={periods.size} durations[h]={np.round(durs*24,2)}")
        res = {}
        for sign, nm in ((1, "dips"), (-1, "bumps (inverted flux)")):
            pw = BoxLeastSquares(t, sign * r, dy).power(periods, durs, objective="snr", method="fast", oversample=5)
            snr = np.asarray(pw.depth_snr); dep = np.asarray(pw.depth); tt = np.asarray(pw.transit_time); du = np.asarray(pw.duration)
            res[nm] = snr
            print(f"  {nm}: max SNR={snr.max():.2f}; median={np.median(snr):.2f}; 99.9th pct={np.percentile(snr,99.9):.2f}")
            for i in top_peaks(periods, snr):
                print(f"      P={periods[i]:.5f} d  SNR={snr[i]:.2f}  depth={dep[i]*1e6:.0f} ppm  dur={du[i]*24:.2f} h  t0={tt[i]:.4f}")
                rows.append(dict(tic=tic, set=sname, kind=nm, P=periods[i], snr=snr[i], depth_ppm=dep[i]*1e6, dur_h=du[i]*24, t0=tt[i]))
            if sign == 1:
                near = np.abs(periods - P) < 0.003
                j = np.where(near)[0][np.argmax(snr[near])]
                rank = int(np.sum(snr > snr[j])) + 1
                ph = ((tt[j] - T0 + P/2) % periods[j]) - P/2
                print(f"  >> claimed P={P}: best nearby P={periods[j]:.5f} SNR={snr[j]:.2f} depth={dep[j]*1e6:.0f} ppm dur={du[j]*24:.2f} h t0={tt[j]:.4f};"
                      f" fraction of all trial periods with higher SNR = {(rank-1)/periods.size:.5f}")
                claimed = snr[j]
                if tic == 253929649:
                    for pa in (4.878611, 4.852472, 4.814564, 4.869591, 4.789319):
                        nr = np.abs(periods - pa) < 0.003; jj = np.where(nr)[0][np.argmax(snr[nr])]
                        print(f"     listed period {pa}: best nearby P={periods[jj]:.5f} SNR={snr[jj]:.2f} depth={dep[jj]*1e6:.0f} ppm")
        print(f"  >> claimed-period peak SNR {claimed:.2f} vs highest dip {res['dips'].max():.2f} vs highest bump {res['bumps (inverted flux)'].max():.2f};"
              f" number of distinct bump peaks above the claimed peak: see list")
pd.DataFrame(rows).to_csv("out_bls_peaks.csv", index=False)
