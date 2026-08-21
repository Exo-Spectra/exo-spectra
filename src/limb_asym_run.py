"""Limb-asymmetry study — per-visit driver.

For every (visit, detector, planet) with a *complete* transit (both
contacts + margin inside the visit): extract the white light curve with
the survey machinery, fit a limb-darkened symmetric transit + ramp
simultaneously (BIC ramp selection), then test ingress-vs-egress radius
asymmetry (Method A) and folded-residual shape asymmetry (Method B).
Errors: residual-permutation bootstrap x beta red-noise inflation.

Outputs: reports/limb_asymmetry/limb_asym_catalog.csv + a 4-panel PNG
per visit-detector.

Usage:
  python limb_asym_run.py [targets...] [--visit jwXXXXXXXXXXX]
                          [--n-boot N] [--mcmc] [--no-plots]

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
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).parent))
from survey_analyze import (BASE, ROOT, GAP_HALF, segment_colflux,       # noqa: E402
                            white_light, predicted_transits, refine_center)
import limb_asym_model as lam                                            # noqa: E402

OUT = ROOT / "reports" / "limb_asymmetry"
GJ1132_DIRS = [ROOT / "data" / "jwst_raw" / "gj1132b_visit1",
               ROOT / "data" / "jwst_raw" / "gj1132b_visit2"]
MARGIN_MIN_D = 30.0 / (24.0 * 60.0)   # full-transit gate margin floor: 30 min
Q_INIT = (0.4, 0.35)                  # Kipping q1/q2 start (M-dwarf-ish)


def extract_white(seg_files):
    """Concatenate segments -> (t, white-light rel flux, segment boundaries).
    Local variant of survey_analyze.extract_visit: skips wavelengths (not
    needed for white light) and records segment boundary indices, which
    extract_visit discards."""
    ts, cfs, ap, bg = [], [], None, None
    bounds = []
    n = 0
    for f in sorted(seg_files):
        t, cf, _prof, ap, bg = segment_colflux(f, ap, bg)
        ts.append(t)
        cfs.append(cf)
        n += len(t)
        bounds.append(n)
    t = np.concatenate(ts)
    colflux = np.vstack(cfs)
    rel = white_light(colflux)
    del colflux
    return t, rel, bounds[:-1]        # last boundary == end of array


def segment_step_check(t, rel, bounds, sigma, oot, npts: int = 20):
    """Flux step across each segment boundary (medians of `npts` points each
    side) in units of its uncertainty. Returns (flag, worst_z, per-segment
    OOT medians for optional renormalization)."""
    worst = 0.0
    for b in bounds:
        a = rel[max(0, b - npts):b]
        c = rel[b:b + npts]
        a, c = a[np.isfinite(a)], c[np.isfinite(c)]
        if len(a) < 5 or len(c) < 5:
            continue
        z = abs(np.median(c) - np.median(a)) / (sigma * np.sqrt(
            np.pi / 2.0) * np.sqrt(1.0 / len(a) + 1.0 / len(c)))
        worst = max(worst, z)
    return worst > 3.0, worst


def renorm_segments(t, rel, bounds, oot):
    """Divide each segment by its own out-of-transit median (fallback:
    global). Applied only when a boundary step is detected."""
    edges = [0] + list(bounds) + [len(rel)]
    out = rel.copy()
    for a, b in zip(edges[:-1], edges[1:]):
        seg_oot = rel[a:b][oot[a:b] & np.isfinite(rel[a:b])]
        if len(seg_oot) >= 20:
            out[a:b] = rel[a:b] / np.median(seg_oot)
    return out / np.nanmedian(out)


def spot_flag(t, res, sigma, t1, t4, run_len: int = 5, thr: float = 4.0):
    """In-transit residual excursion in the run_len-point moving average,
    measured against the average's OWN robust scatter over the whole visit
    (self-calibrated: red noise is absorbed into the reference level).
    thr=4 keeps the per-visit false-alarm rate at the % level given the
    ~hundreds of effectively independent windows."""
    ok = np.isfinite(res)
    tt, r = t[ok], res[ok]
    if ok.sum() < 10 * run_len:
        return False, 0.0
    kern = np.ones(run_len) / run_len
    sm = np.convolve(r, kern, mode="valid")
    tc = tt[run_len // 2:run_len // 2 + len(sm)]
    ref = 1.4826 * np.median(np.abs(sm - np.median(sm)))
    if ref <= 0:
        return False, 0.0
    inn = (tc >= t1) & (tc <= t4)
    if inn.sum() < run_len:
        return False, 0.0
    worst = float(np.max(np.abs(sm[inn] - np.median(sm))) / ref)
    return worst > thr, worst


def collect_visits(only, only_visit):
    """(target, visit, det) -> [segment files]; survey dirs + GJ 1132 pilot."""
    groups = {}

    def add_dir(target, tdir):
        rates = [f for f in tdir.glob("jw*_nrs?_rateints.fits")
                 if f.stat().st_size > 100e6]
        for f in rates:
            m = re.match(r"(jw\d{11})_\d+_\d+-seg\d+_(nrs\d)_rateints", f.name)
            if m:
                groups.setdefault((target, m.group(1), m.group(2)),
                                  []).append(f)
            else:  # GJ 1132 pilot files have no -segNNN part (single segment)
                m = re.match(r"(jw\d{11})_\d+_\d+_(nrs\d)_rateints", f.name)
                if m:
                    groups.setdefault((target, m.group(1), m.group(2)),
                                      []).append(f)

    for tdir in sorted(BASE.iterdir()):
        if tdir.is_dir() and (not only or tdir.name in only):
            add_dir(tdir.name, tdir)
    if not only or "GJ-1132" in only:
        for d in GJ1132_DIRS:
            if d.exists():
                add_dir("GJ-1132", d)
    if only_visit:
        groups = {k: v for k, v in groups.items() if k[1] == only_visit}
    return groups


def process(target, visit, det, seg_files, planets, n_boot, do_plots,
            do_mcmc):
    rows = []
    t, rel, bounds = extract_white(seg_files)
    cands = predicted_transits(planets, t.min(), t.max())
    if not cands:
        print(f"{visit} {det}: no predicted transit — skipped", flush=True)
        return rows

    for name, c0, dur_d in cands:
        p = planets[planets.pl_name == name].iloc[0]
        c, box_depth = refine_center(t, rel, c0, dur_d)
        others = [(oc, od) for on, oc, od in cands if on != name]

        row = {"target": target, "planet": name, "visit": visit, "det": det,
               "n_int": int(np.isfinite(rel).sum()),
               "e": p.pl_orbeccen, "w_deg": p.pl_orblper,
               "eccen_assumed": int(p.eccen_assumed),
               "partial_flag": 0, "spot_flag": 0, "seg_step_flag": 0}

        # ------------------------------------------------ full-transit gate
        inc0 = lam.impact_to_inc(p.pl_imppar, p.pl_ratdor, p.pl_orbeccen,
                                 p.pl_orblper)
        t1, t2, t3, t4 = lam.contacts(c, p.pl_orbper, p.pl_ratdor, inc0,
                                      p.pl_ratror, p.pl_orbeccen, p.pl_orblper)
        if not np.isfinite([t1, t4]).all():
            print(f"{visit} {det} {name}: grazing/no contacts — skipped",
                  flush=True)
            continue
        margin = max(0.3 * (t4 - t1), MARGIN_MIN_D)
        if t1 - margin < t.min() or t4 + margin > t.max():
            print(f"{visit} {det} {name}: PARTIAL transit "
                  f"(T1={t1:.4f}, T4={t4:.4f}, visit [{t.min():.4f}, "
                  f"{t.max():.4f}]) — catalogued, no asymmetry fit",
                  flush=True)
            row["partial_flag"] = 1
            rows.append(row)
            continue

        # -------------------------------- mask other planets' transits
        use = np.isfinite(rel)
        for oc, od in others:
            use &= np.abs(t - oc) >= GAP_HALF * od
        tt, ff = t[use], rel[use]

        # --------------------------------------- segment boundary check
        sig0 = lam.robust_sigma(ff)
        oot = (t < t1) | (t > t4)
        step_bad, step_z = segment_step_check(t, rel, bounds, sig0, oot)
        row["seg_step_flag"], row["seg_step_z"] = int(step_bad), step_z
        if step_bad:
            rel2 = renorm_segments(t, rel, bounds, oot)
            tt, ff = t[use], rel2[use]
            print(f"{visit} {det} {name}: segment step {step_z:.1f} sigma "
                  f"-> per-segment OOT renormalization", flush=True)

        # ------------------------------------------------ symmetric fit
        init = {"t0": c, "k": p.pl_ratror, "aRs": p.pl_ratdor,
                "b": p.pl_imppar, "q1": Q_INIT[0], "q2": Q_INIT[1]}
        priors = {"aRs": (p.pl_ratdor, 0.1 * p.pl_ratdor),
                  "b": (p.pl_imppar, 0.1)}
        try:
            (sym, ts, fs), bics = lam.fit_with_ramp_selection(
                tt, ff, init, priors, dur_d, p.pl_orbper, p.pl_orbeccen,
                p.pl_orblper)
        except Exception as exc:
            print(f"{visit} {det} {name}: symmetric fit FAILED: {exc}",
                  flush=True)
            continue
        # one pass of residual clipping + refit with the selected ramp
        clip = np.abs(sym.residuals) < 5.0 * sym.sigma
        if (~clip).any():
            sym, ts, fs = lam.fit_transit(
                ts[clip], fs[clip], {n: sym[n] for n in sym.names},
                sym.cfg, priors, dur_d, n_restarts=1)[0], ts[clip], fs[clip]

        inc = lam.impact_to_inc(sym["b"], sym["aRs"], p.pl_orbeccen,
                                p.pl_orblper)
        t1, t2, t3, t4 = lam.contacts(sym["t0"], p.pl_orbper, sym["aRs"],
                                      inc, sym["k"], p.pl_orbeccen,
                                      p.pl_orblper)
        row.update({
            "t0_bmjd": sym["t0"], "depth_ppm": sym["k"] ** 2 * 1e6,
            "k": sym["k"], "aRs": sym["aRs"], "b": sym["b"],
            "q1": sym["q1"], "q2": sym["q2"], "ramp_model": sym.cfg.ramp,
            "t14_h": (t4 - t1) * 24.0 if np.isfinite(t4 - t1) else np.nan,
            "sigma_ppm": sym.sigma * 1e6,
            "chi2_red_sym": sym.chi2 / sym.dof,
            "bic_linear": bics.get("linear", np.nan),
            "bic_quad": bics.get("quad", np.nan),
            "bic_quad_exp": bics.get("quad_exp", np.nan),
            "box_depth_ppm": box_depth * 1e6,
        })

        sflag, sworst = spot_flag(ts, sym.residuals, sym.sigma, t1, t4)
        row["spot_flag"], row["spot_z"] = int(sflag), sworst

        # ----------------------------------------------- asymmetric fit
        asym, ta, fa = lam.fit_asymmetric(ts, fs, sym, priors, dur_d)
        dd = lam.delta_depth_ppm(asym)
        # NOTE: the cyclic-shift bootstrap already carries the red-noise
        # structure; beta inflation on top is conservative and may double
        # count — the injection-recovery stage calibrates which error is
        # correct (see limb_asym_inject.py). Both are catalogued.
        dd_err_boot, _samples = lam.bootstrap_delta_depth(
            ta, fa, asym, priors, dur_d, n_boot=n_boot)
        beta = lam.beta_factor(ta, asym.residuals)
        dd_err = dd_err_boot * beta
        dd_sig = dd / dd_err if dd_err > 0 else np.nan
        row.update({
            "dd_ppm": dd, "dd_err_boot_ppm": dd_err_boot,
            "beta": beta, "dd_err_ppm": dd_err, "dd_sigma": dd_sig,
            "p_value": 2.0 * norm.sf(abs(dd_sig)) if np.isfinite(dd_sig)
            else np.nan,
            "dbic_asym": asym.bic - sym.bic,
            "k_in": asym["k_in"], "k_eg": asym["k_eg"],
        })

        # --------------------------------------------- method B + LD sens
        mb = lam.fold_residual_test(ta, sym.residuals if len(sym.residuals)
                                    == len(ta) else asym.residuals,
                                    sym.sigma, sym["t0"], t1, t4)
        row.update(mb)

        frozen = {"q1": (sym["q1"], 1e-6), "q2": (sym["q2"], 1e-6)}
        asym_ld, _, _ = lam.fit_asymmetric(ts, fs, sym, {**priors, **frozen},
                                           dur_d)
        row["dd_ld_sens_ppm"] = lam.delta_depth_ppm(asym_ld) - dd

        escalate = do_mcmc == "force" or (
            do_mcmc and np.isfinite(dd_sig) and abs(dd_sig) > 2.0)
        if escalate:
            print(f"{visit} {det} {name}: emcee escalation "
                  f"(|dd|={abs(dd_sig):.1f} sigma)", flush=True)
            dd_mc, dd_mc_err, _ = lam.mcmc_asym(ta, fa, asym, priors, dur_d)
            row["dd_mcmc_ppm"], row["dd_mcmc_err_ppm"] = dd_mc, dd_mc_err

        print(f"{visit} {det} {name}: depth {row['depth_ppm']:.0f} ppm, "
              f"ramp {sym.cfg.ramp}, beta {beta:.2f}, "
              f"dd {dd:+.0f} +- {dd_err:.0f} ppm ({dd_sig:+.1f} sigma), "
              f"mB z={mb['mb_z']:+.1f}"
              + (" [SPOT?]" if sflag else ""), flush=True)
        rows.append(row)

        # persist the clipped light curve + models so limb_asym_inject.py
        # does not have to re-read the raw 164 GB archive
        lc_dir = OUT / "lightcurves"
        lc_dir.mkdir(exist_ok=True)
        x = (ta - sym.cfg.t_ref) / sym.cfg.span
        i_ramp = len(sym.names) - lam.RAMP_NPAR[sym.cfg.ramp]
        ramp_s = lam.ramp_flux(sym.cfg.ramp, sym.theta[i_ramp:], x,
                               ta - sym.cfg.t_min)
        pd.DataFrame({"t_bmjd": ta, "flux": fa, "model_sym": sym.model,
                      "model_asym": asym.model, "ramp": ramp_s,
                      "res_sym": sym.residuals}).to_csv(
            lc_dir / f"{target}__{visit}_{det}__"
                     f"{name.replace(' ', '_')}.csv", index=False)

        if do_plots:
            plot_visit(target, visit, det, name, ta, fa, sym, asym,
                       t1, t2, t3, t4)
    return rows


def plot_visit(target, visit, det, planet, t, f, sym, asym, t1, t2, t3, t4):
    cfg = sym.cfg
    x = (t - cfg.t_ref) / cfg.span
    xe = t - cfg.t_min
    i_ramp = len(sym.names) - lam.RAMP_NPAR[cfg.ramp]
    ramp = lam.ramp_flux(cfg.ramp, sym.theta[i_ramp:], x, xe)
    h = (t - sym["t0"]) * 24.0

    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=False)
    ax = axes[0]
    ax.plot(h, f, ".", ms=2, alpha=0.4, color="0.5")
    ax.plot(h, sym.model, "-", color="C1", lw=1.2, label=f"ramp={cfg.ramp}")
    ax.set_ylabel("raw rel. flux")
    ax.legend(fontsize=8)
    ax.set_title(f"{planet} — {visit} {det}")

    ax = axes[1]
    ax.plot(h, f / ramp, ".", ms=2, alpha=0.4, color="0.5")
    ax.plot(h, sym.model / ramp, "-", color="C1", lw=1.2, label="symmetric")
    ax.plot(h, asym.model / ramp, "--", color="C3", lw=1.0, label="asymmetric")
    for tc in (t1, t2, t3, t4):
        if np.isfinite(tc):
            ax.axvline((tc - sym["t0"]) * 24.0, color="0.8", lw=0.6)
    ax.set_ylabel("detrended flux")
    ax.legend(fontsize=8)

    ax = axes[2]
    ax.plot(h, sym.residuals * 1e6, ".", ms=2, alpha=0.4, color="0.5")
    ax.axhline(0, color="C1", lw=0.8)
    ax.set_ylabel("residuals [ppm]")
    ax.set_xlabel("hours from mid-transit")

    ax = axes[3]
    ing = (t >= t1) & (t < sym["t0"])
    egr = (t > sym["t0"]) & (t <= t4)
    nb = 40
    for sel, lab, col in ((ing, "ingress side", "C0"),
                          (egr, "egress side (mirrored)", "C3")):
        dt = np.abs(t[sel] - sym["t0"]) * 24.0
        ff = (f / ramp)[sel]
        bins = np.linspace(0, max((t4 - sym["t0"]), (sym["t0"] - t1)) * 24, nb)
        idx = np.digitize(dt, bins)
        bx = [np.mean(dt[idx == i]) for i in range(1, nb) if (idx == i).sum() > 2]
        by = [np.mean(ff[idx == i]) for i in range(1, nb) if (idx == i).sum() > 2]
        ax.plot(bx, by, "o-", ms=3, lw=0.8, color=col, label=lab, alpha=0.8)
    ax.set_xlabel("|hours from mid-transit|")
    ax.set_ylabel("binned flux")
    ax.legend(fontsize=8)

    fig.tight_layout()
    fn = OUT / f"{target}__{visit}_{det}__{planet.replace(' ', '_')}.png"
    fig.savefig(fn, dpi=150)
    plt.close(fig)


def main():
    args = [a for a in sys.argv[1:]]
    do_mcmc = "force" if "--force-mcmc" in args else ("--mcmc" in args)
    do_plots = "--no-plots" not in args
    n_boot = 300
    if "--n-boot" in args:
        n_boot = int(args[args.index("--n-boot") + 1])
        del args[args.index("--n-boot"):args.index("--n-boot") + 2]
    only_visit = None
    if "--visit" in args:
        only_visit = args[args.index("--visit") + 1]
        del args[args.index("--visit"):args.index("--visit") + 2]
    only = {a for a in args if not a.startswith("--")}

    OUT.mkdir(parents=True, exist_ok=True)
    eph = pd.read_csv(ROOT / "data" / "processed" / "survey_ephemerides.csv")
    groups = collect_visits(only, only_visit)
    print(f"{len(groups)} (target, visit, detector) groups to process")

    all_rows = []
    for (target, visit, det), seg_files in sorted(groups.items()):
        planets = eph[eph.mast_target == target]
        if not len(planets):
            print(f"{target}: no ephemerides — skipped", flush=True)
            continue
        try:
            all_rows.extend(process(target, visit, det, seg_files, planets,
                                    n_boot, do_plots, do_mcmc))
        except Exception as exc:
            print(f"{visit} {det}: FAILED: {exc}", flush=True)

    if not all_rows:
        print("no results")
        return
    cat = pd.DataFrame(all_rows)
    out_csv = OUT / "limb_asym_catalog.csv"
    if out_csv.exists() and (only or only_visit):
        old = pd.read_csv(out_csv)
        key = ["target", "planet", "visit", "det"]
        merged = pd.concat([old, cat]).drop_duplicates(key, keep="last")
        merged.to_csv(out_csv, index=False)
        print(f"{len(cat)} rows updated in {out_csv} ({len(merged)} total)")
    else:
        cat.to_csv(out_csv, index=False)
        print(f"{len(cat)} rows -> {out_csv}")


if __name__ == "__main__":
    main()
