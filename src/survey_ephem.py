"""Fetch transit ephemerides (mid-transit epoch, period, duration) for the
mini-survey hosts from the NASA Exoplanet Archive `pscomppars` table.
Used by survey_analyze.py to assign each JWST visit to the right planet
and to set the transit window, and by limb_asym_* for the orbital and
stellar parameters of the transit-shape model (a/Rs, inclination, e, omega,
limb-darkening inputs).
"""
import io
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed" / "survey_ephemerides.csv"
OUT.parent.mkdir(parents=True, exist_ok=True)   # fresh clones lack data/processed

# MAST target_name -> archive hostname
HOSTS = {
    "K2-18": "K2-18",
    "L-98-59": "L 98-59",
    "LHS-1140": "LHS 1140",
    "TOI-1685": "TOI-1685",
    "TOI-776": "TOI-776",
    "GJ-1132": "GJ 1132",
}

hostlist = ",".join(f"'{h}'" for h in HOSTS.values())
query = (
    "select pl_name,hostname,pl_tranmid,pl_orbper,pl_trandur,"
    "pl_trandep,pl_rade,pl_orbeccen,pl_orblper,pl_ratdor,pl_orbincl,"
    "pl_imppar,pl_ratror,st_rad,st_mass,st_teff,st_logg,st_met,"
    "tran_flag from pscomppars "
    f"where hostname in ({hostlist})"
)
resp = requests.get("https://exoplanetarchive.ipac.caltech.edu/TAP/sync",
                    params={"query": query, "format": "csv"}, timeout=120)
resp.raise_for_status()
eph = pd.read_csv(io.StringIO(resp.text))
eph = eph[eph.tran_flag == 1].drop(columns="tran_flag")
eph["mast_target"] = eph.hostname.map({v: k for k, v in HOSTS.items()})

# --- fill gaps needed by the transit-shape model -------------------------
# a/Rs from Kepler's third law where the archive has none
G = 6.674e-11
M_SUN, R_SUN = 1.98892e30, 6.957e8
need = eph.pl_ratdor.isna() & eph.st_mass.notna() & eph.st_rad.notna()
a_m = (G * eph.st_mass * M_SUN
       * (eph.pl_orbper * 86400.0) ** 2 / (4 * np.pi ** 2)) ** (1 / 3)
eph.loc[need, "pl_ratdor"] = a_m[need] / (eph.st_rad[need] * R_SUN)
eph["ratdor_computed"] = need.astype(int)

# inclination from impact parameter (circular approx: b = a/Rs * cos i)
need = eph.pl_orbincl.isna() & eph.pl_imppar.notna() & eph.pl_ratdor.notna()
eph.loc[need, "pl_orbincl"] = np.degrees(
    np.arccos(eph.pl_imppar[need] / eph.pl_ratdor[need]))

# e/omega: default to circular where unknown, flagged
eph["eccen_assumed"] = eph.pl_orbeccen.isna().astype(int)
eph["pl_orbeccen"] = eph.pl_orbeccen.fillna(0.0)
eph["pl_orblper"] = eph.pl_orblper.fillna(90.0)

eph.to_csv(OUT, index=False)
print(eph.to_string(index=False))
print(f"-> {OUT}")
