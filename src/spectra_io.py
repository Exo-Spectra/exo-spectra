"""Reading NASA Exoplanet Archive atmospheric spectra (.tbl, IPAC format)
into a normalized pandas DataFrame.

Normalized columns:
    wave      central wavelength [microns]
    dwave     bandwidth [microns] (NaN if absent)
    value     measured quantity: transit depth [%] / eclipse depth [%] / F_lambda
    err_hi    upper uncertainty (positive)
    err_lo    lower uncertainty (positive; archive stores it negative)
    lim       limit flag (nonzero = upper/lower limit, not a measurement)
    obs_date  observation date [JD, days] where available (mostly Eclipse)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import ascii as astro_ascii

# per spec_type: (value column, err1, err2, limit flag)
VALUE_COLS = {
    "Transmission": ("PL_TRANDEP", "PL_TRANDEPERR1", "PL_TRANDEPERR2", "PL_TRANDEPLIM"),
    "Eclipse": ("ESPECLIPDEP", "ESPECLIPDEPERR1", "ESPECLIPDEPERR2", "ESPECLIPDEPLIM"),
    "Direct Imaging": ("FLAM", "FLAMERR1", "FLAMERR2", "FLAMLIM"),
}


def local_name(spec_path: str) -> str:
    """Map index spec_path to the flattened local filename in data/spectra/."""
    *hash_parts, fname = spec_path.split("/")
    return "-".join(hash_parts) + "__" + fname


def load_spectrum(path: str | Path, spec_type: str) -> pd.DataFrame:
    """Load one .tbl file into the normalized DataFrame (sorted by wave)."""
    tab = astro_ascii.read(str(path), format="ipac")
    df = tab.to_pandas()
    vcol, e1, e2, lim = VALUE_COLS[spec_type]
    out = pd.DataFrame({
        "wave": pd.to_numeric(df["CENTRALWAVELNG"], errors="coerce"),
        "dwave": pd.to_numeric(df.get("BANDWIDTH"), errors="coerce"),
        "value": pd.to_numeric(df.get(vcol), errors="coerce"),
        "err_hi": pd.to_numeric(df.get(e1), errors="coerce"),
        "err_lo": -pd.to_numeric(df.get(e2), errors="coerce"),
        "lim": pd.to_numeric(df.get(lim), errors="coerce").fillna(0),
    })
    if "OBS_DATE" in df.columns:
        out["obs_date"] = pd.to_numeric(df["OBS_DATE"], errors="coerce")
    else:
        out["obs_date"] = np.nan
    return out.sort_values("wave").reset_index(drop=True)


def usable_points(spec: pd.DataFrame) -> pd.DataFrame:
    """Rows that are real measurements: finite value, at least one finite error, not a limit."""
    m = (
        spec["wave"].notna()
        & spec["value"].notna()
        & (spec["err_hi"].notna() | spec["err_lo"].notna())
        & (spec["lim"] == 0)
    )
    return spec[m]
