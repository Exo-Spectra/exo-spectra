"""Mini-survey stage 1 download: RATEINTS + X1DINTS for all PUBLIC
NIRSpec BOTS G395H/M science visits of the top-5 multi-visit exoplanet
hosts (from mast_survey_scout.py). Idempotent: astroquery skips files
already present with the right size.

Layout: data/jwst_raw/survey/<target>/<files, flat>
"""
import socket
import time
from pathlib import Path

from astroquery.mast import Observations

# dead-TCP guard: without this a dropped connection can hang iter_content
# forever (observed twice: 74 min and 20+ min stalls on a 0-byte file)
socket.setdefaulttimeout(120)

RETRIES = 8


def download_with_retry(sel, dest: Path) -> None:
    """Network drops (sleep, Wi-Fi, VPN) abort download_products mid-file;
    already-complete files are cached, so a retry resumes where it stopped."""
    for attempt in range(1, RETRIES + 1):
        try:
            Observations.download_products(sel, download_dir=str(dest), flat=True)
            return
        except Exception as e:  # noqa: BLE001 — any transport error
            print(f"attempt {attempt}/{RETRIES} failed: {e!r}", flush=True)
            if attempt == RETRIES:
                raise
            time.sleep(30 * attempt)

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "jwst_raw" / "survey"

TARGETS = {
    "K2-18": "F290LP;G395H",
    "L-98-59": "F290LP;G395H",
    "LHS-1140": "F290LP;G395M",
    "TOI-1685": "F290LP;G395H",
    "TOI-776": "F290LP;G395H",
}

for target, filt in TARGETS.items():
    dest = BASE / target
    dest.mkdir(parents=True, exist_ok=True)
    obs = Observations.query_criteria(
        obs_collection="JWST", instrument_name="NIRSPEC/SLIT",
        target_name=target, filters=filt,
    )
    df = obs.to_pandas()
    df = df[df["obs_id"].str.match(r"^jw\d{5}-o\d+") & (df["dataRights"] == "PUBLIC")]
    print(f"\n=== {target}: {len(df)} public visits ===", flush=True)
    if not len(df):
        continue
    prods = Observations.get_product_list(obs[[o in set(df["obs_id"]) for o in obs["obs_id"]]])
    want = [d in ("RATEINTS", "X1DINTS") for d in prods["productSubGroupDescription"]]
    sel = prods[want]
    print(f"{len(sel)} files, {sum(sel['size'])/1e9:.1f} GB -> {dest}", flush=True)
    download_with_retry(sel, dest)

print("\nsurvey download complete")
