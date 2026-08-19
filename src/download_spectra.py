"""Download all spectrum .tbl files listed in the NASA Exoplanet Archive
Atmospheric Spectroscopy index (data/raw/spectra_index.csv).

File URL pattern (reverse-engineered from the firefly UI, tab1.js):
    https://exoplanetarchive.ipac.caltech.edu/workspace/<TMP_token>/atmospheres/tab1/data/<spec_path>
The TMP token is session-scoped, so we fetch a fresh one from the firefly page
before downloading.
"""
import csv
import re
import sys
import time
from pathlib import Path

import requests

BASE = "https://exoplanetarchive.ipac.caltech.edu"
FIREFLY = BASE + "/cgi-bin/atmospheres/nph-firefly?atmospheres"
ROOT = Path(__file__).resolve().parents[1]
INDEX = ROOT / "data" / "raw" / "spectra_index.csv"
OUT = ROOT / "data" / "spectra"
FAILED_LOG = ROOT / "data" / "raw" / "download_failures.log"


def get_workspace_token(session: requests.Session) -> str:
    html = session.get(FIREFLY, timeout=60).text
    m = re.search(r"FF_InitPage \('([^']+)'", html)
    if not m:
        raise RuntimeError("Could not find workspace token in firefly page")
    # e.g. /work/TMP_xxx_123/atmospheres/tab1 -> TMP_xxx_123
    return m.group(1).split("/")[2]


def local_name(spec_path: str) -> str:
    # "56/12/36/06/Kepler_20_c_3.101_3665_1.tbl" -> "56-12-36-06__Kepler_20_c_3.101_3665_1.tbl"
    *hash_parts, fname = spec_path.split("/")
    return "-".join(hash_parts) + "__" + fname


def main() -> None:
    with open(INDEX, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    paths = [r["spec_path"] for r in rows if r.get("spec_path")]
    print(f"{len(paths)} spectra listed in index")

    OUT.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    token = get_workspace_token(session)
    print(f"workspace token: {token}")

    failures = []
    done = skipped = 0
    for i, spec_path in enumerate(paths, 1):
        dest = OUT / local_name(spec_path)
        if dest.exists() and dest.stat().st_size > 0:
            skipped += 1
            continue
        url = f"{BASE}/workspace/{token}/atmospheres/tab1/data/{spec_path}"
        try:
            r = session.get(url, timeout=60)
            if r.status_code == 200 and r.content.strip():
                dest.write_bytes(r.content)
                done += 1
            else:
                # token may have expired -> refresh once and retry
                token = get_workspace_token(session)
                r = session.get(f"{BASE}/workspace/{token}/atmospheres/tab1/data/{spec_path}", timeout=60)
                if r.status_code == 200 and r.content.strip():
                    dest.write_bytes(r.content)
                    done += 1
                else:
                    failures.append((spec_path, r.status_code))
        except requests.RequestException as e:
            failures.append((spec_path, str(e)))
        if i % 100 == 0:
            print(f"{i}/{len(paths)} (downloaded {done}, skipped {skipped}, failed {len(failures)})")
            sys.stdout.flush()
        time.sleep(0.05)  # be polite to the archive

    if failures:
        with open(FAILED_LOG, "w", encoding="utf-8") as f:
            for p, err in failures:
                f.write(f"{p}\t{err}\n")
    print(f"DONE: downloaded {done}, skipped {skipped}, failed {len(failures)}")
    if failures:
        print(f"failures logged to {FAILED_LOG}")


if __name__ == "__main__":
    main()
