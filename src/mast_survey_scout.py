"""Reconnaissance for the mini-survey: which targets did JWST observe MORE
THAN ONCE with NIRSpec BOTS G395H/G395M (the mode our extraction handles)?

Read-only MAST metadata query — no downloads. Output: per-target visit
counts, public/proprietary split, rough download size estimate
(~2.4 GB per G395H visit: 2x rateints + x1dints, from GJ 1132 b).
"""
from collections import defaultdict

import pandas as pd
from astroquery.mast import Observations

obs = Observations.query_criteria(
    obs_collection="JWST",
    instrument_name="NIRSPEC/SLIT",
    filters=["F290LP;G395H", "F290LP;G395M"],
)
df = obs.to_pandas()
# keep top-level per-visit observations only (drop per-exposure duplicates)
df = df[df["obs_id"].str.match(r"^jw\d{5}-o\d+")]
print(f"{len(df)} BOTS G395H/G395M visits in MAST\n")

rows = []
for (target, filt), g in df.groupby(["target_name", "filters"]):
    n_public = int((g["dataRights"] == "PUBLIC").sum())
    rows.append({
        "target": target, "grating": filt.split(";")[1],
        "visits": len(g), "public": n_public,
        "proposals": ",".join(sorted(set(g["proposal_id"].astype(str)))),
    })
r = pd.DataFrame(rows)

multi = (r[r.visits >= 2]
         .sort_values(["visits", "target"], ascending=[False, True])
         .reset_index(drop=True))
print("=== targets with >=2 visits (same grating) ===")
print(multi.to_string(index=False))

pub = multi[multi.public >= 2]
est_gb = (pub.public * 2.4).sum()
print(f"\n{len(pub)} targets with >=2 PUBLIC visits; "
      f"{int(pub.public.sum())} public visits total; "
      f"~{est_gb:.0f} GB to download (rateints+x1dints)")

# transiting-exoplanet hosts only (curated from the raw list — the mode is
# also used for SNe, quasars, brown dwarfs and calibration stars)
EXO_HOSTS = {
    "GJ-1132", "GJ-1214", "GJ-3470", "GJ-4102", "GJ3090", "GJ9827",
    "K2-18", "L-98-59", "TOI-175", "LHS-1140", "LTT-1445A",
    "TOI-455-revised", "TOI-1685", "TOI-776", "TOI-134", "TOI-260",
    "TOI-270", "TOI-402", "TOI-836", "TOI-125", "WOLF-437", "TOI-3884b",
    "TOI-824", "TOI-178", "TOI-849", "TOI-1130", "TOI-1801", "TOI-674",
    "TOI-677", "V-V1298-Tau", "HAT-P-11", "HAT-P-26", "WASP-17",
    "WASP-39", "WASP-47", "WD1856+534",
}
exo = multi[multi.target.isin(EXO_HOSTS) & (multi.public >= 2)]
print("\n=== transiting-exoplanet hosts, >=2 public visits ===")
print(exo.to_string(index=False))
print(f"\n{exo.target.nunique()} host stars, {int(exo.public.sum())} public visits, "
      f"~{(exo.public * 2.4).sum():.0f} GB")

multi.to_csv("data/processed/mast_g395_multivisit_targets.csv", index=False)
print("-> data/processed/mast_g395_multivisit_targets.csv")
