#!/usr/bin/env python3
"""
Download the input data. Standard library only.

    python3 -I fetch_data.py                # get the Bouma et al. (2024) age table (required)
    python3 -I fetch_data.py --refresh-nasa # also re-query the NASA Exoplanet Archive

The NASA Exoplanet Archive tables in data/ are the snapshot used in the analysis (queried on
2026-11-26). Re-querying replaces them with the current archive, so numbers can change.

The gyrochronology ages of Bouma et al. (2024) are not redistributed here. This script takes the
table table_allageinfo_agesformatted.csv from the authors' arXiv source (arXiv:2410.06246).
"""
import io
import sys
import tarfile
import urllib.parse
import urllib.request
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data"
TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
NASA_QUERIES = {
    "pscomppars_multi.csv":
        "select hostname,pl_name,pl_orbper,pl_rade,pl_bmasse,st_age,st_ageerr1,st_ageerr2,st_mass,st_met,"
        "st_teff,sy_pnum,discoverymethod,disc_facility,disc_year,pl_controv_flag from pscomppars where sy_pnum>1",
    "young_systems_archive.csv":
        "select hostname,pl_name,pl_orbper,tran_flag,discoverymethod,disc_year,disc_facility,st_age,st_ageerr1,"
        "st_ageerr2,st_age_reflink,st_teff,st_met,st_mass,sy_pnum,pl_rade,pl_controv_flag from pscomppars "
        "where sy_pnum>1 and st_age<1.0",
    "kepler_koi_cumulative.csv":
        "select kepid,kepoi_name,kepler_name,koi_disposition,koi_pdisposition,koi_period,koi_steff,koi_smet,"
        "koi_smass,koi_model_snr from cumulative where koi_pdisposition='CANDIDATE'",
}
BOUMA_EPRINT = "https://arxiv.org/e-print/2410.06246"
BOUMA_MEMBER = "table_allageinfo_agesformatted.csv"
HEADERS = {"User-Agent": "young-planets-resonance-data-fetch/1.0"}


def get(url, timeout=180):
    with urllib.request.urlopen(urllib.request.Request(url, headers=HEADERS), timeout=timeout) as r:
        return r.read()


def fetch_nasa():
    for name, query in NASA_QUERIES.items():
        url = TAP + "?" + urllib.parse.urlencode({"query": query, "format": "csv"})
        data = get(url)
        (DATA / name).write_bytes(data)
        rows = data.count(b"\n") - 1
        print(f"NASA Exoplanet Archive -> data/{name} ({rows} rows)")


def fetch_bouma():
    raw = get(BOUMA_EPRINT)
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:*") as tar:
        member = tar.getmember(BOUMA_MEMBER)          # only this file is read; nothing else is extracted
        data = tar.extractfile(member).read()
    (DATA / "bouma2024_koi_ages.csv").write_bytes(data)
    rows = data.count(b"\n") - 1
    print(f"Bouma et al. 2024 (arXiv:2410.06246) -> data/bouma2024_koi_ages.csv ({rows} rows; expected 2461)")


def main():
    DATA.mkdir(exist_ok=True)
    if "--refresh-nasa" in sys.argv or not all((DATA / n).exists() for n in NASA_QUERIES):
        fetch_nasa()
    fetch_bouma()


if __name__ == "__main__":
    main()
