#!/usr/bin/env python3
"""Analysis B step 2: probe candidate GEO datasets for a usable PA-vs-control contrast.

Downloads file listings and, where small enough, the series matrix / supplementary
counts, then prints sample-level metadata so the group structure can be read off.

Usage:
    python 11_analysisB_fetch_candidates.py list
    python 11_analysisB_fetch_candidates.py pull GSE12245
"""

from __future__ import annotations

import gzip
import io
import os
import re
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "analysisB_data")
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/"

CANDIDATES = {
    "GSE201465": "GSE201nnn",
    "GSE12245": "GSE12nnn",
}


def listing(accession: str) -> list[str]:
    prefix = CANDIDATES[accession]
    found = []
    for sub in ("suppl/", "matrix/"):
        url = "%s%s/%s/%s" % (FTP, prefix, accession, sub)
        try:
            html = urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            found.append("!! %s -> %s" % (sub, exc))
            continue
        names = re.findall(r'href="([^"?/][^"]*)"', html)
        for name in names:
            if name.endswith("/"):
                continue
            found.append("%s%s" % (sub, name))
    return found


def download(url: str, target: str) -> str:
    os.makedirs(os.path.dirname(target), exist_ok=True)
    if not os.path.exists(target):
        with urllib.request.urlopen(url, timeout=300) as response, open(target, "wb") as handle:
            handle.write(response.read())
    return target


def show_groups(accession: str) -> None:
    prefix = CANDIDATES[accession]
    name = "%s_series_matrix.txt.gz" % accession
    url = "%s%s/%s/matrix/%s" % (FTP, prefix, accession, name)
    path = download(url, os.path.join(CACHE, name))
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title") or line.startswith("!Sample_characteristics_ch1") \
                    or line.startswith("!Sample_source_name_ch1") or line.startswith("!Sample_geo_accession"):
                print(line.rstrip()[:600])
            if line.startswith("!series_matrix_table_begin"):
                break
    print("cached:", path)


def main() -> None:
    mode = sys.argv[1] if len(sys.argv) > 1 else "list"
    if mode == "list":
        for accession in CANDIDATES:
            print("=" * 80)
            print(accession)
            for entry in listing(accession):
                print("  ", entry)
    elif mode == "pull":
        show_groups(sys.argv[2])
    else:
        raise SystemExit("unknown mode: %s" % mode)


if __name__ == "__main__":
    main()
