#!/usr/bin/env python3
"""Analysis B: screen GEO for additional independent human PA infection datasets.

Purpose: the human stratum of the TLR5 meta-analysis currently contains only two
studies and its CI crosses zero. This script lists candidate human datasets with
their design metadata so a decision to add one can be made on evidence.

Outputs:
  analysisB_geo_candidates.tsv
"""

from __future__ import annotations

import json
import os
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
HERE = os.path.dirname(os.path.abspath(__file__))

SEARCH = (
    '(Pseudomonas aeruginosa[All Fields]) AND (Homo sapiens[Organism]) AND '
    '(expression profiling by high throughput sequencing[DataSet Type] OR '
    'expression profiling by array[DataSet Type])'
)

# Datasets already used in the submitted manuscript.
ALREADY_USED = {
    "GSE192528",   # human lung explants, PA/IAV/BCG
    "GSE334176",   # human corneal epithelium
    "GSE276816",   # human iPSC macrophages
    "GSE262433",   # Calu-3 persistent infection (human cell line)
}

# Known human cohort / patient datasets worth a look even if the query misses them.
EXTRA_CANDIDATES = ["GSE121871", "GSE201465", "GSE311629", "GSE325885", "GSE199424",
                    "GSE163095", "GSE218554", "GSE109118"]


def _get(endpoint: str, params: dict) -> str:
    params = dict(params)
    params.setdefault("retmode", "json")
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read().decode("utf-8", "replace")


def main() -> None:
    payload = json.loads(_get("esearch.fcgi", {"db": "gds", "term": SEARCH, "retmax": "200"}))
    ids = payload["esearchresult"]["idlist"]
    print("GEO records matching the human-PA expression query: %s"
          % payload["esearchresult"]["count"])

    rows = []
    for start in range(0, len(ids), 100):
        summary = json.loads(
            _get("esummary.fcgi", {"db": "gds", "id": ",".join(ids[start:start + 100])})
        )
        for uid in summary["result"].get("uids", []):
            record = summary["result"][uid]
            accession = record.get("accession", "")
            if not accession.startswith("GSE"):
                continue
            rows.append(
                {
                    "accession": accession,
                    "title": (record.get("title") or "").strip(),
                    "n_samples": record.get("n_samples", ""),
                    "gdstype": (record.get("gdstype") or "").strip(),
                    "taxon": (record.get("taxon") or "").strip(),
                    "summary": (record.get("summary") or "").strip(),
                    "already_used": accession in ALREADY_USED,
                }
            )
        time.sleep(0.4)

    seen = {row["accession"] for row in rows}
    for accession in EXTRA_CANDIDATES:
        if accession in seen:
            continue
        summary = json.loads(_get("esummary.fcgi", {"db": "gds", "id": accession}))
        for uid in summary["result"].get("uids", []):
            record = summary["result"][uid]
            rows.append(
                {
                    "accession": record.get("accession", accession),
                    "title": (record.get("title") or "").strip(),
                    "n_samples": record.get("n_samples", ""),
                    "gdstype": (record.get("gdstype") or "").strip(),
                    "taxon": (record.get("taxon") or "").strip(),
                    "summary": (record.get("summary") or "").strip(),
                    "already_used": False,
                }
            )
        time.sleep(0.4)

    # Heuristic priority: mentions PAO1/PA and has replication and is not used yet.
    def priority(row: dict) -> tuple:
        text = (row["title"] + " " + row["summary"]).lower()
        pa = 2 if ("pseudomonas" in text or "pao1" in text or "aeruginosa" in text) else 0
        human = 1 if "homo sapiens" in row["taxon"].lower() else 0
        try:
            samples = int(row["n_samples"])
        except (TypeError, ValueError):
            samples = 0
        replication = 1 if samples >= 6 else 0
        return (-(pa + human + replication), -samples)

    rows.sort(key=priority)

    out = os.path.join(HERE, "analysisB_geo_candidates.tsv")
    with open(out, "w", encoding="utf-8") as handle:
        handle.write("accession\tn_samples\ttaxon\talready_used\ttitle\tsummary\n")
        for row in rows:
            handle.write("\t".join(
                str(row[key]).replace("\t", " ").replace("\n", " ")
                for key in ("accession", "n_samples", "taxon", "already_used", "title", "summary")
            ) + "\n")

    print("candidates written: %d -> %s" % (len(rows), out))
    print()
    print("Top candidates not yet used:")
    for row in [r for r in rows if not r["already_used"]][:15]:
        print("  %-12s n=%-4s %s" % (row["accession"], row["n_samples"], row["title"][:88]))


if __name__ == "__main__":
    main()
