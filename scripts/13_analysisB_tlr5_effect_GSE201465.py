#!/usr/bin/env python3
"""Analysis B step 4: TLR5 effect size from GSE201465.

Design: human nasopharyngeal airway organoids from 6 donors (3 COPD, 3 non-diseased),
each donor assayed uninfected and PA-infected at 6 h -> six paired contrasts.
The count files are per-sample featureCounts-style tables inside GSE201465_RAW.tar.

Outputs:
  analysisB_GSE201465_tlr5_pairs.tsv
  analysisB_GSE201465_tlr5_effect.tsv
"""

from __future__ import annotations

import gzip
import io
import math
import os
import re
import statistics
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "analysisB_data")
TAR_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE201nnn/GSE201465/suppl/GSE201465_RAW.tar"
TAR_PATH = os.path.join(CACHE, "GSE201465_RAW.tar")

TLR5_ENSEMBL = "ENSG00000187554"


def ensure_tar() -> str:
    os.makedirs(CACHE, exist_ok=True)
    if not os.path.exists(TAR_PATH):
        with urllib.request.urlopen(TAR_URL, timeout=600) as response, open(TAR_PATH, "wb") as handle:
            handle.write(response.read())
    return TAR_PATH


def read_counts(blob: bytes) -> dict[str, float]:
    """Return gene -> count from a one- or two-column count table."""
    text = gzip.decompress(blob).decode("utf-8", "replace")
    counts: dict[str, float] = {}
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            continue
        gene = fields[0].split(".")[0]
        try:
            counts[gene] = float(fields[1])
        except ValueError:
            continue
    return counts


def main() -> None:
    tar_path = ensure_tar()
    samples: dict[str, dict[str, float]] = {}
    with tarfile.open(tar_path) as tar:
        for member in tar.getmembers():
            if not member.name.endswith(".txt.gz"):
                continue
            match = re.search(r"(NR(?:MK|PA)\d+)", member.name)
            if not match:
                continue
            samples[match.group(1)] = read_counts(tar.extractfile(member).read())

    print("samples loaded: %d -> %s" % (len(samples), sorted(samples)))
    probe = next(iter(samples))
    example = sorted(samples[probe])[:3]
    print("example gene ids:", example)

    # resolve the TLR5 identifier used by this quantification
    key = None
    for candidate in (TLR5_ENSEMBL, "TLR5"):
        for sample in samples.values():
            if candidate in sample:
                key = candidate
                break
        if key:
            break
    if key is None:
        # look for a gene-id whose row name contains TLR5
        for gene in list(samples[probe])[:50000]:
            if "TLR5" in gene.upper():
                key = gene
                break
    print("TLR5 identifier:", key)
    if key is None:
        raise SystemExit("TLR5 not found in the count tables")

    def cpm(sample: dict[str, float]) -> float:
        total = sum(sample.values())
        return (sample.get(key, 0.0) + 0.5) / (total + 1.0) * 1e6

    rows = []
    deltas = []
    for donor in ("4", "5", "6", "10", "11", "12"):
        mock, infected = samples.get("NRMK" + donor), samples.get("NRPA" + donor)
        if not mock or not infected:
            continue
        delta = math.log2(cpm(infected)) - math.log2(cpm(mock))
        deltas.append(delta)
        rows.append((donor, "COPD" if donor in ("4", "5", "6") else "non-diseased",
                     round(cpm(mock), 3), round(cpm(infected), 3), round(delta, 4)))

    with open(os.path.join(HERE, "analysisB_GSE201465_tlr5_pairs.tsv"), "w",
              encoding="utf-8") as handle:
        handle.write("donor\tdisease_state\tcontrol_CPM\tPA_CPM\tlog2FC\n")
        for row in rows:
            handle.write("\t".join(str(v) for v in row) + "\n")

    mean = statistics.mean(deltas)
    sd = statistics.stdev(deltas)
    se = sd / math.sqrt(len(deltas))
    print("paired donors: %d" % len(deltas))
    print("TLR5 log2FC per donor: %s" % ", ".join("%.2f" % d for d in deltas))
    print("mean log2FC = %.3f | SD = %.3f | SE = %.3f" % (mean, sd, se))

    with open(os.path.join(HERE, "analysisB_GSE201465_tlr5_effect.tsv"), "w",
              encoding="utf-8") as handle:
        handle.write("dataset\tidentifier\tn_paired\tdesign\tmean_log2FC\tsd\tse\n")
        handle.write("GSE201465\t%s\t%d\t%s\t%.4f\t%.4f\t%.4f\n"
                     % (key, len(deltas),
                        "human nasopharyngeal organoids, donor-paired uninfected vs PA",
                        mean, sd, se))
    print("wrote analysisB_GSE201465_tlr5_pairs.tsv and _effect.tsv")


if __name__ == "__main__":
    main()
