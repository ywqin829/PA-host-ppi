#!/usr/bin/env python3
"""Write the GSE201465 count matrix and sample sheet consumed by the DESeq2 step."""

from __future__ import annotations

import gzip
import os
import re
import tarfile

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "analysisB_data")
TAR_PATH = os.path.join(CACHE, "GSE201465_RAW.tar")


def read_counts(blob: bytes) -> dict[str, int]:
    counts: dict[str, int] = {}
    for line in gzip.decompress(blob).decode("utf-8", "replace").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            continue
        try:
            counts[fields[0].split(".")[0]] = int(float(fields[1]))
        except ValueError:
            continue
    return counts


def main() -> None:
    samples: dict[str, dict[str, int]] = {}
    with tarfile.open(TAR_PATH) as tar:
        for member in tar.getmembers():
            if not member.name.endswith(".txt.gz"):
                continue
            match = re.search(r"(NR(?:MK|PA)\d+)", member.name)
            if not match:
                continue
            samples[match.group(1)] = read_counts(tar.extractfile(member).read())

    genes = sorted(set.intersection(*(set(v) for v in samples.values())))
    order = sorted(samples, key=lambda s: (int(re.search(r"(\d+)", s).group(1)), s))

    with open(os.path.join(HERE, "analysisB_GSE201465_counts.tsv"), "w", encoding="utf-8") as handle:
        handle.write("gene\t" + "\t".join(order) + "\n")
        for gene in genes:
            handle.write(gene + "\t" + "\t".join(str(samples[s][gene]) for s in order) + "\n")

    with open(os.path.join(HERE, "analysisB_GSE201465_samples.tsv"), "w", encoding="utf-8") as handle:
        handle.write("sample\tdonor\tcondition\tdisease_state\n")
        for sample in order:
            donor = re.search(r"(\d+)", sample).group(1)
            condition = "PA" if sample.startswith("NRPA") else "control"
            state = "COPD" if donor in ("4", "5", "6") else "non-diseased"
            handle.write("%s\t%s\t%s\t%s\n" % (sample, "D" + donor, condition, state))

    print("genes: %d | samples: %d" % (len(genes), len(order)))
    print("wrote analysisB_GSE201465_counts.tsv and analysisB_GSE201465_samples.tsv")


if __name__ == "__main__":
    main()
