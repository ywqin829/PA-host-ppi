#!/usr/bin/env python3
"""Analysis B step 3: TLR5 effect size from GSE12245 (human alveolar macrophages).

Design: paired donor samples, control vs laboratory PA strain (PAK).
Computes a per-donor log2 fold change for TLR5 and a paired mean/SE, which is
the quantity the cross-study meta-analysis consumes.

Outputs:
  analysisB_GSE12245_tlr5_pairs.tsv
  analysisB_GSE12245_tlr5_effect.tsv
"""

from __future__ import annotations

import gzip
import math
import os
import re
import statistics
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
MATRIX = os.path.join(HERE, "analysisB_data", "GSE12245_series_matrix.txt.gz")
CACHE = os.path.join(HERE, "analysisB_data")


def read_matrix(path: str) -> tuple[list[str], dict[str, list[str]], list[str]]:
    titles: list[str] = []
    rows: dict[str, list[str]] = {}
    platforms: list[str] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("!Sample_title"):
                titles = [t.strip('"') for t in line.split("\t")[1:]]
            elif line.startswith("!Series_platform_id"):
                platforms = [t.strip('"') for t in line.split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
            elif line.startswith("!series_matrix_table_end"):
                in_table = False
            elif in_table:
                fields = [f.strip('"') for f in line.split("\t")]
                if fields and fields[0] != "ID_REF":
                    rows[fields[0]] = fields[1:]
    return titles, rows, platforms


def fetch_platform_annotation(gpl: str) -> dict[str, str]:
    """Return probe -> gene symbol for a GEO platform."""
    target = os.path.join(CACHE, "%s.annot.gz" % gpl)
    digits = re.sub(r"\D", "", gpl)
    # GEO FTP groups platforms by digit count, with a minimum of three: GPLnnn / GPLnnnn / ...
    folder = "GPL" + "n" * max(3, len(digits))
    if not os.path.exists(target):
        url = ("https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/annot/%s.annot.gz"
               % (folder, gpl, gpl))
        try:
            with urllib.request.urlopen(url, timeout=300) as response, open(target, "wb") as handle:
                handle.write(response.read())
        except Exception as exc:  # noqa: BLE001
            url = ("https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/soft/%s_family.soft.gz"
                   % (folder, gpl, gpl))
            target = os.path.join(CACHE, "%s_family.soft.gz" % gpl)
            if not os.path.exists(target):
                with urllib.request.urlopen(url, timeout=600) as response, open(target, "wb") as handle:
                    handle.write(response.read())
            print("annot unavailable (%s); using SOFT family file" % exc)
    mapping: dict[str, str] = {}
    opener = gzip.open
    with opener(target, "rt", encoding="utf-8", errors="replace") as handle:
        header: list[str] = []
        for line in handle:
            if line.startswith("!platform_table_begin"):
                header = []
                continue
            if line.startswith("!platform_table_end"):
                break
            if line.startswith("ID\t") or (not header and line.startswith("ID")):
                header = [c.strip() for c in line.rstrip("\n").split("\t")]
                continue
            if not header:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            record = dict(zip(header, fields))
            symbol = (record.get("Gene symbol") or record.get("GENE_SYMBOL")
                      or record.get("gene_assignment") or "")
            if symbol:
                mapping[record["ID"]] = symbol
    return mapping


def main() -> None:
    titles, rows, platforms = read_matrix(MATRIX)
    print("samples: %d | platform: %s" % (len(titles), platforms))

    gpl = platforms[0] if platforms else "GPL96"
    mapping = fetch_platform_annotation(gpl)
    tlr5_probes = [probe for probe, symbol in mapping.items()
                   if re.search(r"\bTLR5\b", symbol)]
    print("TLR5 probes on %s: %s" % (gpl, tlr5_probes))
    if not tlr5_probes:
        raise SystemExit("no TLR5 probe found in platform annotation")

    # group samples: control vs lab strain PAK, keyed by donor label
    def donor(title: str) -> str:
        match = re.search(r"(non-smoker|smoker)\s+(\d+)", title)
        return "%s %s" % (match.group(1), match.group(2)) if match else title

    def condition(title: str) -> str:
        low = title.lower()
        if "lab strain pak" in low:
            return "PA"
        if "control" in low:
            return "control"
        return "other"

    pairs = []
    for index, title in enumerate(titles):
        if condition(title) not in ("PA", "control"):
            continue
        pairs.append((donor(title), condition(title), index, title))

    by_donor: dict[str, dict[str, int]] = {}
    for name, cond, index, _ in pairs:
        by_donor.setdefault(name, {})[cond] = index

    out_rows = []
    deltas = []
    for probe in tlr5_probes:
        values = [float(v) for v in rows[probe]]
        deltas = []
        for name, conditions in sorted(by_donor.items()):
            if "PA" not in conditions or "control" not in conditions:
                continue
            ctrl = values[conditions["control"]]
            treat = values[conditions["PA"]]
            if ctrl <= 0 or treat <= 0:
                continue
            delta = math.log2(treat) - math.log2(ctrl)
            deltas.append(delta)
            out_rows.append((probe, name, ctrl, treat, round(delta, 4)))
        if deltas:
            break

    with open(os.path.join(HERE, "analysisB_GSE12245_tlr5_pairs.tsv"), "w",
              encoding="utf-8") as handle:
        handle.write("probe\tdonor\tcontrol_value\tPA_value\tlog2FC\n")
        for row in out_rows:
            handle.write("\t".join(str(v) for v in row) + "\n")

    if not deltas:
        raise SystemExit("no usable paired samples")

    mean = statistics.mean(deltas)
    sd = statistics.stdev(deltas) if len(deltas) > 1 else float("nan")
    se = sd / math.sqrt(len(deltas)) if len(deltas) > 1 else float("nan")
    print("paired donors: %d" % len(deltas))
    print("TLR5 log2FC per donor: %s" % ", ".join("%.2f" % d for d in deltas))
    print("mean log2FC = %.3f | SD = %.3f | SE = %.3f" % (mean, sd, se))

    with open(os.path.join(HERE, "analysisB_GSE12245_tlr5_effect.tsv"), "w",
              encoding="utf-8") as handle:
        handle.write("dataset\tprobe\tn_paired\tdesign\tmean_log2FC\tsd\tse\n")
        handle.write("GSE12245\t%s\t%d\t%s\t%.4f\t%.4f\t%.4f\n"
                     % (probes[0], len(deltas),
                        "human alveolar macrophages, donor-paired control vs PAK",
                        mean, sd, se))
    print("wrote analysisB_GSE12245_tlr5_pairs.tsv and _effect.tsv")


if __name__ == "__main__":
    main()
