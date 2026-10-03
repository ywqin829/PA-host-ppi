#!/usr/bin/env python3
"""Rebuild Table 1 (TLR5 multi-context) with the revision-added dataset.

Source of truth for the existing rows: submission_upload_v5/Table1_TLR5_multi_context.tsv
The revision adds GSE201465 (human nasopharyngeal airway organoids, donor-paired
DESeq2) as effect estimate 5, and renumbers the meta-analysis inputs to 1-8.

Outputs (TSV + XLSX) to the revision workspace and to the Zenodo revision folder.
"""

from __future__ import annotations

import os
import shutil

PROJECT = r"C:\Users\yw_qi\Downloads\pa-host-ppi"
REVISION = os.path.join(PROJECT, "revision_v7_2026-10-03")
ZENODO = os.path.join(PROJECT, "Zenodo_upload_ready",
                      "12_revision_analyses_Revision_Analyses")
ZENODO_ZH = os.path.join(PROJECT, "补充材料_Zenodo",
                         "12_返修新增分析_Revision_Analyses")

HEADER = ["Context", "Dataset", "Species/cell", "Evidence tier", "Effect estimate #",
          "TLR5 change (log2FC)", "Statistical annotation"]

ROWS = [
    ["Mouse lung SRP439193", "SRP439193", "mouse lung", "Tier 1", "1", "-2.27",
     "p=0.0006 (per-study)"],
    ["Mouse lung SRP353174", "SRP353174", "mouse lung", "Tier 1", "2", "-1.11",
     "p=0.22"],
    ["Mouse lung SP-AKI", "GSE324794", "mouse lung", "Tier 1", "3", "-0.90",
     "p=0.0009"],
    ["Human lung explants", "GSE192528", "human lung tissue", "Tier 1", "4", "-0.55",
     "p=0.090 (DESeq2, donor-blocked); 4/5 donors negative"],
    ["Human nasopharyngeal airway organoids (added in revision)", "GSE201465",
     "human airway organoids", "Tier 1", "5", "-0.016",
     "p=0.86 (DESeq2, donor-paired, ~ donor + condition, 6 donors, 6 h); added during revision"],
    ["Human corneal epithelium", "GSE334176", "human epithelial cells", "Tier 2", "6",
     "-2.08", "p=0.0005 (DESeq2)"],
    ["Mouse MDSC (M-MDSC)", "GSE287469", "mouse myeloid cells", "Tier 2", "7", "-1.41",
     "p=0.25 (limma-voom)"],
    ["Mouse MDSC (PMN-MDSC)", "GSE287469", "mouse myeloid cells", "Tier 2", "8", "-0.14",
     "p=0.81 (limma-voom); same dataset as M-MDSC (one dataset, two effect estimates)"],
    ["Human iPSC macrophages", "GSE276816", "human macrophages",
     "Tier 3 (single donor per genotype)", "-", "-2.97",
     "450 min vs 0; technical replicates; not in meta-analysis"],
    ["Calu-3 ALI (persistence)", "GSE262433", "human airway epithelial line",
     "Boundary context (normalized)", "-", "+0.90/+0.73",
     "24/120 h; DESeq2 BH q=0.0012/q=0.0023; tobramycin-maintained persistence; excluded from primary meta-analysis"],
]


def write_tsv(path: str) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(HEADER) + "\n")
        for row in ROWS:
            handle.write("\t".join(row) + "\n")


def write_xlsx(path: str) -> bool:
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Font
    except ImportError:
        print("openpyxl unavailable; xlsx not written")
        return False
    book = Workbook()
    sheet = book.active
    sheet.title = "Table1_TLR5_multi_context"
    sheet.append(HEADER)
    for cell in sheet[1]:
        cell.font = Font(bold=True)
        cell.alignment = Alignment(wrap_text=True, vertical="top")
    for row in ROWS:
        sheet.append(row)
    for column, width in zip("ABCDEFG", (46, 12, 24, 30, 16, 18, 62)):
        sheet.column_dimensions[column].width = width
    for row in sheet.iter_rows(min_row=2):
        for cell in row:
            cell.alignment = Alignment(wrap_text=True, vertical="top")
    book.save(path)
    return True


def main() -> None:
    os.makedirs(ZENODO, exist_ok=True)
    os.makedirs(ZENODO_ZH, exist_ok=True)
    tsv = os.path.join(REVISION, "Table1_TLR5_multi_context_v7.tsv")
    xlsx = os.path.join(REVISION, "Table1_TLR5_multi_context_v7.xlsx")
    write_tsv(tsv)
    ok = write_xlsx(xlsx)
    print("wrote", tsv)
    if ok:
        print("wrote", xlsx)
    for target in (ZENODO, ZENODO_ZH):
        shutil.copyfile(tsv, os.path.join(target, os.path.basename(tsv)))
        if ok:
            shutil.copyfile(xlsx, os.path.join(target, os.path.basename(xlsx)))
        print("copied to", target)
    print("rows in table:", len(ROWS))


if __name__ == "__main__":
    main()
