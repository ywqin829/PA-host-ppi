from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy import stats


base = Path(__file__).resolve().parents[2]
v9 = base / "data" / "v9"
results = v9 / "reproduced"

meta = pd.read_csv(results / "meta_primary_9datasets_input.tsv", sep="\t")
panel = pd.read_csv(results / "prr_panel_corrected_per_study.tsv", sep="\t")
contexts = pd.read_csv(
    v9 / "combined_interface_contexts.tsv",
    sep="\t",
)
persistent = contexts[
    (contexts["dataset"] == "GSE262433")
    & (contexts["time"].astype(str) == "120 h")
].iloc[0]
ipsc = contexts[
    (contexts["dataset"] == "GSE276816")
    & (contexts["model_context"].str.contains("WT"))
    & (contexts["time"].astype(str) == "450")
].iloc[0]
gse97036 = pd.read_csv(
    results / "gse97036_six_donor_level_effects_all_panel_genes.tsv",
    sep="\t",
)
gse97036_tlr5 = gse97036[gse97036["gene"] == "TLR5"].iloc[0]

rows = []
for _, row in meta.iterrows():
    study = row["study"]
    if study == "GSE233206 / SRP439193":
        source = results / "srp_reconstructed_tlr5_effects.tsv"
        srp = pd.read_csv(source, sep="\t")
        p_value = srp[
            (srp["dataset"] == "GSE233206")
            & srp["model"].str.contains("DESeq2", regex=False)
        ].iloc[0]["p_value"]
    elif study == "GSE192890 / SRP353174":
        srp = pd.read_csv(results / "srp_reconstructed_tlr5_effects.tsv", sep="\t")
        p_value = srp[
            (srp["dataset"] == "GSE192890")
            & srp["model"].str.contains("acute", regex=False)
        ].iloc[0]["p_value"]
    elif study == "GSE324794 (mouse lung SP-AKI)":
        p_value = 0.000903
    elif study == "GSE192528 (human lung)":
        p_value = 0.091055
    elif study == "GSE334176 (human cornea)":
        p_value = 0.000469
    elif study == "GSE287469 (M+PMN collapsed)":
        p_value = 2 * stats.norm.sf(abs(row["yi"] / row["sei"]))
    elif study == "GSE201465 (human airway organoids)":
        p_value = 0.855879
    elif study == "GSE97036 primary HBE":
        p_value = gse97036_tlr5["p_value"]
    elif study == "GSE262433 persistent Calu-3, 120 h":
        p_value = 9.389468e-05
    else:
        p_value = None
    rows.append(
        {
            "study": study,
            "species": row["species"],
            "tissue_or_model": row["tissue"],
            "design": row["design"],
            "log2FC": row["yi"],
            "SE": row["sei"],
            "p_value_if_available": p_value,
            "included_in_primary_meta": True,
        }
    )

for row_index, row in enumerate(rows):
    if row["study"] == "GSE97036 primary HBE":
        rows[row_index]["design"] = (
            "6 donors x 4 arsenic levels; six donor-level means"
        )

rows.extend(
    [
        {
            "study": "GSE276816 iPSC macrophages, WT 450 min",
            "species": "human",
            "tissue_or_model": "iPSC-derived macrophages",
            "design": "single donor per genotype; descriptive",
            "log2FC": ipsc["log2FC.TLR5"],
            "SE": ipsc["SE.TLR5"],
            "p_value_if_available": ipsc["p_value.TLR5"],
            "included_in_primary_meta": False,
        }
    ]
)

table1 = pd.DataFrame(rows)
table1.to_csv(
    results / "Table1_TLR5_multi_context_v9.tsv",
    sep="\t",
    index=False,
)
print(table1.to_string(index=False))
