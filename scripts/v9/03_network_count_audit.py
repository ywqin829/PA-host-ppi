from __future__ import annotations

from pathlib import Path

import pandas as pd


project_root = Path(__file__).resolve().parents[2]
data_dir = project_root / "data" / "v9"
out_dir = data_dir / "reproduced"
out_dir.mkdir(parents=True, exist_ok=True)

curated = pd.read_csv(data_dir / "pa_human_positives.tsv", sep="\t")
hub = pd.read_csv(
    data_dir / "hub_convergence_pa_enrichment.tsv",
    sep="\t",
)
source_sensitivity = pd.read_csv(
    data_dir / "robustness_leave_database_out.tsv",
    sep="\t",
)

curated_proteins = set(curated["human_uniprot"].dropna())
hub["is_pa"] = hub["is_pa"].astype(str).str.lower().eq("true")
pa_graph = hub[hub["is_pa"]].copy()
pa_proteins = set(pa_graph["human_uniprot"].dropna())

represented = curated_proteins & pa_proteins
absent = curated_proteins - pa_proteins

n_1_to_5 = int(((pa_graph["n_sp"] >= 1) & (pa_graph["n_sp"] <= 5)).sum())
n_6_or_more = int((pa_graph["n_sp"] >= 6).sum())
n_top_two_deciles = int((pa_graph["decile"] >= 9).sum())

summary = pd.DataFrame(
    [
        {
            "quantity": "curated unique human proteins (all evidence)",
            "n": len(curated_proteins),
            "definition": "Distinct human UniProt accessions in the 132 curated PA-human pairs.",
        },
        {
            "quantity": "curated proteins represented by >=1 non-PA bacterial partner",
            "n": len(represented),
            "definition": "Distinct curated human accessions with at least one non-PA cross-species interaction record.",
        },
        {
            "quantity": "curated proteins with no non-PA partner in the graph",
            "n": len(absent),
            "definition": "Absence from this particular non-PA graph; this is not evidence of PA-exclusive biology.",
        },
        {
            "quantity": "represented proteins linked to 1-5 non-PA species",
            "n": n_1_to_5,
            "definition": "Distinct non-PA bacterial species degree between 1 and 5.",
        },
        {
            "quantity": "represented proteins linked to >=6 non-PA species",
            "n": n_6_or_more,
            "definition": "High-recurrence set; not a biological pan-pathogen class.",
        },
        {
            "quantity": "represented proteins in the top two preassigned degree deciles",
            "n": n_top_two_deciles,
            "definition": "Used for the Fisher enrichment analysis; not a mutually exclusive target class.",
        },
    ]
)

curated_targets = (
    curated[["human_uniprot", "human_gene"]]
    .dropna(subset=["human_uniprot"])
    .assign(human_gene=lambda frame: frame["human_gene"].fillna(""))
    .groupby("human_uniprot", as_index=False)
    .agg(human_gene=("human_gene", lambda values: ";".join(sorted(set(values)))))
)
detail = (
    curated_targets
    .merge(
        pa_graph[
            [
                "human_uniprot",
                "gene",
                "n_sp",
                "n_pairs",
                "decile",
            ]
        ],
        on="human_uniprot",
        how="left",
        validate="many_to_one",
    )
)
detail["represented_in_non_pa_graph"] = detail["n_sp"].notna()
detail["degree_band"] = pd.cut(
    detail["n_sp"],
    bins=[-1, 0, 5, float("inf")],
    labels=["none", "1-5", ">=6"],
).astype(str)

summary.to_csv(out_dir / "network_count_audit.tsv", sep="\t", index=False)
detail.to_csv(
    out_dir / "network_target_degree_audit.tsv",
    sep="\t",
    index=False,
)
source_sensitivity.to_csv(
    out_dir / "network_source_sensitivity_from_prior_run.tsv",
    sep="\t",
    index=False,
)

missing_labels = pd.DataFrame(
    [
        {
            "withdrawn_label": "7 pan-pathogen hubs",
            "replacement": "7 represented proteins have non-PA degree >=6",
        },
        {
            "withdrawn_label": "18 multi-pathogen targets",
            "replacement": "No reconstructable mutually exclusive class with this definition",
        },
        {
            "withdrawn_label": "65 PA-exclusive targets",
            "replacement": "31 curated proteins are absent from this non-PA graph; absence is not biological PA exclusivity",
        },
    ]
)
missing_labels.to_csv(
    out_dir / "network_withdrawn_label_mapping.tsv",
    sep="\t",
    index=False,
)

print(summary.to_string(index=False))
print("\nWithdrawn labels:")
print(missing_labels.to_string(index=False))
