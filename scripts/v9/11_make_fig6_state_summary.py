from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


project = Path(__file__).resolve().parents[2]
data_dir = project / "data" / "v9"
out_dir = project / "figures" / "reproduced"
results_dir = data_dir / "reproduced"

paired = pd.read_csv(
    data_dir / "GSE97036_within_donor_paired_deltas.tsv",
    sep="\t",
)
contexts = pd.read_csv(
    data_dir / "combined_interface_contexts.tsv",
    sep="\t",
)
gse_donor_effects = pd.read_csv(
    results_dir / "gse97036_six_donor_level_effects_all_panel_genes.tsv",
    sep="\t",
).set_index("gene")

genes = ["TLR5", "MUC1", "NFKB2", "CXCL8", "IL6"]
summary_rows = []
for gene in genes:
    values = paired.loc[paired["gene"] == gene, "paired_delta_vst"].to_numpy()
    donor_means = (
        paired.loc[paired["gene"] == gene]
        .groupby("donor", as_index=False)["paired_delta_vst"]
        .mean()["paired_delta_vst"]
        .to_numpy()
    )
    mean = donor_means.mean()
    half_width = (
        stats.t.ppf(0.975, len(donor_means) - 1)
        * stats.sem(donor_means)
    )
    summary_rows.append(
        {
            "gene": gene,
            "median": np.median(values),
            "donor_mean": mean,
            "ci_low": mean - half_width,
            "ci_high": mean + half_width,
            "n_donors_down": int((donor_means < 0).sum()),
            "n_donors": len(donor_means),
            "n_pairs_down": int((values < 0).sum()),
            "n_pairs": len(values),
        }
    )
summary = pd.DataFrame(summary_rows).set_index("gene")

fig, axes = plt.subplots(
    2,
    1,
    figsize=(7.0, 8.2),
    gridspec_kw={"height_ratios": [0.85, 1.15]},
    constrained_layout=True,
)

ax = axes[0]
y_positions = np.arange(len(genes))[::-1]
for y, gene in zip(y_positions, genes):
    gene_pairs = paired.loc[paired["gene"] == gene]
    values = gene_pairs["paired_delta_vst"].to_numpy()
    donor_means = (
        gene_pairs.groupby("donor", as_index=False)["paired_delta_vst"]
        .mean()["paired_delta_vst"]
        .to_numpy()
    )
    ax.scatter(
        values,
        np.full(len(values), y),
        s=16,
        color="#B7C9D6",
        edgecolor="none",
        alpha=0.8,
        zorder=1,
    )
    ax.scatter(
        donor_means,
        np.full(len(donor_means), y),
        s=24,
        color="#8EA9BF",
        edgecolor="white",
        linewidth=0.4,
        zorder=2,
    )
    row = summary.loc[gene]
    ax.errorbar(
        row["donor_mean"],
        y,
        xerr=[
            [row["donor_mean"] - row["ci_low"]],
            [row["ci_high"] - row["donor_mean"]],
        ],
        fmt="o",
        color="#1F4E79",
        ecolor="#1F4E79",
        capsize=3,
        markersize=5,
        linewidth=1.2,
        zorder=3,
    )
    ax.text(
        3.45,
        y,
        (
            f"{int(row['n_pairs_down'])}/{int(row['n_pairs'])} pairs down\n"
            f"{int(row['n_donors_down'])}/{int(row['n_donors'])} donor means down"
        ),
        va="center",
        fontsize=6.5,
    )

ax.axvline(0, color="#666666", linewidth=0.8, linestyle="--")
ax.set_yticks(y_positions)
ax.set_yticklabels(genes)
ax.set_xlim(-1.6, 6.2)
ax.set_xlabel("Within-donor PA minus control (VST expression)")
ax.set_title(
    "A  GSE97036: six donors, four arsenic levels (24 repeated pairs)",
    loc="left",
    fontsize=10,
)
ax.text(
    0.01,
    0.02,
    (
        "Grey: 24 repeated donor–arsenic pairs; blue points: six donor means; "
        "interval: 95% t interval across the six donors"
    ),
    transform=ax.transAxes,
    fontsize=6.8,
    color="#444444",
)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

selected = contexts[
    (
        (contexts["dataset"] == "GSE192528")
        & contexts["model_context"].str.contains("PA vs uninfected", regex=False)
    )
    | (
        (contexts["dataset"] == "GSE334176")
        & contexts["model_context"].str.contains("live PA", regex=False)
    )
    | (contexts["dataset"] == "GSE262433")
    | (
        (contexts["dataset"] == "GSE276816")
        & (contexts["time"].astype(str).isin(["150", "450"]))
    )
].copy()
selected = selected[selected["dataset"] != "GSE97036"].copy()
selected = pd.concat(
    [
        selected,
        pd.DataFrame(
            [
                {
                    "dataset": "GSE97036",
                    "model_context": "six donor-level means",
                    "stimulus": "live PA",
                    "time": "6 h",
                    "log2FC.TLR5": gse_donor_effects.loc["TLR5", "mean_log2FC"],
                    "log2FC.MUC1": gse_donor_effects.loc["MUC1", "mean_log2FC"],
                }
            ]
        ),
    ],
    ignore_index=True,
)

label_map = {
    ("GSE97036", "6 h"): "Primary HBE",
    ("GSE192528", "24 h"): "Lung explant",
    ("GSE334176", "infection assay"): "Corneal epithelium",
    ("GSE262433", "24 h"): "Persistent Calu-3, 24 h",
    ("GSE262433", "120 h"): "Persistent Calu-3, 120 h",
}

ax = axes[1]
ax.axhline(0, color="#777777", linewidth=0.8)
ax.axvline(0, color="#777777", linewidth=0.8)

colors = {
    "GSE97036": "#1F77B4",
    "GSE192528": "#2CA02C",
    "GSE334176": "#D17B0F",
    "GSE262433": "#D62728",
    "GSE276816": "#17BECF",
}

for _, row in selected.iterrows():
    dataset = row["dataset"]
    ax.scatter(
        row["log2FC.TLR5"],
        row["log2FC.MUC1"],
        s=52 if dataset != "GSE276816" else 70,
        color=colors[dataset],
        edgecolor="white",
        linewidth=0.5,
    )
    key = (dataset, str(row["time"]))
    if key in label_map:
        label = label_map[key]
    elif dataset == "GSE276816":
        genotype = "WT" if "WT" in row["model_context"] else "CF"
        label = f"Macrophage {genotype}, {row['time']} min"
    else:
        label = dataset
    ax.annotate(
        label,
        (row["log2FC.TLR5"], row["log2FC.MUC1"]),
        xytext=(4, 4),
        textcoords="offset points",
        fontsize=6.5,
    )

limit = max(
    abs(selected[["log2FC.TLR5", "log2FC.MUC1"]].to_numpy()).max(),
    3.3,
) * 1.12
ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("TLR5 effect (log2 fold-change)")
ax.set_ylabel("MUC1 effect (log2 fold-change)")
ax.set_title(
    "B  TLR5 and MUC1 occupy different states rather than one fixed module",
    loc="left",
    fontsize=10,
)
ax.text(
    0.01,
    0.01,
    "Descriptive cross-model comparison; cell type, species, time and platform differ",
    transform=ax.transAxes,
    fontsize=6.8,
    color="#444444",
)

out_dir.mkdir(parents=True, exist_ok=True)
fig.savefig(out_dir / "Fig6_state_decoupling.pdf", bbox_inches="tight")
fig.savefig(out_dir / "Fig6_state_decoupling.png", dpi=350, bbox_inches="tight")
plt.close(fig)

summary.to_csv(
    results_dir / "GSE97036_fig6_paired_summary.tsv",
    sep="\t",
)
print(summary.to_string())
