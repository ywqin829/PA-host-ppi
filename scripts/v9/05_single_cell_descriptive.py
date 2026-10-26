from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


project_root = Path(__file__).resolve().parents[2]
data_dir = project_root / "data" / "v9"
out_dir = data_dir / "reproduced"
out_dir.mkdir(parents=True, exist_ok=True)

celltype = pd.read_csv(
    data_dir / "gse274823_tlr5_celltype.tsv",
    sep="\t",
)

selected = celltype[celltype["celltype"].isin(["AlvMac", "cDC"])].copy()
selected = selected[
    [
        "celltype",
        "genotype",
        "n_sham",
        "n_PA",
        "TLR5_mean_sham",
        "TLR5_mean_PA",
        "TLR5_pct_sham",
        "TLR5_pct_PA",
    ]
]
selected = selected.rename(
    columns={
        "n_sham": "n_cells_sham",
        "n_PA": "n_cells_PA",
        "TLR5_mean_sham": "mean_lognorm_sham",
        "TLR5_mean_PA": "mean_lognorm_PA",
        "TLR5_pct_sham": "detected_percent_sham",
        "TLR5_pct_PA": "detected_percent_PA",
    }
)
selected.to_csv(
    out_dir / "single_cell_tlr5_descriptive_summary.tsv",
    sep="\t",
    index=False,
)

plot_data = selected.copy()
plot_data["cell_label"] = plot_data["celltype"].replace(
    {"AlvMac": "AlvMac", "cDC": "cDC"}
)
plot_data["WT_sham"] = np.where(
    plot_data["genotype"] == "WT",
    plot_data["detected_percent_sham"],
    np.nan,
)
plot_data["WT_PA"] = np.where(
    plot_data["genotype"] == "WT",
    plot_data["detected_percent_PA"],
    np.nan,
)
plot_data["DKO_sham"] = np.where(
    plot_data["genotype"] == "DKO",
    plot_data["detected_percent_sham"],
    np.nan,
)
plot_data["DKO_PA"] = np.where(
    plot_data["genotype"] == "DKO",
    plot_data["detected_percent_PA"],
    np.nan,
)

categories = ["AlvMac", "cDC"]
groups = [
    ("WT_sham", "WT sham", "#9ecae1"),
    ("WT_PA", "WT PA", "#08519c"),
    ("DKO_sham", "DKO sham", "#fdae6b"),
    ("DKO_PA", "DKO PA", "#a63603"),
]
values = {}
for column, _, _ in groups:
    values[column] = [
        float(
            plot_data.loc[
                plot_data["celltype"] == cell,
                column,
            ].dropna().iloc[0]
        )
        for cell in categories
    ]

fig, ax = plt.subplots(figsize=(3.35, 2.25))
x = np.arange(len(categories))
width = 0.18
for index, (column, label, color) in enumerate(groups):
    ax.bar(
        x + (index - 1.5) * width,
        values[column],
        width=width,
        label=label,
        color=color,
        edgecolor="white",
        linewidth=0.4,
    )

ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.set_ylabel("Tlr5-detected cells (%)")
ax.set_ylim(0, 14)
ax.set_title("Pooled BALF libraries", fontsize=8.5, pad=4)
ax.text(
    0.01,
    0.98,
    "One library per condition; descriptive only",
    transform=ax.transAxes,
    ha="left",
    va="top",
    fontsize=6.5,
    color="#444444",
)
ax.legend(
    frameon=False,
    fontsize=5.8,
    ncol=2,
    loc="upper right",
    handlelength=1.1,
    columnspacing=0.7,
)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="both", labelsize=7)
fig.tight_layout(pad=0.4)
fig.savefig(
    out_dir / "Fig4D_single_cell_descriptive.pdf",
    bbox_inches="tight",
)
fig.savefig(
    out_dir / "Fig4D_single_cell_descriptive.png",
    dpi=350,
    bbox_inches="tight",
)
plt.close(fig)

print(selected.to_string(index=False))
