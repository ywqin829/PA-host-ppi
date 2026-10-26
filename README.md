# PA–Host PPI / PA-HIR

**A bias-audited Pseudomonas aeruginosa host-interaction resource for
context-aware analysis of infection transcriptomes**

This repository accompanies the PA-HIR manuscript. It contains the curated
PA–human interaction set, cross-species network analyses, classifier-bias
benchmarks, infection transcriptome reanalyses and the v9 evidence-correction
outputs.

> **Status:** v9 evidence-corrected release candidate. The analysis is
> computationally reproducible from public inputs, but no protein-level or
> functional validation is included.

## Current version

The current manuscript is:

- `docs/manuscript_v9_expert_review.md`
- `docs/Cover_Letter_v9_expert_review.md`
- `docs/Supplemental_Material_v9.pdf`
- `docs/Supplementary_Tables_v9.xlsx`

The earlier v6 manuscript is retained only in Git history.

## Corrected v9 results

- **Curated interactome:** 132 PA–human protein pairs covering 35 PA proteins
  and 92 human proteins.
- **Cross-species resource:** 10,491 bacteria–human pairs from 102 bacterial
  taxa.
- **Network denominator:** 61 of the 92 curated human proteins occur in the
  non-PA cross-species graph; 31 do not. Among the 61, 54 have 1–5 non-PA
  species and 7 have at least six.
- **Network enrichment:** top-two-degree-decile Fisher OR 3.47; after removing
  the dominant curation source, OR 2.67 (p=0.263).
- **Classifier bias:** ESM2 AUC is 0.777 against random negatives, 0.413
  against top-300 promiscuous-hub negatives and 0.247 against degree ≥8
  negatives. The decline is monotonic with negative-pool host degree
  (Spearman rho=−0.976).
- **TLR5 meta-analysis:** nine eligible datasets; REML log2FC −0.86
  (95% CI −1.50 to −0.22), I²=96.1%, prediction interval −2.81 to 1.09.
  Excluding both reanalysed SRP-linked datasets gives −0.57 (−1.20 to 0.06).
- **GSE97036 repeated measures:** six donors × four arsenic levels ×
  PA/control. TLR5 decreases in all 24 within-donor arsenic-level pairs.
  The 24 pairs are repeated conditions; the six donors are the independent
  biological units. The donor-level TLR5 effect is −1.224 (SE 0.147), which is
  the conservative input used in the cross-study meta-analysis.
- **TLR5/MUC1 context analysis:** transcripts encoding the previously linked
  proteins do not form a fixed co-regulated module across acute primary
  epithelium, persistent Calu-3 infection and iPSC-derived macrophages.

These are transcript-level findings. They do not establish protein abundance,
surface receptor availability or functional signalling.

## v9 analysis files

```text
data/v9/
├── meta_primary_9datasets_input.tsv
├── meta_primary_and_sensitivity_summary.tsv
├── meta_primary_leave_one_out.tsv
├── meta_input_tlr5_base7datasets.tsv
├── prr_panel_corrected_per_study.tsv
├── prr_panel_corrected_meta.tsv
├── network_count_audit.tsv
├── network_target_degree_audit.tsv
├── srp_reconstructed_tlr5_effects.tsv
├── gse287469_tost_corrected_df12.tsv
├── single_cell_tlr5_descriptive_summary.tsv
├── gse97036_design_rank_check.tsv
├── gse97036_mixed_model_average_effects.tsv
├── gse97036_fixed_donor_average_effects.tsv
├── gse97036_arsenic_interaction_lrt.tsv
├── gse97036_six_donor_level_effects_all_panel_genes.tsv
├── gse97036_six_donor_tlr5_minus_muc1.tsv
├── gse97036_six_donor_correlations.tsv
├── Table1_TLR5_multi_context_v9.tsv
├── GSE97036_Sample_Gene_Count_Table.csv.gz
├── GSE97036_series_matrix.txt.gz
└── input/                         public source matrices and sample metadata

scripts/v9/
├── 01_reconstruct_srp_effects.R
├── 02_corrected_meta_analysis.R
├── 03_network_count_audit.py
├── 04_correct_tost.R
├── 05_single_cell_descriptive.py
├── 06_prr_panel_corrected.R
├── 12_make_table1_v9.py
├── 15_gse97036_repeated_measures.R
├── 20_assemble_v9_figure_pdfs.py
├── 09_make_fig5_corrected_panels.R
├── 10_make_fig4f_tost_corrected.R
└── 11_make_fig6_state_summary.py
```

## Key corrections relative to v8

1. GSE262433 is included in the primary meta-analysis at its 120-hour endpoint;
   exclusion of persistent infection and the 24-hour endpoint are sensitivity
   analyses.
2. GSE233206/SRP439193 is reanalysed from raw counts with 3 vs 3 samples.
3. GSE192890/SRP353174 is reanalysed from deposited FPKM with 2 vs 2 samples.
4. GSE97036 is modelled with donor and arsenic structure. The three-df
   arsenic × infection test is reported, and donor means are used for
   conservative biological inference.
5. The unreconstructable 90-target and 7/18/65-class network statements are
   replaced by traceable 92/61/31 and 7/54 degree counts.
6. The single-cell cell-level OR and p-value are removed; the analysis is
   descriptive because infection and library are confounded.
7. TOST is recalculated with residual df=12.
8. Egger regression is not used for inference with fewer than ten studies.

## Reproduction

The v9 scripts use the compact source inputs in `data/v9/input/`. Public GEO and
SRA accession numbers and original publications are cited in the manuscript.

Run from the repository root:

```powershell
Rscript scripts/v9/01_reconstruct_srp_effects.R
Rscript scripts/v9/15_gse97036_repeated_measures.R
Rscript scripts/v9/02_corrected_meta_analysis.R
python scripts/v9/03_network_count_audit.py
Rscript scripts/v9/04_correct_tost.R
python scripts/v9/05_single_cell_descriptive.py
Rscript scripts/v9/06_prr_panel_corrected.R
python scripts/v9/12_make_table1_v9.py
Rscript scripts/v9/09_make_fig5_corrected_panels.R
Rscript scripts/v9/10_make_fig4f_tost_corrected.R
python scripts/v9/11_make_fig6_state_summary.py
python scripts/v9/20_assemble_v9_figure_pdfs.py
```

Corrected composite figure PDFs and previews are in `figures/final/`.

The current input copies are public datasets. Required R packages include
`DESeq2`, `glmmTMB`, `emmeans`, `GEOquery`, `lme4`, `limma`, `metafor`,
`ggplot2`, `dplyr`, `patchwork`, `ggsignif` and `ggrepel`. Source datasets
remain subject to their original repository terms.

## Release status

The repository now contains the v9 analysis inputs, outputs and scripts.
The existing Zenodo record `10.5281/zenodo.22025912` is the prior v1.1
release and is not represented here as a version-matched v9 archive. A new
Zenodo version must be created from this release candidate before a manuscript
cites it as the v9 data archive.

## License

The Zenodo data release is licensed under CC BY 4.0. The software in this
repository is released under the MIT License.

## Contact

Yongwei Qin

Nantong University

yw_qin@ntu.edu.cn
