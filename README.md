# PA–Host PPI

**A bias-aware host–pathogen interaction framework for _Pseudomonas aeruginosa_ infection**

铜绿假单胞菌（_Pseudomonas aeruginosa_, PA）–宿主蛋白互作与感染转录组整合分析项目。该仓库用于整理数据、运行统计分析、生成论文图表，并提出可证伪的宿主感染感知假说。

> **Project status:** computational, hypothesis-generating, and under active preparation for submission. No protein-level or functional validation is included in the current release.

## Overview

This project combines:

- a curated PA–human protein–protein interaction (PPI) set;
- a cross-species bacteria–human interactome for convergence analysis and transfer benchmarking;
- bulk and single-cell infection transcriptomics across mouse and human models;
- persistence and bacterial-state analyses;
- drug–target annotation with direction-aware prioritization; and
- reproducible figure and manuscript-generation workflows.

The central question is whether PA targets host nodes that are repeatedly used by different bacterial pathogens, and whether infection-response data can identify a directionally interpretable and experimentally testable host-sensing candidate.

## Current results

The current manuscript version (`docs/manuscript_v6_submission.md`) reports:

- **Curated interactome:** 132 PA–human protein pairs covering 35 PA proteins and 92 human proteins.
- **Cross-species resource:** 10,491 curated bacteria–human pairs spanning 102 bacterial taxa.
- **Convergence:** PA targets are enriched among highly co-targeted host nodes after bias-aware null-model and robustness analyses (top-two deciles, odds ratio 3.5; direct-physical subset odds ratio 4.5).
- **Protein-language-model benchmark:** ESM2 is used as a genome-wide candidate-coverage screen, not as proof of interaction specificity. Protein-disjoint cross-validation AUC is approximately 0.744 and PA transfer AUC is approximately 0.740; performance falls against hub-like hard negatives.
- **TLR5 candidate:** a random-effects meta-analysis across six independent datasets gives a pooled TLR5 log2 fold-change of −1.18 (95% CI −1.81 to −0.54), with substantial heterogeneity (I² = 75.1%). The direction is reproducible across included acute-infection contexts but is not universal.
- **Therapeutic interpretation:** TLR5 agonism/restoration, neutrophil-elastase inhibition, and IFN-γ supplementation are retained only as direction-aware, testable hypotheses. They are not validated treatments.

## Evidence boundary

The project should be interpreted as a computational prioritization framework:

1. TLR5 transcriptional attenuation is a candidate phenotype, not evidence of reduced surface receptor abundance or impaired signaling.
2. The pooled effect is heterogeneous and context-dependent; the prediction interval crosses zero.
3. The FliC experiments do not establish FliC independence.
4. The single-cell analysis uses pooled samples and is descriptive rather than biological-replicate inference.
5. ESM2 predictions provide candidate coverage and ranking; they do not establish physical binding.
6. Drug directions require receptor-level, infection-model, toxicity, and efficacy validation.

## Repository layout

```text
data/
├── raw/                  Downloaded source data (ignored by Git)
├── processed/            Curated tables and analysis outputs (ignored by Git)
└── sources.md            Data provenance, access status, and limitations

scripts/                  R/Python data-processing and statistical-analysis scripts
figure_reproduction/      Final figure scripts and figure input tables
m8_robustness/            Leave-one-pathogen-out and robustness analyses
docs/                     Manuscripts, analysis plans, audits, and data inventory
submission_upload_v5/    Manuscript, figures, tables, and submission materials
补充材料_Zenodo/           Release package for supplementary/raw-data archiving
PLM-interact/             PPI model code used as a project dependency/reference
PROJECT_STATE.md          Detailed project state and provenance log
HANDOFF_CLAUDE.md         Extended handoff and reproducibility notes
```

Large raw and processed files are intentionally excluded from ordinary Git commits by `.gitignore`. See [`data/sources.md`](data/sources.md), [`docs/DATA_INVENTORY_FOR_REVIEW.md`](docs/DATA_INVENTORY_FOR_REVIEW.md), and the Zenodo release package for provenance and data-release details.

## Reproduction

The analyses were developed on Windows with R 4.5.1 and Python 3.13. The main computational dependencies include `DESeq2`, `limma`, `metafor`, `pandas`, `numpy`, `scipy`, `statsmodels`, and `scanpy`.

The final figures can be regenerated from the figure-reproduction package:

```powershell
cd figure_reproduction/scripts
Rscript run_all.R
```

The editorial statistical analyses are organized as follows:

```text
scripts/51_editorial_analyses.R    Donor-blocked DESeq2 and FliC analyses
scripts/52_meta_collapse.R         Collapse MDSC effects to six datasets
scripts/53_meta_reml.R              REML meta-analysis and sensitivity checks
scripts/54_editorial_mcda.py       Therapeutic-priority uncertainty analysis
scripts/56_editorial_network.py    Study-intensity-adjusted network models
```

Most analysis scripts currently contain a local project-root path near the top of the file. Update that path before running on another machine. Reproduction also requires the corresponding raw inputs and processed tables; these are not bundled into a normal source checkout because of their size and redistribution constraints.

## Data sources

The project uses public resources including UniProt, IntAct, PHI-base, BioGRID, ChEMBL, GEO, and published supplementary tables. The source status and known limitations are documented in [`data/sources.md`](data/sources.md). In particular, HoPaCI-DB and HPIDB were not used as silently assumed complete sources when their current data could not be reliably retrieved.

## Manuscript and release materials

- Main computational manuscript: [`docs/manuscript_v6_submission.md`](docs/manuscript_v6_submission.md)
- Data inventory: [`docs/DATA_INVENTORY_FOR_REVIEW.md`](docs/DATA_INVENTORY_FOR_REVIEW.md)
- Reproducibility handoff: [`docs/HANDOFF_2026-08-18.md`](docs/HANDOFF_2026-08-18.md)
- Supplementary-table index: [`docs/SUPPLEMENTARY_TABLES_INDEX.md`](docs/SUPPLEMENTARY_TABLES_INDEX.md)
- Submission package: [`submission_upload_v5/README_SUBMISSION.md`](submission_upload_v5/README_SUBMISSION.md)

## Citation

The manuscript and citation metadata will be updated after authorship, target-journal, and archival-DOI details are finalized. Until then, please cite the repository URL and identify the manuscript version used.

## License

No top-level software or data license has been selected yet. Source datasets remain subject to their original licenses and terms of use. A project-level license should be added by the maintainer before broad redistribution.

## Contact

Open an issue for reproducibility questions, data-provenance corrections, or methodological discussion.
