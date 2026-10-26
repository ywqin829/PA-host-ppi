# Cover letter

**Manuscript:** PA-HIR: a bias-audited Pseudomonas aeruginosa host-interaction resource for context-aware analysis of infection transcriptomes

**Article type:** Technical Resource

**Target journal:** Microbial Genomics

Dear Editors,

We submit our manuscript, “PA-HIR: a bias-audited Pseudomonas aeruginosa host-interaction resource for context-aware analysis of infection transcriptomes,” for consideration as a Technical Resource in Microbial Genomics.

Studies of Pseudomonas aeruginosa host interactions are limited by a small curated evidence base and by network properties that mix biological recurrence with study intensity and high-degree host-node behaviour. Infection transcriptomes provide the required direction information, but they come from different tissues, infection models, species and disease states. Researchers therefore have few reproducible ways to determine which parts of a host-target hypothesis are supported by interaction data, which are artefacts of database composition, and which depend on infection context.

PA-HIR addresses this problem by keeping interaction evidence, cross-pathogen context and infection-response evidence as separate fields. The resource combines a curated PA–human interaction map with a cross-pathogen interactome and public infection transcriptomes, and retains source, assay class, replication unit, direction and uncertainty for each layer. It also records how network enrichment changes when the dominant curation source is removed and how sequence-embedding performance changes against degree-stratified host-negative pools. These outputs allow users to query a candidate without collapsing heterogeneous evidence into a single score.

The methodological contribution challenges current practice in host–pathogen interaction benchmarking. The ESM2 classifier achieved an AUC of 0.78 against random negatives but fell to 0.41 against the top-300 high-degree host-negative pool and 0.25 against the highest-degree stratum. The apparent predictive performance was therefore highly sensitive to negative-set composition and host-node degree, while pathogen-specific interaction discrimination remains unresolved. This result has practical implications for studies that use interaction classifiers to nominate host targets.

We demonstrate the resource with TLR5 and its previously described protein-level association with MUC1, without claiming a new interaction. In six-donor primary human bronchial epithelial cultures, TLR5 decreased in all 24 infected-versus-control comparisons across four arsenic levels; the six donor-level means were all negative, and NF-κB-associated and inflammatory-response transcripts increased. Across nine eligible datasets, the average TLR5 effect was negative (−0.86, 95% CI −1.50 to −0.22) but highly heterogeneous, and exclusion of the two reanalysed SRP-linked datasets rendered the confidence interval compatible with no average effect. Persistent Calu-3 infection showed TLR5/MUC1 induction, whereas iPSC macrophages showed TLR5 attenuation with MUC1 induction. The main finding is not a universal direction of TLR5 change, but context-dependent differences in TLR5 and MUC1 transcript abundance across infection models. These observations define testable hypotheses about a previously described regulatory association without establishing a new causal mechanism.

The contribution is intended for researchers who integrate microbial interaction data, infection transcriptomes and host-directed hypotheses. The resource is distributed with versioned data tables, analysis code and a structured release directory. It is designed for reuse beyond Pseudomonas aeruginosa through the same separation of interaction, context and response evidence.

The processed data, analysis outputs and analysis code are publicly available at https://github.com/ywqin829/PA-host-ppi in release `v9.0.0`, under data/v9 and scripts/v9. The study used previously published public datasets and generated no new human or animal data. All authors have approved the manuscript, the work is not under consideration elsewhere, and the authors declare no competing interests.

We would be grateful for external peer review.

Sincerely,

Yongwei Qin  
Nantong University Medical School  
yw_qin@ntu.edu.cn  
On behalf of all authors: Jianhao Xu, Cuiqi Yang, Liangqiong Chen, Jian Zhang, Jiahui Tang, Yongwei Qin
