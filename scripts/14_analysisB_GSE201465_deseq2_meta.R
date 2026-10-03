#!/usr/bin/env Rscript
# Analysis B step 5: DESeq2 paired model for GSE201465, then re-fit the
# cross-study TLR5 meta-analysis with the added human dataset.
#
# Inputs (written by 13_analysisB_tlr5_effect_GSE201465.py companion step):
#   analysisB_GSE201465_counts.tsv   gene x sample raw counts
#   analysisB_GSE201465_samples.tsv  sample -> donor / condition
#   ../data/processed/editorial/meta_input_tlr5_6ds.tsv  existing 6-dataset input
#
# Outputs:
#   analysisB_GSE201465_deseq2.tsv
#   analysisB_meta_with_GSE201465.tsv

suppressPackageStartupMessages({
  library(DESeq2)
  library(metafor)
})

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- getwd()

counts_path <- file.path(here, "analysisB_GSE201465_counts.tsv")
samples_path <- file.path(here, "analysisB_GSE201465_samples.tsv")
meta_path <- file.path(here, "..", "data", "processed", "editorial", "meta_input_tlr5_6ds.tsv")

counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)
meta <- read.delim(samples_path, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample
meta <- meta[colnames(counts), , drop = FALSE]

meta$donor <- factor(meta$donor)
meta$condition <- factor(meta$condition, levels = c("control", "PA"))

dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts)),
                              colData = meta,
                              design = ~ donor + condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds, quiet = TRUE)

gene <- "ENSG00000187554"  # TLR5
res <- results(dds, contrast = c("condition", "PA", "control"))
tlr5 <- res[gene, ]

cat("GSE201465 DESeq2 TLR5: log2FC =", round(tlr5$log2FoldChange, 4),
    "| SE =", round(tlr5$lfcSE, 4),
    "| p =", signif(tlr5$pvalue, 3),
    "| padj =", signif(tlr5$padj, 3), "\n")

write.table(data.frame(dataset = "GSE201465",
                       contrast = "PA vs control (donor-paired)",
                       log2FC = tlr5$log2FoldChange,
                       SE = tlr5$lfcSE,
                       pvalue = tlr5$pvalue,
                       padj = tlr5$padj),
            file.path(here, "analysisB_GSE201465_deseq2.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---- re-fit the meta-analysis -------------------------------------------------
base <- read.delim(meta_path, stringsAsFactors = FALSE)
base$study <- trimws(base$study)
base <- base[base$study != "", ]

# collapse the two GSE287469 cell types to one dataset-level estimate
mdsc <- base[grepl("GSE287469", base$study), ]
if (nrow(mdsc) > 1) {
  w <- 1 / mdsc$sei^2
  collapsed <- data.frame(study = "GSE287469 (mouse MDSC)",
                           species = "mouse", tissue = "myeloid", tier = 2,
                           yi = sum(w * mdsc$yi) / sum(w),
                           sei = sqrt(1 / sum(w)))
  base <- rbind(base[!grepl("GSE287469", base$study), ], collapsed)
}

added <- data.frame(study = "GSE201465 (human airway organoids)",
                    species = "human", tissue = "nasopharynx", tier = 1,
                    yi = tlr5$log2FoldChange, sei = tlr5$lfcSE)

fit_and_report <- function(dat, label) {
  fit <- rma(yi = yi, sei = sei, data = dat, method = "REML")
  pred <- predict(fit)
  cat(sprintf("%-28s k=%d  pooled=%.3f (95%% CI %.3f to %.3f)  I2=%.1f%%  PI %.3f to %.3f\n",
              label, nrow(dat), fit$b, fit$ci.lb, fit$ci.ub,
              fit$I2, pred$pi.lb, pred$pi.ub))
  data.frame(set = label, k = nrow(dat), pooled = fit$b,
             ci_lb = fit$ci.lb, ci_ub = fit$ci.ub, p = fit$pval,
             I2 = fit$I2, tau2 = fit$tau2,
             pi_lb = pred$pi.lb, pi_ub = pred$pi.ub)
}

rows <- rbind(
  fit_and_report(base, "6 datasets (submitted)"),
  fit_and_report(rbind(base, added), "7 datasets (with GSE201465)"),
  fit_and_report(base[base$species == "human", ], "human stratum, submitted"),
  fit_and_report(rbind(base[base$species == "human", ], added), "human stratum, +GSE201465")
)

print(rows, row.names = FALSE)
write.table(rows, file.path(here, "analysisB_meta_with_GSE201465.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nwrote analysisB_GSE201465_deseq2.tsv and analysisB_meta_with_GSE201465.tsv\n")
