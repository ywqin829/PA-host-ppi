#!/usr/bin/env Rscript
# GSE201465: stratify the TLR5 response by disease state.
#
# The overall donor-paired estimate is near zero, but that number is the sum of
# two opposite effects. This script fits the paired model separately in the COPD
# donors and in the non-diseased donors, and tests the interaction.
#
# Outputs:
#   analysisD_GSE201465_stratified.tsv
#   analysisD_GSE201465_interaction.tsv

suppressPackageStartupMessages({ library(DESeq2); library(metafor) })

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- getwd()

counts <- read.delim(file.path(here, "analysisB_GSE201465_counts.tsv"),
                     row.names = 1, check.names = FALSE)
smeta <- read.delim(file.path(here, "analysisB_GSE201465_samples.tsv"),
                    stringsAsFactors = FALSE)
rownames(smeta) <- smeta$sample
smeta <- smeta[colnames(counts), , drop = FALSE]
smeta$donor <- factor(smeta$donor)
smeta$condition <- factor(smeta$condition, levels = c("control", "PA"))
smeta$disease_state <- factor(smeta$disease_state,
                              levels = c("non-diseased", "COPD"))

TLR5 <- "ENSG00000187554"
fit_stratum <- function(state) {
  keep <- smeta$disease_state == state
  if (sum(keep) < 4) return(NULL)
  dds <- DESeqDataSetFromMatrix(round(as.matrix(counts[, keep])),
                                smeta[keep, , drop = FALSE],
                                design = ~ donor + condition)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  dds <- DESeq(dds, quiet = TRUE)
  r <- results(dds, contrast = c("condition", "PA", "control"))[TLR5, ]
  list(log2FC = r$log2FoldChange, SE = r$lfcSE, p = r$pvalue, padj = r$padj,
       n = sum(keep))
}

copd <- fit_stratum("COPD")
healthy <- fit_stratum("non-diseased")

cat("GSE201465 TLR5, donor-paired DESeq2 by disease state:\n")
cat(sprintf("  COPD (n=%d):         log2FC %+.3f (SE %.3f, p=%.3f)\n",
            copd$n, copd$log2FC, copd$SE, copd$p))
cat(sprintf("  non-diseased (n=%d): log2FC %+.3f (SE %.3f, p=%.3f)\n",
            healthy$n, healthy$log2FC, healthy$SE, healthy$p))

# full model with interaction: condition x disease_state
dds_all <- DESeqDataSetFromMatrix(round(as.matrix(counts)), smeta,
                                  design = ~ disease_state + condition + disease_state:condition)
dds_all <- dds_all[rowSums(counts(dds_all) >= 10) >= 3, ]
dds_all <- DESeq(dds_all, quiet = TRUE)
res_all <- results(dds_all)
coef_names <- resultsNames(dds_all)
cat("model coefficients:", paste(coef_names, collapse = " | "), "\n")
candidates <- coef_names[grepl("disease", coef_names, ignore.case = TRUE) &
                         grepl("condition", coef_names, ignore.case = TRUE)]
inter <- NULL
if (length(candidates)) {
  inter <- results(dds_all, name = candidates[1])[TLR5, ]
}

if (!is.null(inter)) {
  cat(sprintf("interaction (COPD x PA): log2FC %+.3f (SE %.3f, p=%.3f, padj=%.3f)\n",
              inter$log2FoldChange, inter$lfcSE, inter$pvalue, inter$padj))
}

strat <- data.frame(
  dataset = "GSE201465",
  stratum = c("COPD", "non-diseased", "overall (unstratified)"),
  n_donors = c(copd$n / 2, healthy$n / 2, 6),
  log2FC = c(copd$log2FC, healthy$log2FC,
             read.delim(file.path(here, "analysisB_GSE201465_deseq2.tsv"))$log2FC),
  SE = c(copd$SE, healthy$SE,
         read.delim(file.path(here, "analysisB_GSE201465_deseq2.tsv"))$SE),
  p = c(copd$p, healthy$p,
        read.delim(file.path(here, "analysisB_GSE201465_deseq2.tsv"))$pvalue)
)
write.table(strat, file.path(here, "analysisD_GSE201465_stratified.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

if (!is.null(inter)) {
  write.table(data.frame(dataset = "GSE201465",
                         test = "condition x disease_state interaction",
                         log2FC = inter$log2FoldChange, SE = inter$lfcSE,
                         p = inter$pvalue, padj = inter$padj),
              file.path(here, "analysisD_GSE201465_interaction.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("\nwrote analysisD_GSE201465_stratified.tsv",
    if (!is.null(inter)) "and analysisD_GSE201465_interaction.tsv" else "", "\n")
