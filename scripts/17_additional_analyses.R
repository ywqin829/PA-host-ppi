#!/usr/bin/env Rscript
# Additional analyses for the Microbiology Spectrum resubmission.
#
#  C1  PRR/adaptor panel re-computed in the added independent human dataset (GSE201465)
#  C2  meta-analysis with the added dataset: leave-one-out, Egger regression,
#      species subgroups
#
# Outputs:
#   analysisC_GSE201465_prr_panel.tsv
#   analysisC_meta7_sensitivity.tsv
#   analysisC_meta7_loo.tsv

suppressPackageStartupMessages({
  library(DESeq2)
  library(metafor)
})

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- getwd()

panel <- c(
  TLR2   = "ENSG00000137462", TLR4  = "ENSG00000136869", TLR5  = "ENSG00000187554",
  TLR6   = "ENSG00000174130", TLR9  = "ENSG00000239732", NOD1  = "ENSG00000106100",
  NOD2   = "ENSG00000167207", NLRC4 = "ENSG00000091106", NLRP3 = "ENSG00000162711",
  MYD88  = "ENSG00000172936", TIRAP = "ENSG00000150455", IRAK1 = "ENSG00000184216",
  IRAK4  = "ENSG00000198001"
)

# ---- C1: panel in GSE201465 ---------------------------------------------------
counts <- read.delim(file.path(here, "analysisB_GSE201465_counts.tsv"),
                     row.names = 1, check.names = FALSE)
smeta <- read.delim(file.path(here, "analysisB_GSE201465_samples.tsv"), stringsAsFactors = FALSE)
rownames(smeta) <- smeta$sample
smeta <- smeta[colnames(counts), , drop = FALSE]
smeta$donor <- factor(smeta$donor)
smeta$condition <- factor(smeta$condition, levels = c("control", "PA"))

dds <- DESeqDataSetFromMatrix(round(as.matrix(counts)), smeta, design = ~ donor + condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
dds <- DESeq(dds, quiet = TRUE)

panel_res <- lapply(names(panel), function(symbol) {
  r <- results(dds, contrast = c("condition", "PA", "control"))[panel[[symbol]], ]
  data.frame(gene = symbol, ensembl = panel[[symbol]],
             log2FC = r$log2FoldChange, SE = r$lfcSE,
             pvalue = r$pvalue, padj = r$padj)
})
panel_df <- do.call(rbind, panel_res)
panel_df$q_BH <- p.adjust(panel_df$pvalue, method = "BH")
panel_df <- panel_df[order(panel_df$log2FC), ]

cat("C1  GSE201465 PRR/adaptor panel (DESeq2, donor-paired):\n")
print(panel_df, row.names = FALSE, digits = 4)
write.table(panel_df, file.path(here, "analysisC_GSE201465_prr_panel.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---- C2: meta-analysis with the added dataset ---------------------------------
meta_path <- file.path(here, "..", "data", "processed", "editorial", "meta_input_tlr5_6ds.tsv")
base <- read.delim(meta_path, stringsAsFactors = FALSE)
base$study <- trimws(base$study)
base <- base[base$study != "", ]

mdsc <- base[grepl("GSE287469", base$study), ]
if (nrow(mdsc) > 1) {
  w <- 1 / mdsc$sei^2
  base <- rbind(
    base[!grepl("GSE287469", base$study), ],
    data.frame(study = "GSE287469 (mouse MDSC)", species = "mouse",
               tissue = "myeloid", tier = 2,
               yi = sum(w * mdsc$yi) / sum(w), sei = sqrt(1 / sum(w)))
  )
}

tlr5 <- read.delim(file.path(here, "analysisB_GSE201465_deseq2.tsv"))
added <- data.frame(study = "GSE201465 (human airway organoids)", species = "human",
                    tissue = "nasopharynx", tier = 1,
                    yi = tlr5$log2FC, sei = tlr5$SE)
full <- rbind(base, added)

fit <- rma(yi = yi, sei = sei, data = full, method = "REML")

# Egger-type regression test for funnel asymmetry
egg <- regtest(fit, model = "rma")

# leave-one-out
loo <- leave1out(fit)
loo_df <- data.frame(omitted = full$study, estimate = loo$estimate,
                     ci_lb = loo$ci.lb, ci_ub = loo$ci.ub,
                     I2 = loo$I2, tau2 = loo$tau2)

# species subgroups
sub <- rma(yi = yi, sei = sei, data = full, method = "REML", subset = species == "mouse")
subh <- rma(yi = yi, sei = sei, data = full, method = "REML", subset = species == "human")
mod <- rma(yi = yi, sei = sei, mods = ~ species, data = full, method = "REML")

rows <- data.frame(
  analysis = c("7-dataset REML", "Egger test (z)", "Egger test (p)",
               "mouse subgroup", "human subgroup", "species moderator (p)"),
  value = c(sprintf("%.3f (95%% CI %.3f to %.3f); I2=%.1f%%, PI %.3f to %.3f",
                    fit$b, fit$ci.lb, fit$ci.ub, fit$I2,
                    predict(fit)$pi.lb, predict(fit)$pi.ub),
            sprintf("%.3f", egg$zval), sprintf("%.4f", egg$pval),
            sprintf("%.3f (%.3f to %.3f)", sub$b, sub$ci.lb, sub$ci.ub),
            sprintf("%.3f (%.3f to %.3f)", subh$b, subh$ci.lb, subh$ci.ub),
            sprintf("%.4f", mod$QMp))
)

cat("\nC2  meta-analysis with the added dataset:\n")
print(rows, row.names = FALSE)
cat("\nleave-one-out (7 datasets):\n")
print(loo_df, row.names = FALSE, digits = 3)

write.table(rows, file.path(here, "analysisC_meta7_sensitivity.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(loo_df, file.path(here, "analysisC_meta7_loo.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nwrote analysisC_GSE201465_prr_panel.tsv, analysisC_meta7_sensitivity.tsv, analysisC_meta7_loo.tsv\n")
