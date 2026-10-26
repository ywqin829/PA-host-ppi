suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(metafor)
})

args <- commandArgs(FALSE)
script_path <- sub("--file=", "", grep("--file=", args, value = TRUE)[1])
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
data_dir <- file.path(project_root, "data", "v9")
out_dir <- file.path(data_dir, "reproduced")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel <- data.frame(
  human = c(
    "TLR2", "TLR4", "TLR5", "TLR6", "TLR9", "NOD1", "NOD2",
    "NLRC4", "NLRP3", "MYD88", "TIRAP", "IRAK1", "IRAK4"
  ),
  mouse = c(
    "Tlr2", "Tlr4", "Tlr5", "Tlr6", "Tlr9", "Nod1", "Nod2",
    "Nlrc4", "Nlrp3", "Myd88", "Tirap", "Irak1", "Irak4"
  ),
  stringsAsFactors = FALSE
)

effect_frame <- function(dataset, gene, log2FC, SE) {
  data.frame(
    dataset = dataset,
    gene = gene,
    yi = log2FC,
    sei = SE,
    stringsAsFactors = FALSE
  )
}

# GSE233206: raw mouse lung counts, three replicates/group.
counts_233 <- read.delim(
  file.path(data_dir, "GSE233206_Tissue_RNAseq_rawCounts.txt.gz"),
  check.names = FALSE
)
rownames(counts_233) <- counts_233$Geneid
counts_233 <- as.matrix(counts_233[, c("TC1", "TC2", "TC3", "TP1", "TP2", "TP3")])
mode(counts_233) <- "numeric"
metadata_233 <- data.frame(
  condition = factor(c("CTRL", "CTRL", "CTRL", "PA", "PA", "PA"), levels = c("CTRL", "PA")),
  row.names = colnames(counts_233)
)
dds_233 <- DESeqDataSetFromMatrix(
  counts_233[rowSums(counts_233 >= 10) >= 3, , drop = FALSE],
  metadata_233,
  design = ~ condition
)
dds_233 <- DESeq(dds_233, quiet = TRUE)
res_233 <- as.data.frame(results(dds_233, contrast = c("condition", "PA", "CTRL")))
res_233$ensembl <- rownames(res_233)
mouse_map <- read.delim(file.path(data_dir, "mouse_ensembl_map.tsv"), check.names = FALSE)
res_233 <- merge(
  res_233,
  mouse_map,
  by.x = "ensembl",
  by.y = "ensembl_gene_id",
  all.x = FALSE,
  all.y = FALSE
)

# GSE192890: deposited bulk FPKM, two mice/group.
fpkm <- read.delim(file.path(data_dir, "GSE192890_gene_fpkm.txt.gz"), check.names = FALSE)
expression <- as.matrix(fpkm[, c("CTRL_1", "CTRL_2", "AI_1", "AI_2")])
mode(expression) <- "numeric"
rownames(expression) <- fpkm$gene_id
group <- factor(c("CTRL", "CTRL", "ACUTE", "ACUTE"), levels = c("CTRL", "ACUTE"))
design <- model.matrix(~ group)
fit_192 <- eBayes(lmFit(log2(expression + 1), design))
res_192 <- topTable(fit_192, coef = "groupACUTE", number = Inf, sort.by = "none")
res_192$gene <- fpkm$gene_name
res_192$SE <- res_192$logFC / res_192$t

# Existing independent count-model effects.
res_324 <- read.delim(file.path(data_dir, "gse324794_ttest.tsv"), check.names = FALSE)

counts_528 <- read.delim(
  file.path(data_dir, "gse192528_counts.tsv"),
  row.names = 1,
  check.names = FALSE
)
metadata_528 <- read.delim(
  file.path(data_dir, "gse192528_coldata.tsv"),
  row.names = 1,
  check.names = FALSE
)
metadata_528$patient <- factor(metadata_528$patient)
metadata_528$pathogen <- factor(
  metadata_528$pathogen,
  levels = c("CTRL", "PA", "IAV", "BCG")
)
dds_528 <- DESeqDataSetFromMatrix(
  as.matrix(counts_528[rowSums(counts_528 > 0) > 0, , drop = FALSE]),
  metadata_528,
  design = ~ patient + pathogen
)
dds_528 <- DESeq(dds_528, quiet = TRUE)
res_528 <- as.data.frame(results(dds_528, contrast = c("pathogen", "PA", "CTRL")))
res_528$gene <- rownames(res_528)

counts_334 <- read.delim(
  file.path(data_dir, "gse334176_counts.tsv"),
  row.names = 1,
  check.names = FALSE
)
counts_334 <- rowsum(as.matrix(counts_334), rownames(counts_334))
metadata_334 <- data.frame(
  condition = factor(c("CTRL", "CTRL", "CTRL", "PA", "PA", "PA"), levels = c("CTRL", "PA")),
  row.names = colnames(counts_334)
)
dds_334 <- DESeqDataSetFromMatrix(
  counts_334[rowSums(counts_334 >= 10) >= 2, , drop = FALSE],
  metadata_334,
  design = ~ condition
)
dds_334 <- DESeq(dds_334, quiet = TRUE)
res_334 <- as.data.frame(results(dds_334, contrast = c("condition", "PA", "CTRL")))
res_334$gene <- rownames(res_334)

res_287 <- read.delim(file.path(data_dir, "gse287469_panel_classes.tsv"), check.names = FALSE)
gse97036_effects <- read.delim(
  file.path(out_dir, "gse97036_six_donor_level_effects_all_panel_genes.tsv"),
  check.names = FALSE
)
study_rows <- list()
for (index in seq_len(nrow(panel))) {
  human_gene <- panel$human[index]
  mouse_gene <- panel$mouse[index]

  row_233 <- res_233[res_233$mgi_symbol == mouse_gene, , drop = FALSE]
  if (nrow(row_233) == 1) {
    study_rows[[length(study_rows) + 1]] <- effect_frame(
      "GSE233206 / SRP439193", human_gene,
      row_233$log2FoldChange, row_233$lfcSE
    )
  }

  row_192 <- res_192[res_192$gene == mouse_gene, , drop = FALSE]
  if (nrow(row_192) == 1) {
    study_rows[[length(study_rows) + 1]] <- effect_frame(
      "GSE192890 / SRP353174", human_gene,
      row_192$logFC, row_192$SE
    )
  }

  row_324 <- res_324[tolower(res_324$gene) == tolower(mouse_gene), , drop = FALSE]
  if (nrow(row_324) == 1) {
    study_rows[[length(study_rows) + 1]] <- effect_frame(
      "GSE324794", human_gene,
      row_324$lfc, row_324$se
    )
  }

  row_528 <- res_528[res_528$gene == human_gene, , drop = FALSE]
  if (nrow(row_528) == 1) {
    study_rows[[length(study_rows) + 1]] <- effect_frame(
      "GSE192528", human_gene,
      row_528$log2FoldChange, row_528$lfcSE
    )
  }

  row_334 <- res_334[res_334$gene == human_gene, , drop = FALSE]
  if (nrow(row_334) == 1) {
    study_rows[[length(study_rows) + 1]] <- effect_frame(
      "GSE334176", human_gene,
      row_334$log2FoldChange, row_334$lfcSE
    )
  }

  row_287 <- res_287[res_287$gene == mouse_gene, , drop = FALSE]
  if (nrow(row_287) > 0) {
    weights <- 1 / row_287$WT_vs_NT_se^2
    collapsed_effect <- sum(weights * row_287$WT_vs_NT_lfc) / sum(weights)
    collapsed_se <- sqrt(1 / sum(weights))
    study_rows[[length(study_rows) + 1]] <- effect_frame(
      "GSE287469", human_gene,
      collapsed_effect, collapsed_se
    )
  }

  row_97036 <- gse97036_effects[gse97036_effects$gene == human_gene, , drop = FALSE]
  if (nrow(row_97036) == 1) {
    study_rows[[length(study_rows) + 1]] <- effect_frame(
      "GSE97036", human_gene,
      row_97036$mean_log2FC, row_97036$SE
    )
  }
}

per_study <- do.call(rbind, study_rows)
write.table(
  per_study,
  file.path(out_dir, "prr_panel_corrected_per_study.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

meta_rows <- lapply(panel$human, function(human_gene) {
  data <- per_study[per_study$gene == human_gene, , drop = FALSE]
  if (nrow(data) < 2) {
    return(NULL)
  }
  fit <- rma(yi = yi, sei = sei, data = data, method = "REML")
  data.frame(
    gene = human_gene,
    k = fit$k,
    pooled = fit$b[1],
    se = fit$se,
    p = fit$pval,
    ci_lo = fit$ci.lb,
    ci_hi = fit$ci.ub,
    I2 = fit$I2,
    tau2 = fit$tau2,
    stringsAsFactors = FALSE
  )
})
meta_result <- do.call(rbind, meta_rows)
meta_result$q_BH <- p.adjust(meta_result$p, method = "BH")
meta_result$significant_FDR05 <- meta_result$q_BH < 0.05
meta_result <- meta_result[order(meta_result$pooled), ]
meta_result$rank_negativity <- seq_len(nrow(meta_result))
write.table(
  meta_result,
  file.path(out_dir, "prr_panel_corrected_meta.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("Corrected 13-gene panel meta-analysis:\n")
print(meta_result, row.names = FALSE)
