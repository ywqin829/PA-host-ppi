suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
})

args <- commandArgs(FALSE)
script_path <- sub("--file=", "", grep("--file=", args, value = TRUE)[1])
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
data_dir <- file.path(project_root, "data", "v9")
out_dir <- file.path(data_dir, "reproduced")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

extract_one <- function(result, gene_id, dataset, model, control_n, case_n) {
  result <- as.data.frame(result)
  row <- result[rownames(result) == gene_id, , drop = FALSE]
  if (nrow(row) != 1) {
    stop(sprintf("%s: expected one row for %s, found %d", dataset, gene_id, nrow(row)))
  }
  data.frame(
    dataset = dataset,
    model = model,
    gene_id = gene_id,
    log2FC = row$log2FoldChange,
    SE = row$lfcSE,
    statistic = row$stat,
    p_value = row$pvalue,
    p_adj = row$padj,
    n_control = control_n,
    n_case = case_n,
    stringsAsFactors = FALSE
  )
}

# GSE233206 / SRP439193: raw lung RNA-seq counts, three biological replicates/group.
gse233206 <- read.delim(
  file.path(data_dir, "GSE233206_Tissue_RNAseq_rawCounts.txt.gz"),
  check.names = FALSE
)
rownames(gse233206) <- gse233206$Geneid
counts <- as.matrix(gse233206[, c("TC1", "TC2", "TC3", "TP1", "TP2", "TP3")])
mode(counts) <- "numeric"
metadata <- data.frame(
  condition = factor(
    c("CTRL", "CTRL", "CTRL", "PA", "PA", "PA"),
    levels = c("CTRL", "PA")
  ),
  row.names = colnames(counts)
)
keep <- rowSums(counts >= 10) >= 3
dds <- DESeqDataSetFromMatrix(counts[keep, , drop = FALSE], metadata, design = ~ condition)
dds <- DESeq(dds, quiet = TRUE)
gse233206_result <- results(dds, contrast = c("condition", "PA", "CTRL"))

# GSE192890 / SRP353174: deposited bulk FPKM; two mice/group.
gse192890 <- read.delim(
  file.path(data_dir, "GSE192890_gene_fpkm.txt.gz"),
  check.names = FALSE
)
expression <- as.matrix(
  gse192890[, c("CTRL_1", "CTRL_2", "AI_1", "AI_2", "CI_1", "CI_2")]
)
mode(expression) <- "numeric"
rownames(expression) <- gse192890$gene_id
log_expression <- log2(expression + 1)

group <- factor(
  c("CTRL", "CTRL", "ACUTE", "ACUTE", "CHRONIC", "CHRONIC"),
  levels = c("CTRL", "ACUTE", "CHRONIC")
)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- eBayes(lmFit(log_expression, design))
acute_contrast <- makeContrasts(ACUTE - CTRL, levels = design)
chronic_contrast <- makeContrasts(CHRONIC - CTRL, levels = design)
acute_fit <- eBayes(contrasts.fit(fit, acute_contrast))
chronic_fit <- eBayes(contrasts.fit(fit, chronic_contrast))

tidy_limma <- function(fit_object, coefficient, dataset, model, n_control, n_case) {
  table <- topTable(fit_object, coef = coefficient, number = Inf, sort.by = "none")
  row <- table["ENSMUSG00000079164", , drop = FALSE]
  data.frame(
    dataset = dataset,
    model = model,
    gene_id = "ENSMUSG00000079164",
    log2FC = row$logFC,
    SE = row$logFC / row$t,
    statistic = row$t,
    p_value = row$P.Value,
    p_adj = row$adj.P.Val,
    n_control = n_control,
    n_case = n_case,
    stringsAsFactors = FALSE
  )
}

effects <- rbind(
  extract_one(
    gse233206_result,
    "ENSMUSG00000079164",
    "GSE233206",
    "lung tissue, PA vs control, raw-count DESeq2",
    3,
    3
  ),
  tidy_limma(
    acute_fit,
    "ACUTE - CTRL",
    "GSE192890",
    "bulk lung, acute infection vs control, FPKM limma",
    2,
    2
  ),
  tidy_limma(
    chronic_fit,
    "CHRONIC - CTRL",
    "GSE192890",
    "bulk lung, chronic infection vs control, FPKM limma",
    2,
    2
  )
)

write.table(
  effects,
  file.path(out_dir, "srp_reconstructed_tlr5_effects.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("Reconstructed SRP effects:\n")
print(effects, row.names = FALSE)
