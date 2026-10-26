suppressPackageStartupMessages({
  library(metafor)
})

args <- commandArgs(FALSE)
script_path <- sub("--file=", "", grep("--file=", args, value = TRUE)[1])
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
data_dir <- file.path(project_root, "data", "v9")
out_dir <- file.path(data_dir, "reproduced")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

original <- read.delim(
  file.path(data_dir, "meta_input_tlr5_base7datasets.tsv"),
  check.names = FALSE
)
reconstructed <- read.delim(
  file.path(out_dir, "srp_reconstructed_tlr5_effects.tsv"),
  check.names = FALSE
)
persistent <- read.delim(
  file.path(data_dir, "gse262433_tlr5_axis.tsv"),
  check.names = FALSE
)
persistent <- persistent[persistent$gene == "TLR5", , drop = FALSE]
gse97036_effects <- read.delim(
  file.path(out_dir, "gse97036_six_donor_level_effects_all_panel_genes.tsv"),
  check.names = FALSE
)
gse97036_tlr5 <- gse97036_effects[
  gse97036_effects$gene == "TLR5",
  ,
  drop = FALSE
]

original <- original[
  !original$study %in% c("SRP439193 (mouse lung)", "SRP353174 (mouse lung)"),
  ,
  drop = FALSE
]
original$design <- "published count-model effect"
srp <- data.frame(
  study = c("GSE233206 / SRP439193", "GSE192890 / SRP353174"),
  species = "mouse",
  tissue = "lung",
  tier = 1,
  design = c("3 vs 3; raw-count DESeq2", "2 vs 2; deposited-FPKM limma"),
  yi = c(
    reconstructed$log2FC[reconstructed$dataset == "GSE233206"],
    reconstructed$log2FC[
      reconstructed$dataset == "GSE192890" &
        grepl("acute", reconstructed$model, ignore.case = TRUE)
    ]
  ),
  sei = c(
    reconstructed$SE[reconstructed$dataset == "GSE233206"],
    reconstructed$SE[
      reconstructed$dataset == "GSE192890" &
        grepl("acute", reconstructed$model, ignore.case = TRUE)
    ]
  ),
  stringsAsFactors = FALSE
)

main <- rbind(srp, original)
main <- rbind(
  main,
  data.frame(
    study = "GSE97036 primary HBE",
    species = "human",
    tissue = "airway",
    tier = 1,
    design = "6 donors x 4 arsenic levels; six donor-level means",
    yi = gse97036_tlr5$mean_log2FC,
    sei = gse97036_tlr5$SE,
    stringsAsFactors = FALSE
  )
)
main <- main[
  ,
  c("study", "species", "tissue", "tier", "design", "yi", "sei"),
  drop = FALSE
]

persistent_120 <- data.frame(
  study = "GSE262433 persistent Calu-3, 120 h",
  species = "human",
  tissue = "airway",
  tier = 2,
  design = "3 vs 3; DESeq2 at persistent endpoint",
  yi = persistent$log2FC_120H,
  sei = persistent$se_120H,
  stringsAsFactors = FALSE
)
persistent_24 <- data.frame(
  study = "GSE262433 persistent Calu-3, 24 h",
  species = "human",
  tissue = "airway",
  tier = 2,
  design = "2 vs 3; DESeq2 at 24 h",
  yi = persistent$log2FC_24H,
  sei = persistent$se_24H,
  stringsAsFactors = FALSE
)

primary_input <- rbind(main, persistent_120)
write.table(
  primary_input,
  file.path(out_dir, "meta_primary_9datasets_input.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

fit_summary <- function(data, label, test = "z") {
  fit <- rma(yi = yi, sei = sei, data = data, method = "REML", test = test)
  prediction <- if (test == "z") predict(fit) else NULL
  data.frame(
    analysis = label,
    k = fit$k,
    pooled = fit$b[1],
    ci_lo = fit$ci.lb,
    ci_hi = fit$ci.ub,
    p = fit$pval,
    I2 = fit$I2,
    tau2 = fit$tau2,
    pi_lo = if (is.null(prediction)) NA_real_ else prediction$pi.lb,
    pi_hi = if (is.null(prediction)) NA_real_ else prediction$pi.ub,
    stringsAsFactors = FALSE
  )
}

summary_rows <- list(
  fit_summary(primary_input, "Primary: all eligible datasets, persistent 120 h"),
  fit_summary(
    primary_input,
    "Primary: Hartung-Knapp inference",
    test = "knha"
  ),
  fit_summary(main, "Sensitivity: exclude persistent GSE262433"),
  fit_summary(
    primary_input[primary_input$study != "GSE97036 primary HBE", , drop = FALSE],
    "Sensitivity: exclude GSE97036"
  ),
  fit_summary(
    rbind(main, persistent_24),
    "Sensitivity: persistent GSE262433 at 24 h"
  )
)

mdsc_rows <- original[original$study == "GSE287469 (M+PMN collapsed)", , drop = FALSE]
mdsc_components <- read.delim(
  file.path(data_dir, "gse287469_flic_equiv.tsv"),
  check.names = FALSE
)
# The WT-vs-NT cell-type effects are maintained in the current panel output.
mdsc_panel <- read.delim(
  file.path(data_dir, "gse287469_panel_classes.tsv"),
  check.names = FALSE
)
mdsc_tlr5 <- mdsc_panel[mdsc_panel$gene == "Tlr5", , drop = FALSE]
if (nrow(mdsc_tlr5) == 2) {
  sensitivity_9 <- primary_input[
    primary_input$study != "GSE287469 (M+PMN collapsed)",
    ,
    drop = FALSE
  ]
  sensitivity_9 <- rbind(
    sensitivity_9,
    data.frame(
      study = paste("GSE287469", mdsc_tlr5$celltype),
      species = "mouse",
      tissue = "myeloid",
      tier = 2,
      design = "within-cell-type limma-voom",
      yi = mdsc_tlr5$WT_vs_NT_lfc,
      sei = mdsc_tlr5$WT_vs_NT_se,
      stringsAsFactors = FALSE
    )
  )
  summary_rows[[length(summary_rows) + 1]] <- fit_summary(
    sensitivity_9,
    "Sensitivity: split GSE287469 MDSC cell types"
  )
}

fit_without_srp <- rma(
  yi = yi,
  sei = sei,
  data = primary_input[
    !grepl("GSE233206|GSE192890", primary_input$study),
    ,
    drop = FALSE
  ],
  method = "REML"
)
summary_rows[[length(summary_rows) + 1]] <- data.frame(
  analysis = "Sensitivity: exclude both reanalysed SRP-linked datasets",
  k = fit_without_srp$k,
  pooled = fit_without_srp$b[1],
  ci_lo = fit_without_srp$ci.lb,
  ci_hi = fit_without_srp$ci.ub,
  p = fit_without_srp$pval,
  I2 = fit_without_srp$I2,
  tau2 = fit_without_srp$tau2,
  pi_lo = predict(fit_without_srp)$pi.lb,
  pi_hi = predict(fit_without_srp)$pi.ub,
  stringsAsFactors = FALSE
)

summary_table <- do.call(rbind, summary_rows)
write.table(
  summary_table,
  file.path(out_dir, "meta_primary_and_sensitivity_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

primary_fit <- rma(
  yi = yi,
  sei = sei,
  data = primary_input,
  method = "REML"
)
loo <- leave1out(primary_fit)
loo_table <- data.frame(
  removed = primary_input$study,
  estimate = loo$estimate,
  ci_lo = loo$ci.lb,
  ci_hi = loo$ci.ub,
  I2 = loo$I2,
  tau2 = loo$tau2,
  stringsAsFactors = FALSE
)
write.table(
  loo_table,
  file.path(out_dir, "meta_primary_leave_one_out.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("Corrected meta-analysis summary:\n")
print(summary_table, row.names = FALSE)
cat("\nLeave-one-dataset-out:\n")
print(loo_table, row.names = FALSE)
