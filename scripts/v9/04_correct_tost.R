args <- commandArgs(FALSE)
script_path <- sub("--file=", "", grep("--file=", args, value = TRUE)[1])
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
data_dir <- file.path(project_root, "data", "v9")
out_dir <- file.path(data_dir, "reproduced")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
input <- read.delim(
  file.path(data_dir, "gse287469_flic_equiv.tsv"),
  check.names = FALSE
)

# Six treatment/cell-type groups with three samples each give 12 residual df.
residual_df <- 18 - 6
alpha <- 0.05

tost_p <- function(effect, se, margin, df) {
  lower_test <- pt((effect - margin) / se, df = df, lower.tail = TRUE)
  upper_test <- pt((effect + margin) / se, df = df, lower.tail = FALSE)
  max(lower_test, upper_test)
}

result <- input[, c(
  "gene",
  "celltype",
  "WT_vs_dFliC_lfc",
  "WT_vs_dFliC_se"
)]
result$df <- residual_df
result$ci_lo <- with(
  result,
  WT_vs_dFliC_lfc - qt(0.975, residual_df) * WT_vs_dFliC_se
)
result$ci_hi <- with(
  result,
  WT_vs_dFliC_lfc + qt(0.975, residual_df) * WT_vs_dFliC_se
)
result$tost_p_margin_0_5 <- mapply(
  tost_p,
  result$WT_vs_dFliC_lfc,
  result$WT_vs_dFliC_se,
  0.5,
  residual_df
)
result$tost_p_margin_1_0 <- mapply(
  tost_p,
  result$WT_vs_dFliC_lfc,
  result$WT_vs_dFliC_se,
  1.0,
  residual_df
)
result$equivalence_established_0_5 <-
  result$tost_p_margin_0_5 < alpha
result$equivalence_established_1_0 <-
  result$tost_p_margin_1_0 < alpha
result$interpretation <- ifelse(
  result$equivalence_established_1_0,
  "Equivalence established within +/-1.0 log2FC",
  "Equivalence not established; confidence interval too wide"
)

write.table(
  result,
  file.path(out_dir, "gse287469_tost_corrected_df12.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("Corrected TOST analysis with residual df = 12:\n")
print(result, row.names = FALSE)
