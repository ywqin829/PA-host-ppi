suppressPackageStartupMessages({
  library(DESeq2)
  library(glmmTMB)
  library(emmeans)
  library(GEOquery)
  library(lme4)
})

args <- commandArgs(FALSE)
script_path <- sub("--file=", "", grep("--file=", args, value = TRUE)[1])
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
data_dir <- file.path(project_root, "data", "v9")
out_dir <- file.path(data_dir, "reproduced")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.csv(
  file.path(data_dir, "GSE97036_Sample_Gene_Count_Table.csv.gz"),
  row.names = 1,
  check.names = FALSE
)
eset <- getGEO(filename = file.path(data_dir, "GSE97036_series_matrix.txt.gz"))
metadata <- pData(eset)
metadata <- metadata[match(colnames(counts), metadata$title), ]
rownames(metadata) <- colnames(counts)
metadata$condition <- factor(
  ifelse(trimws(metadata[["pseudomonas aeruginosa infection:ch1"]]) == "+", "PA", "CTRL"),
  levels = c("CTRL", "PA")
)
metadata$donor <- factor(sub("donor: ", "", metadata[["donor:ch1"]]))
metadata$arsenic <- factor(
  sub("arsenic \\(ug/l\\): ", "", metadata[["arsenic (ug/l):ch1"]])
)

# Use one common library-size normalization across all 48 samples.
dds_size <- DESeqDataSetFromMatrix(
  as.matrix(counts),
  metadata,
  design = ~ donor + arsenic + condition
)
dds_size <- estimateSizeFactors(dds_size)
size_factors <- sizeFactors(dds_size)

random_fixed_matrix <- model.matrix(~ condition * arsenic, data = metadata)
fixed_donor_matrix <- model.matrix(~ donor + condition * arsenic, data = metadata)
rank_check <- rbind(
  data.frame(
    model = "random-intercept GLMM fixed effects",
    formula = "~ condition * arsenic",
    n_samples = nrow(random_fixed_matrix),
    n_coefficients = ncol(random_fixed_matrix),
    rank = qr(random_fixed_matrix)$rank,
    full_rank = qr(random_fixed_matrix)$rank == ncol(random_fixed_matrix)
  ),
  data.frame(
    model = "fixed-donor negative-binomial sensitivity model",
    formula = "~ donor + condition * arsenic",
    n_samples = nrow(fixed_donor_matrix),
    n_coefficients = ncol(fixed_donor_matrix),
    rank = qr(fixed_donor_matrix)$rank,
    full_rank = qr(fixed_donor_matrix)$rank == ncol(fixed_donor_matrix)
  )
)

genes <- c("TLR5", "MUC1", "NFKB1", "NFKB2", "RELA", "IL6", "CXCL8")

extract_average_contrast <- function(model, gene) {
  emm <- emmeans(model, ~ condition * arsenic, type = "link")
  grid <- as.data.frame(emm@grid)
  contrast_vector <- numeric(nrow(grid))
  for (level in levels(metadata$arsenic)) {
    pa_index <- which(grid$condition == "PA" & grid$arsenic == level)
    control_index <- which(grid$condition == "CTRL" & grid$arsenic == level)
    contrast_vector[pa_index] <- 1 / length(levels(metadata$arsenic))
    contrast_vector[control_index] <- -1 / length(levels(metadata$arsenic))
  }
  contrast_result <- contrast(emm, method = list(average_PA_vs_control = contrast_vector))
  summary_result <- as.data.frame(summary(contrast_result, infer = c(TRUE, TRUE)))
  data.frame(
    gene = gene,
    log2FC = summary_result$estimate / log(2),
    SE = summary_result$SE / log(2),
    statistic = summary_result$z.ratio,
    p_value = summary_result$p.value,
    ci_low = summary_result$asymp.LCL / log(2),
    ci_high = summary_result$asymp.UCL / log(2),
    stringsAsFactors = FALSE
  )
}

mixed_results <- list()
fixed_results <- list()
interaction_results <- list()
donor_results <- list()
donor_deltas_all <- list()

for (gene in genes) {
  values <- as.numeric(counts[gene, ])
  model_data <- metadata
  model_data$count <- values
  model_data$log_size_factor <- log(size_factors[rownames(model_data)])

  full_model <- glmmTMB(
    count ~ condition * arsenic + offset(log_size_factor) + (1 | donor),
    data = model_data,
    family = nbinom2
  )
  reduced_model <- glmmTMB(
    count ~ condition + arsenic + offset(log_size_factor) + (1 | donor),
    data = model_data,
    family = nbinom2
  )

  contrast_result <- extract_average_contrast(full_model, gene)
  mixed_results[[gene]] <- contrast_result

  fixed_model <- glmmTMB(
    count ~ donor + condition * arsenic +
      offset(log_size_factor),
    data = model_data,
    family = nbinom2
  )
  fixed_results[[gene]] <- extract_average_contrast(fixed_model, gene)

  lrt <- anova(reduced_model, full_model)
  interaction_results[[gene]] <- data.frame(
    gene = gene,
    chisq = lrt$Chisq[2],
    df = lrt$`Chi Df`[2],
    p_value = lrt$`Pr(>Chisq)`[2],
    stringsAsFactors = FALSE
  )

  normalized <- counts(dds_size, normalized = TRUE)[gene, ]
  delta_rows <- list()
  for (donor in levels(metadata$donor)) {
    for (arsenic in levels(metadata$arsenic)) {
      selected <- metadata$donor == donor & metadata$arsenic == arsenic
      control_sample <- rownames(metadata)[selected & metadata$condition == "CTRL"]
      pa_sample <- rownames(metadata)[selected & metadata$condition == "PA"]
      delta_rows[[length(delta_rows) + 1]] <- data.frame(
        donor = donor,
        arsenic = arsenic,
        gene = gene,
        log2_delta = log2(normalized[pa_sample] + 1) - log2(normalized[control_sample] + 1),
        stringsAsFactors = FALSE
      )
    }
  }
  delta_table <- do.call(rbind, delta_rows)
  donor_deltas_all[[gene]] <- delta_table
  donor_means <- aggregate(log2_delta ~ donor, data = delta_table, FUN = mean)
  test <- t.test(donor_means$log2_delta, mu = 0)
  donor_results[[gene]] <- data.frame(
    gene = gene,
    n_donors = nrow(donor_means),
    mean_log2FC = mean(donor_means$log2_delta),
    SE = sd(donor_means$log2_delta) / sqrt(nrow(donor_means)),
    t = unname(test$statistic),
    df = unname(test$parameter),
    p_value = test$p.value,
    ci_low = test$conf.int[1],
    ci_high = test$conf.int[2],
    n_down = sum(donor_means$log2_delta < 0),
    stringsAsFactors = FALSE
  )
}

mixed_result <- do.call(rbind, mixed_results)
fixed_result <- do.call(rbind, fixed_results)
interaction_result <- do.call(rbind, interaction_results)
donor_result <- do.call(rbind, donor_results)
paired_deltas <- do.call(rbind, donor_deltas_all)

panel_genes <- c(
  "TLR2", "TLR4", "TLR5", "TLR6", "TLR9", "NOD1", "NOD2",
  "NLRC4", "NLRP3", "MYD88", "TIRAP", "IRAK1", "IRAK4"
)
missing_panel_genes <- setdiff(panel_genes, donor_result$gene)
additional_donor_results <- list()
for (gene in missing_panel_genes) {
  normalized <- counts(dds_size, normalized = TRUE)[gene, ]
  delta_rows <- list()
  for (donor in levels(metadata$donor)) {
    for (arsenic in levels(metadata$arsenic)) {
      selected <- metadata$donor == donor & metadata$arsenic == arsenic
      control_sample <- rownames(metadata)[selected & metadata$condition == "CTRL"]
      pa_sample <- rownames(metadata)[selected & metadata$condition == "PA"]
      delta_rows[[length(delta_rows) + 1]] <- data.frame(
        donor = donor,
        arsenic = arsenic,
        gene = gene,
        log2_delta = log2(normalized[pa_sample] + 1) - log2(normalized[control_sample] + 1),
        stringsAsFactors = FALSE
      )
    }
  }
  delta_table <- do.call(rbind, delta_rows)
  donor_means <- aggregate(log2_delta ~ donor, data = delta_table, FUN = mean)
  test <- t.test(donor_means$log2_delta, mu = 0)
  additional_donor_results[[gene]] <- data.frame(
    gene = gene,
    n_donors = nrow(donor_means),
    mean_log2FC = mean(donor_means$log2_delta),
    SE = sd(donor_means$log2_delta) / sqrt(nrow(donor_means)),
    t = unname(test$statistic),
    df = unname(test$parameter),
    p_value = test$p.value,
    ci_low = test$conf.int[1],
    ci_high = test$conf.int[2],
    n_down = sum(donor_means$log2_delta < 0),
    stringsAsFactors = FALSE
  )
}
if (length(additional_donor_results)) {
  donor_result <- rbind(donor_result, do.call(rbind, additional_donor_results))
}

wide <- reshape(
  paired_deltas,
  idvar = c("donor", "arsenic"),
  timevar = "gene",
  direction = "wide"
)
wide$inflammatory_axis <- rowMeans(
  wide[, c(
    "log2_delta.NFKB1",
    "log2_delta.NFKB2",
    "log2_delta.RELA",
    "log2_delta.IL6",
    "log2_delta.CXCL8"
  )]
)

donor_aggregate <- aggregate(
  cbind(
    log2_delta.TLR5,
    log2_delta.MUC1,
    inflammatory_axis
  ) ~ donor,
  data = wide,
  FUN = mean
)

tlr5_muc1_difference <- donor_aggregate$log2_delta.TLR5 -
  donor_aggregate$log2_delta.MUC1
tlr5_muc1_test <- t.test(tlr5_muc1_difference, mu = 0)
tlr5_muc1_summary <- data.frame(
  comparison = "six donor means: TLR5 minus MUC1",
  mean_difference = mean(tlr5_muc1_difference),
  SE = sd(tlr5_muc1_difference) / sqrt(length(tlr5_muc1_difference)),
  t = unname(tlr5_muc1_test$statistic),
  df = unname(tlr5_muc1_test$parameter),
  p_value = tlr5_muc1_test$p.value,
  ci_low = tlr5_muc1_test$conf.int[1],
  ci_high = tlr5_muc1_test$conf.int[2],
  stringsAsFactors = FALSE
)

cor_tlr5_muc1 <- cor.test(
  donor_aggregate$log2_delta.TLR5,
  donor_aggregate$log2_delta.MUC1,
  method = "spearman",
  exact = FALSE
)
cor_tlr5_inflammation <- cor.test(
  donor_aggregate$log2_delta.TLR5,
  donor_aggregate$inflammatory_axis,
  method = "spearman",
  exact = FALSE
)
correlations <- data.frame(
  comparison = c(
    "donor means: TLR5 vs MUC1",
    "donor means: TLR5 vs inflammatory transcript axis"
  ),
  n_donors = 6,
  rho = c(
    unname(cor_tlr5_muc1$estimate),
    unname(cor_tlr5_inflammation$estimate)
  ),
  p_value = c(cor_tlr5_muc1$p.value, cor_tlr5_inflammation$p.value),
  stringsAsFactors = FALSE
)

write.table(
  rank_check,
  file.path(out_dir, "gse97036_design_rank_check.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  mixed_result,
  file.path(out_dir, "gse97036_mixed_model_average_effects.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
write.table(
  fixed_result,
  file.path(out_dir, "gse97036_fixed_donor_average_effects.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
write.table(
  interaction_result,
  file.path(out_dir, "gse97036_arsenic_interaction_lrt.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
write.table(
  donor_result,
  file.path(out_dir, "gse97036_six_donor_level_effects_all_panel_genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
write.table(
  tlr5_muc1_summary,
  file.path(out_dir, "gse97036_six_donor_tlr5_minus_muc1.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
write.table(
  paired_deltas,
  file.path(out_dir, "gse97036_paired_deltas_condition_level.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
write.table(
  correlations,
  file.path(out_dir, "gse97036_six_donor_correlations.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("Design rank check:\n")
print(rank_check)
cat("\nMixed-model average infection effects:\n")
print(mixed_result, row.names = FALSE)
cat("\nFixed-donor-model average infection effects:\n")
print(fixed_result, row.names = FALSE)
cat("\nArsenic x infection LRT:\n")
print(interaction_result, row.names = FALSE)
cat("\nConservative six-donor inference:\n")
print(donor_result, row.names = FALSE)
cat("\nSix-donor TLR5-minus-MUC1 test:\n")
print(tlr5_muc1_summary, row.names = FALSE)
cat("\nSix-donor correlations:\n")
print(correlations, row.names = FALSE)
