suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(FALSE)
script_path <- sub("--file=", "", grep("--file=", args, value = TRUE)[1])
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
results_dir <- file.path(project_root, "data", "v9", "reproduced")
fig_dir <- file.path(project_root, "figures", "reproduced")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

input <- read.delim(file.path(results_dir, "meta_primary_9datasets_input.tsv"), check.names = FALSE)
summary <- read.delim(file.path(results_dir, "meta_primary_and_sensitivity_summary.tsv"), check.names = FALSE)
loo <- read.delim(file.path(results_dir, "meta_primary_leave_one_out.tsv"), check.names = FALSE)
panel <- read.delim(file.path(results_dir, "prr_panel_corrected_meta.tsv"), check.names = FALSE)

label_map <- c(
  "GSE233206 / SRP439193" = "Mouse lung, GSE233206",
  "GSE192890 / SRP353174" = "Mouse lung, GSE192890 (acute)",
  "GSE324794 (mouse lung SP-AKI)" = "Mouse lung, SP-AKI",
  "GSE192528 (human lung)" = "Human lung explant",
  "GSE334176 (human cornea)" = "Human cornea",
  "GSE287469 (M+PMN collapsed)" = "Mouse MDSC",
  "GSE201465 (human airway organoids)" = "Human airway organoids",
  "GSE97036 primary HBE" = "Primary HBE (GSE97036)",
  "GSE262433 persistent Calu-3, 120 h" = "Persistent Calu-3, 120 h"
)
input$label <- unname(label_map[input$study])
if (anyNA(input$label)) {
  stop("Unmapped study labels in Fig5A input.")
}
input$lo <- input$yi - qnorm(0.975) * input$sei
input$hi <- input$yi + qnorm(0.975) * input$sei

primary <- summary[summary$analysis == "Primary: all eligible datasets, persistent 120 h", ]
prediction_lo <- primary$pi_lo
prediction_hi <- primary$pi_hi

forest <- rbind(
  data.frame(
    label = input$label,
    yi = input$yi,
    lo = input$lo,
    hi = input$hi,
    row_type = "study",
    stringsAsFactors = FALSE
  ),
  data.frame(
    label = c("Pooled (REML)", "Prediction interval"),
    yi = c(primary$pooled, (prediction_lo + prediction_hi) / 2),
    lo = c(primary$ci_lo, prediction_lo),
    hi = c(primary$ci_hi, prediction_hi),
    row_type = c("pooled", "pi"),
    stringsAsFactors = FALSE
  )
)
forest$label <- factor(forest$label, levels = rev(forest$label))

plot_a <- ggplot(forest, aes(y = label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_segment(
    aes(x = lo, xend = hi, y = label, yend = label, colour = row_type),
    linewidth = 0.7
  ) +
  geom_point(aes(x = yi, colour = row_type), size = 2.3) +
  geom_text(
    aes(x = 1.17, label = sprintf("%.2f (%.2f, %.2f)", yi, lo, hi)),
    hjust = 0,
    size = 3
  ) +
  scale_colour_manual(
    values = c(study = "#C65D3A", pooled = "#8E2A18", pi = "#7C4A3A"),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(-3.4, 2.8), clip = "off") +
  labs(
    x = "TLR5 log2FC (PA vs control)",
    y = NULL,
    title = "A  TLR5 meta-analysis, nine datasets",
    subtitle = sprintf(
      "Pooled %.2f (95%% CI %.2f to %.2f); I2=%.1f%%; PI %.2f to %.2f",
      primary$pooled, primary$ci_lo, primary$ci_hi, primary$I2,
      prediction_lo, prediction_hi
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.margin = margin(5, 110, 5, 5),
    axis.text.y = element_text(size = 8.5)
  )

panel_plot <- panel[order(panel$pooled), ]
panel_plot$gene <- factor(panel_plot$gene, levels = panel_plot$gene)
panel_plot$colour <- ifelse(panel_plot$gene == "TLR5", "#B03A2E", "#B58B6A")

plot_b <- ggplot(panel_plot, aes(pooled, gene)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_segment(
    aes(x = ci_lo, xend = ci_hi, y = gene, yend = gene),
    colour = panel_plot$colour,
    linewidth = 0.65
  ) +
  geom_point(aes(size = k), colour = panel_plot$colour) +
  scale_size_continuous(range = c(1.8, 3.8), guide = "none") +
  geom_text(
    aes(label = ifelse(q_BH < 0.05, "*", "")),
    hjust = -0.5,
    size = 4
  ) +
  coord_cartesian(xlim = c(-4, 6), clip = "off") +
  labs(
    x = "Pooled log2FC (REML)",
    y = NULL,
    title = "B  Fixed receptor/adaptor panel",
    subtitle = "TLR5 is the only FDR-significant negative effect"
  ) +
  theme_classic(base_size = 9) +
  theme(plot.margin = margin(5, 25, 5, 5))

loo$removed_label <- unname(label_map[loo$removed])
loo$removed_label <- factor(loo$removed_label, levels = rev(loo$removed_label))

plot_c <- ggplot(loo, aes(estimate, removed_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_segment(
    aes(x = ci_lo, xend = ci_hi, y = removed_label, yend = removed_label),
    colour = "#C65D3A",
    linewidth = 0.7
  ) +
  geom_point(colour = "#C65D3A", size = 2.4) +
  coord_cartesian(xlim = c(-2.1, 0.35), clip = "off") +
  labs(
    x = "Pooled TLR5 log2FC with dataset excluded",
    y = NULL,
    title = "C  Leave-one-dataset-out"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 7),
    plot.margin = margin(5, 5, 5, 5)
  )

selected_analyses <- c(
  "Primary: all eligible datasets, persistent 120 h",
  "Primary: Hartung-Knapp inference",
  "Sensitivity: exclude persistent GSE262433",
  "Sensitivity: exclude GSE97036",
  "Sensitivity: persistent GSE262433 at 24 h",
  "Sensitivity: split GSE287469 MDSC cell types",
  "Sensitivity: exclude both reanalysed SRP-linked datasets"
)
short_names <- c(
  "Primary: all eligible datasets, persistent 120 h" = "Primary (9 datasets)",
  "Primary: Hartung-Knapp inference" = "Hartung-Knapp",
  "Sensitivity: exclude persistent GSE262433" = "Exclude persistent",
  "Sensitivity: exclude GSE97036" = "Exclude GSE97036",
  "Sensitivity: persistent GSE262433 at 24 h" = "Persistent endpoint 24 h",
  "Sensitivity: split GSE287469 MDSC cell types" = "MDSC effects split",
  "Sensitivity: exclude both reanalysed SRP-linked datasets" = "Exclude reanalysed SRPs"
)
sensitivity <- summary[match(selected_analyses, summary$analysis), ]
sensitivity$label <- unname(short_names[selected_analyses])
sensitivity$label <- factor(sensitivity$label, levels = rev(sensitivity$label))

plot_d <- ggplot(sensitivity, aes(pooled, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_segment(
    aes(x = ci_lo, xend = ci_hi, y = label, yend = label),
    colour = "#D35400",
    linewidth = 0.75
  ) +
  geom_point(colour = "#D35400", size = 2.5) +
  coord_cartesian(xlim = c(-2.2, 0.55), clip = "off") +
  labs(
    x = "Pooled TLR5 log2FC",
    y = NULL,
    title = "D  Sensitivity analyses",
    subtitle = "No Egger inference with fewer than 10 studies"
  ) +
  theme_classic(base_size = 9) +
  theme(plot.margin = margin(5, 20, 5, 5))

save_panel <- function(plot, name, width, height) {
  ggsave(
    file.path(fig_dir, paste0(name, ".pdf")),
    plot,
    width = width,
    height = height,
    device = "pdf",
    useDingbats = FALSE,
    limitsize = FALSE
  )
  ggsave(
    file.path(fig_dir, paste0(name, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 350,
    limitsize = FALSE
  )
}

save_panel(plot_a, "Fig5A_updated_v9", 9.2, 4.1)
save_panel(plot_b, "Fig5B_updated_v9", 4.3, 4.1)
save_panel(plot_c, "Fig5C_updated_v9", 4.3, 4.1)
save_panel(plot_d, "Fig5D_updated_v9", 4.5, 4.1)

cat("Wrote corrected Fig5 panels to", fig_dir, "\n")
