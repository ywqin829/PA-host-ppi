suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(FALSE)
script_path <- sub("--file=", "", grep("--file=", args, value = TRUE)[1])
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
results_dir <- file.path(project_root, "data", "v9", "reproduced")
fig_dir <- file.path(project_root, "figures", "reproduced")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

tost <- read.delim(file.path(results_dir, "gse287469_tost_corrected_df12.tsv"), check.names = FALSE)
tost$label <- paste(tost$gene, tost$celltype, sep = " ")
tost$label <- factor(tost$label, levels = rev(tost$label))
tost$status <- ifelse(tost$equivalence_established_1_0, "equivalent", "not established")

plot <- ggplot(tost, aes(WT_vs_dFliC_lfc, label)) +
  geom_vline(xintercept = c(-1, -0.5, 0, 0.5, 1), linetype = c(3, 2, 1, 2, 3), colour = "grey70") +
  geom_segment(
    aes(x = ci_lo, xend = ci_hi, y = label, yend = label, colour = status),
    linewidth = 0.8
  ) +
  geom_point(aes(colour = status), size = 2.8) +
  scale_colour_manual(
    values = c("not established" = "#C65D3A", "equivalent" = "#2C7FB8"),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(-3.5, 6), clip = "off") +
  labs(
    x = "WT vs dFliC log2FC contrast (limma-voom; df=12)",
    y = NULL,
    title = "F  WT-vs-dFliC equivalence",
    subtitle = "Gray lines: ±0.5 and ±1.0 log2FC margins; TLR5 equivalence not established"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.margin = margin(5, 18, 5, 5),
    axis.text.y = element_text(size = 7)
  )

ggsave(
  file.path(fig_dir, "Fig4F_updated_v9.pdf"),
  plot,
  width = 4.2,
  height = 3.1,
  device = "pdf",
  useDingbats = FALSE
)
ggsave(
  file.path(fig_dir, "Fig4F_updated_v9.png"),
  plot,
  width = 4.2,
  height = 3.1,
  dpi = 350
)
cat("Wrote Fig4F_updated_v9\n")
