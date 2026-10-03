# ============================================================
# Fig5_meta.R  (final v6 figures)
# Cross-study random-effects meta-analysis (REML primary)
# ?????? setwd() ?? figure_reproduction ???
# ?????? PDF??????????
# ============================================================
setwd("C:/Users/yw_qi/Downloads/pa-host-ppi/figure_reproduction")  # <- ???????
ED  <- file.path(getwd(), "data_editorial")
FIG <- file.path(getwd(), "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(ED)) stop("??? data_editorial/ ??????? setwd() ??")

read_tab <- function(f) {
  if (!file.exists(f)) stop("???????: ", f)
  read.delim(f, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})
theme_clean <- theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold", colour = "black"),
        panel.border = element_blank(),
        axis.text = element_text(size = 12), axis.title = element_text(size = 13),
        legend.text = element_text(size = 11), legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))
theme_set(theme_clean)

# ---------- ?? ----------
d   <- read_tab(file.path(ED, "meta_input_tlr5_6ds.tsv"))      # 6 independent datasets
prim<- read_tab(file.path(ED, "meta_primary_reml.tsv"))        # REML primary
st  <- read_tab(file.path(ED, "meta_strata_tlr5_6ds_reml.tsv"))
lo  <- read_tab(file.path(ED, "meta_loo_tlr5_6ds_reml.tsv"))
sens<- read_tab(file.path(ED, "meta_sensitivity_summary.tsv"))
pr  <- read_tab(file.path(ED, "meta_prr_panel_collapsed.tsv"))

# ---------- A: forest (6 datasets, REML primary) ----------
d$lab <- sub(" \\(.*\\)$", "", d$study)
d <- d[order(d$yi), ]; d$lab <- factor(d$lab, levels = d$lab)
d$w <- 1/d$sei
A <- ggplot(d, aes(yi, lab)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  geom_segment(aes(x = yi - 1.96*sei, xend = yi + 1.96*sei, y = lab, yend = lab), linewidth = 0.7, colour = "#C65D3A") +
  geom_point(aes(size = w), colour = "#C65D3A", alpha = 0.9) +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  annotate("segment", x = sens$pi_lo[sens$analysis == "Primary: REML, 6 datasets"],
           xend = sens$pi_hi[sens$analysis == "Primary: REML, 6 datasets"], y = -0.55, yend = -0.55,
           linewidth = 0.8, colour = "#8E2A18") +
  annotate("segment", x = sens$pi_lo[sens$analysis == "Primary: REML, 6 datasets"],
           xend = sens$pi_lo[sens$analysis == "Primary: REML, 6 datasets"], y = -0.62, yend = -0.48,
           linewidth = 0.8, colour = "#8E2A18") +
  annotate("segment", x = sens$pi_hi[sens$analysis == "Primary: REML, 6 datasets"],
           xend = sens$pi_hi[sens$analysis == "Primary: REML, 6 datasets"], y = -0.62, yend = -0.48,
           linewidth = 0.8, colour = "#8E2A18") +
  annotate("text", x = sens$pi_hi[sens$analysis == "Primary: REML, 6 datasets"], y = -0.55, label = "95% prediction interval", size = 3.5,
           colour = "#8E2A18", hjust = 0, vjust = -0.6) +
  geom_polygon(data = data.frame(x = c(prim$ci_lo, prim$pooled, prim$ci_hi, prim$pooled),
                                 y = c(0, 0.18, 0, -0.18)),
               aes(x = x, y = y), inherit.aes = FALSE, fill = "#B03A2E", colour = "#B03A2E") +
  annotate("text", x = prim$pooled, y = 0.34, label = paste0("pooled ", round(prim$pooled, 2),
           " (95% CI ", round(prim$ci_lo, 2), " to ", round(prim$ci_hi, 2), ")"),
           size = 3.7, colour = "#B03A2E", fontface = "bold") +
  coord_cartesian(xlim = c(-3.2, 2.6), ylim = c(-1.0, 6.8), clip = "off") +
  labs(x = "TLR5 log2FC (PA vs control)", y = NULL,
       title = "A  TLR5 meta-analysis (six datasets)",
       subtitle = "REML; I2 = 75.1%; prediction interval crosses 0")

# ---------- B: PRR panel ----------
pr <- pr[order(pr$pooled), ]
pr$gene <- factor(pr$gene, levels = pr$gene)
pr$col <- ifelse(pr$gene == "TLR5", "#B03A2E", "#C89B7B")
pr$star <- ifelse(pr$q_BH < 0.05, "*", "")
B <- ggplot(pr, aes(pooled, gene)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, y = gene, yend = gene), linewidth = 0.6, colour = pr$col) +
  geom_point(aes(size = k), colour = pr$col) +
  scale_size_continuous(range = c(2, 4.5), guide = "none") +
  geom_text(aes(label = sprintf("%.2f%s", pooled, star)), hjust = ifelse(pr$pooled < 0, 1.2, -0.4), size = 3.5) +
  coord_cartesian(xlim = c(-2.1, 3.6)) +
  labs(x = "Pooled log2FC (REML; * BH q<0.05)", y = NULL,
       title = "B  Fixed 13-gene receptor/adaptor panel",
       subtitle = "TLR5 is the only FDR-significant negative pooled effect")

# ---------- C: leave-one-dataset-out ----------
lo$removed <- sub(" \\(.*\\)$", "", lo$removed)
loo_levels <- rev(c("GSE287469", "GSE192528", "GSE324794", "SRP353174", "GSE334176", "SRP439193"))
lo$removed <- factor(lo$removed, levels = loo_levels)
C <- ggplot(lo, aes(pooled, removed)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  geom_vline(xintercept = prim$pooled, colour = "#B03A2E", linetype = 2) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, y = removed, yend = removed), linewidth = 0.7, colour = "#C65D3A") +
  geom_point(colour = "#C65D3A", size = 3.2) +
  annotate("text", x = prim$pooled, y = 6.55,
           label = sprintf("Full pooled estimate = %.2f", prim$pooled),
           colour = "#B03A2E", size = 3.6, hjust = 0.5) +
  coord_cartesian(xlim = c(-2.5, 0.2), ylim = c(0.4, 6.8)) +
  labs(x = "Pooled TLR5 log2FC after exclusion of the indicated dataset", y = NULL,
       title = "C  Leave-one-dataset-out analysis",
       subtitle = "All pooled estimates remain negative; 95% CIs exclude 0")

# ---------- D: estimator / dependence / SE robustness ----------
sens_d <- sens %>%
  mutate(
    method = sub("^Primary: ", "", analysis),
    method = sub("^Sensitivity: ", "", method),
    method = ifelse(method == "excluding reconstructed-SE studies (4 datasets)",
                    "REML excluding reconstructed-SE studies\n(4 datasets)", method),
    is_primary = grepl("^Primary:", analysis),
    pt_size = ifelse(is_primary, 4.0, 3.1),
    lnw = ifelse(is_primary, 1.1, 0.7),
    effect_txt = sprintf("%.2f (%.2f, %.2f)", pooled, ci_lo, ci_hi)
  ) %>%
  mutate(method = factor(method, levels = rev(method)))

D <- ggplot(sens_d, aes(pooled, method)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, y = method, yend = method, linewidth = lnw),
               colour = "#D35400", lineend = "round") +
  geom_point(aes(size = pt_size), colour = "#D35400") +
  geom_text(aes(label = effect_txt), nudge_y = 0.23, size = 3.4) +
  scale_size_identity() +
  scale_linewidth_identity() +
  coord_cartesian(xlim = c(-3.2, 0.8), ylim = c(0.4, 4.7), clip = "off") +
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 11),
        plot.margin = margin(6, 18, 10, 18)) +
  labs(x = "Pooled TLR5 log2FC", y = NULL) +
  labs(title = "D  Sensitivity analyses",
       subtitle = "All specified analyses retain a negative pooled effect")

# ---------- ???? PDF? ----------
fig5 <- (A | B) / (C | D) + plot_layout(heights = c(1.18, 1))
ggsave(file.path(FIG, "Fig5.pdf"), fig5, width = 12.5, height = 11, device = "pdf", useDingbats = FALSE)
cat("Fig5 saved ->", file.path(FIG, "Fig5.pdf"), "\n")
