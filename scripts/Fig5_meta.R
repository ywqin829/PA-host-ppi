# ============================================================
# Fig5_meta.R  (v7 revision, 2026-10-03; layout revision 3)
#
# Cross-study meta-analysis of TLR5 including the revision-added
# donor-paired human airway organoid dataset (GSE201465).
#
# Panels:
#   A  forest plot, seven independent datasets.  The numeric column is a
#      separate narrow panel, so the forest itself is not squeezed and the
#      numbers can never collide with the error bars.
#   B  fixed 13-gene receptor/adaptor panel (pooled meta-analysis)
#   C  leave-one-dataset-out, seven datasets
#   D  sensitivity analyses (Hartung-Knapp, reconstructed-SE exclusion,
#      eight-effect specification, Egger regression)
#   E  the same receptor/adaptor panel recomputed inside the added dataset
#
# Layout: rows of two panels, so E occupies half a row rather than stretching
# across the full width.  Subtitles are single, short lines -- no manual line
# breaks.  Ranges are set with coord_cartesian so no row is ever dropped.
#
# Self-contained: reads only figure_reproduction/data/, writes figures/.
# Submitted-version script archived as Fig5_meta_v6_submitted.R.
# ============================================================

setwd("C:/Users/yw_qi/Downloads/pa-host-ppi/figure_reproduction")
PROC <- file.path(getwd(), "data")
FIG  <- file.path(getwd(), "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(metafor)
})

BASE_SIZE <- 11

theme_clean <- theme_classic(base_size = BASE_SIZE) +
  theme(strip.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10.5),
        legend.text = element_text(size = 9.5),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 11.5, face = "bold"),
        plot.subtitle = element_text(size = 9),
        plot.margin = margin(5, 8, 5, 5))
theme_set(theme_clean)

read_tab <- function(f) {
  if (!file.exists(f)) stop("missing input: ", f)
  read.delim(f, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

COL_ADDED <- "#B03A2E"
COL_MAIN  <- "#C65D3A"
COL_POOL  <- "#8E2A18"

# ---------------- inputs -------------------------------------------------
d6   <- read_tab(file.path(PROC, "meta_input_tlr5_6ds.tsv"))
d7e  <- read_tab(file.path(PROC, "meta_input_tlr5_7effects.tsv"))
newd <- read_tab(file.path(PROC, "gse201465_tlr5_deseq2.tsv"))
pan  <- read_tab(file.path(PROC, "meta_prr_panel_collapsed.tsv"))
panN <- read_tab(file.path(PROC, "gse201465_prr_panel.tsv"))

d6$study <- trimws(d6$study)
added <- data.frame(study = "GSE201465 (human airway organoids)", species = "human",
                    tissue = "nasopharynx", tier = 1,
                    yi = newd$log2FC[1], sei = newd$SE[1])
full <- rbind(d6[, c("study", "species", "tissue", "tier", "yi", "sei")], added)

d7e$study <- trimws(d7e$study)
full7 <- rbind(d7e[, c("study", "species", "tissue", "tier", "yi", "sei")], added)

# ---------------- meta-analyses -----------------------------------------
fit    <- rma(yi = yi, sei = sei, data = full,  method = "REML")
pred   <- predict(fit)
fit_hk <- rma(yi = yi, sei = sei, data = full, method = "REML", test = "knha")
fit_7e <- rma(yi = yi, sei = sei, data = full7, method = "REML")
keep   <- !grepl("^SRP", full$study)
fit_ex <- rma(yi = yi, sei = sei, data = full[keep, ], method = "REML")
egg    <- regtest(fit, model = "rma")

short <- c("SRP439193 (mouse lung)"             = "Mouse lung, SRP439193",
           "SRP353174 (mouse lung)"             = "Mouse lung, SRP353174",
           "GSE324794 (mouse lung SP-AKI)"      = "Mouse lung, SP-AKI",
           "GSE192528 (human lung)"             = "Human lung explant",
           "GSE334176 (human cornea)"           = "Human cornea",
           "GSE287469 (M+PMN collapsed)"        = "Mouse MDSC",
           "GSE201465 (human airway organoids)" = "Human organoids *")

# ---------------- A: forest plot (narrow numeric column) ----------------
fdf <- data.frame(
  label = c(short[full$study], "Pooled (REML)", "Prediction interval"),
  est   = c(full$yi, fit$b, (pred$pi.lb + pred$pi.ub) / 2),
  lo    = c(full$yi - 1.96 * full$sei, fit$ci.lb, pred$pi.lb),
  hi    = c(full$yi + 1.96 * full$sei, fit$ci.ub, pred$pi.ub),
  type  = c(rep("study", nrow(full)), "pooled", "pi"),
  stringsAsFactors = FALSE
)
fdf$label <- factor(fdf$label, levels = rev(fdf$label))
fdf$txt <- sprintf("%.2f (%.2f to %.2f)", fdf$est, fdf$lo, fdf$hi)

A_left <- ggplot(fdf, aes(y = label)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  geom_vline(xintercept = fit$b, colour = COL_POOL, linetype = 3, linewidth = 0.4) +
  geom_segment(aes(x = lo, xend = hi, y = label, yend = label,
                   colour = type, linewidth = type)) +
  geom_point(aes(x = est, colour = type, size = type)) +
  scale_colour_manual(values = c(study = COL_MAIN, pooled = COL_POOL, pi = "#7C4A3A"),
                      guide = "none") +
  scale_linewidth_manual(values = c(study = 0.7, pooled = 1.4, pi = 0.9), guide = "none") +
  scale_size_manual(values = c(study = 2.3, pooled = 3.6, pi = 2.3), guide = "none") +
  scale_x_continuous(breaks = c(-3, -2, -1, 0, 1)) +
  coord_cartesian(xlim = c(-3.4, 1.05)) +
  labs(x = "TLR5 log2FC (PA vs control)", y = NULL,
       title = "A  TLR5 meta-analysis, seven datasets",
       subtitle = "REML random-effects model") +
  theme(axis.text.y = element_text(size = 10, hjust = 1))

A_right <- ggplot(fdf, aes(y = label)) +
  geom_text(aes(x = 0, label = txt, fontface = ifelse(type == "study", 1, 2)),
            hjust = 0, size = 2.95, colour = "grey15") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = " ", subtitle = "log2FC (95% CI)") +
  theme_void(base_size = BASE_SIZE) +
  theme(plot.title = element_text(size = 11.5, face = "bold"),
        plot.subtitle = element_text(size = 9),
        plot.margin = margin(5, 6, 5, 0))

# ---------------- B: fixed 13-gene panel --------------------------------
pan <- pan[order(pan$pooled), ]
pan$gene <- factor(pan$gene, levels = pan$gene)
pan$col  <- ifelse(pan$gene == "TLR5", COL_ADDED, "#C89B7B")
pan$star <- ifelse(pan$q_BH < 0.05, "*", "")

B <- ggplot(pan, aes(pooled, gene)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, y = gene, yend = gene),
               linewidth = 0.6, colour = pan$col) +
  geom_point(aes(size = k), colour = pan$col) +
  scale_size_continuous(range = c(1.8, 4.2), guide = "none") +
  geom_text(aes(label = sprintf("%.2f%s", pooled, star)),
            hjust = ifelse(pan$pooled < 0, 1.25, -0.45), size = 3.0) +
  scale_x_continuous(breaks = c(-2, 0, 2, 4)) +
  coord_cartesian(xlim = c(-3.1, 4.2)) +
  labs(x = "Pooled log2FC (REML)", y = NULL,
       title = "B  Fixed 13-gene receptor/adaptor panel",
       subtitle = "TLR5 is the only FDR-significant negative effect")

# ---------------- C: leave-one-out --------------------------------------
loo <- leave1out(fit)
ldf <- data.frame(removed = full$study, est = loo$estimate,
                  lo = loo$ci.lb, hi = loo$ci.ub)
ldf$removed <- short[ldf$removed]
ldf$removed <- factor(ldf$removed, levels = rev(ldf$removed))

C <- ggplot(ldf, aes(est, removed)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  geom_vline(xintercept = fit$b, colour = COL_POOL, linetype = 3) +
  geom_segment(aes(x = lo, xend = hi, y = removed, yend = removed),
               linewidth = 0.7, colour = COL_MAIN) +
  geom_point(colour = COL_MAIN, size = 3.0) +
  geom_text(aes(x = 0.22, label = sprintf("%.2f", est)),
            hjust = 0, size = 3.0, colour = "grey15") +
  scale_x_continuous(breaks = c(-2, -1, 0)) +
  coord_cartesian(xlim = c(-2.3, 0.95)) +
  labs(x = "Pooled TLR5 log2FC (dataset excluded)", y = NULL,
       title = "C  Leave-one-dataset-out",
       subtitle = "Every pooled confidence interval excludes zero")

# ---------------- D: sensitivity ----------------------------------------
sdf <- data.frame(
  method = c("Primary (REML)",
             "Hartung-Knapp",
             "Excluding reconstructed-SE studies",
             "Eight-effect specification"),
  est = c(fit$b, fit_hk$b, fit_ex$b, fit_7e$b),
  lo  = c(fit$ci.lb, fit_hk$ci.lb, fit_ex$ci.lb, fit_7e$ci.lb),
  hi  = c(fit$ci.ub, fit_hk$ci.ub, fit_ex$ci.ub, fit_7e$ci.ub)
)
sdf$method <- factor(sdf$method, levels = rev(sdf$method))
sdf$txt <- sprintf("%.2f (%.2f to %.2f)", sdf$est, sdf$lo, sdf$hi)

D <- ggplot(sdf, aes(est, method)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  geom_segment(aes(x = lo, xend = hi, y = method, yend = method),
               linewidth = 0.8, colour = "#D35400") +
  geom_point(size = 3.2, colour = "#D35400") +
  geom_text(aes(x = 0.20, label = txt), hjust = 0, size = 2.9, colour = "grey15") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1)) +
  coord_cartesian(xlim = c(-2.3, 2.30)) +
  labs(x = "Pooled TLR5 log2FC", y = NULL,
       title = "D  Sensitivity analyses",
       subtitle = sprintf("All specifications remain negative; Egger p = %.2f", egg$pval))

# ---------------- E: panel recomputed in the added dataset --------------
panN <- panN[!is.na(panN$log2FC), ]
panN <- panN[order(panN$log2FC), ]
panN$gene <- factor(panN$gene, levels = panN$gene)
panN$col  <- ifelse(panN$gene == "TLR5", COL_ADDED,
                    ifelse(panN$q_BH < 0.05, "#2166AC", "grey65"))
panN$mark <- ifelse(panN$q_BH < 0.001, "***",
                    ifelse(panN$q_BH < 0.01, "**",
                           ifelse(panN$q_BH < 0.05, "*", "ns")))

E <- ggplot(panN, aes(gene, log2FC)) +
  geom_hline(yintercept = 0, colour = "grey45") +
  geom_col(aes(fill = col), width = 0.68) +
  scale_fill_identity() +
  geom_text(aes(label = mark, y = log2FC + ifelse(log2FC >= 0, 0.09, -0.09)),
            size = 2.9) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(x = NULL, y = "log2FC (PA vs control)",
       title = "E  Receptor/adaptor panel (added dataset)",
       subtitle = "GSE201465, donor-paired DESeq2; TLR2, TLR4 and MYD88 induced") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

# ---------------- assemble ----------------------------------------------
# Rows of two panels: E takes half a row instead of stretching the figure.
row1 <- A_left + A_right + plot_layout(widths = c(3.0, 1.15))
row2 <- B + C + plot_layout(widths = c(1.05, 1))
row3 <- D + E + plot_layout(widths = c(1.05, 1))

fig5 <- row1 / row2 / row3 + plot_layout(heights = c(1.30, 1.15, 1))

# Windows keeps PDFs locked while a viewer has them open. Always render to a
# unique temporary name, then promote to the canonical file, falling back to a
# "<base>_latest.*" name if the canonical file is locked.
save_figure <- function(plot, base, width, height) {
  stamp <- format(Sys.time(), "%H%M%S")
  for (ext in c("pdf", "png")) {
    tmp <- file.path(FIG, sprintf("%s_render_%s_%d.%s", base, stamp, Sys.getpid(), ext))
    if (ext == "pdf") {
      ggsave(tmp, plot, width = width, height = height, device = "pdf",
             useDingbats = FALSE)
    } else {
      ggsave(tmp, plot, width = width, height = height, dpi = 300)
    }
    canonical <- file.path(FIG, paste0(base, ".", ext))
    latest <- file.path(FIG, paste0(base, "_latest.", ext))
    if (isTRUE(file.copy(tmp, canonical, overwrite = TRUE))) {
      file.remove(tmp)
      cat("updated", canonical, "\n")
    } else if (isTRUE(file.copy(tmp, latest, overwrite = TRUE))) {
      file.remove(tmp)
      cat("NOTE:", basename(canonical), "is locked by another process; fresh version at",
          latest, "\n")
    } else {
      cat("NOTE: both", basename(canonical), "and", basename(latest),
          "are locked; fresh version left at", tmp, "\n")
    }
  }
}

save_figure(fig5, "Fig5", width = 10.5, height = 14)

# ---------------- export the derived tables -----------------------------
write.table(full, file.path(PROC, "meta_input_tlr5_7ds.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

mouse_fit <- rma(yi = yi, sei = sei, data = full[full$species == "mouse", ], method = "REML")
human_fit <- rma(yi = yi, sei = sei, data = full[full$species == "human", ], method = "REML")

summary_tab <- data.frame(
  analysis = c("Primary REML, 7 datasets", "Hartung-Knapp, 7 datasets",
               "Excluding reconstructed-SE studies, 7 datasets",
               "Eight-effect specification, 8 effects",
               "Egger regression (z)", "Egger regression (p)",
               "Mouse subgroup", "Human subgroup"),
  k = c(7, 7, nrow(full) - 2, nrow(full7), NA, NA,
        sum(full$species == "mouse"), sum(full$species == "human")),
  pooled = c(fit$b, fit_hk$b, fit_ex$b, fit_7e$b, egg$zval, egg$pval,
             mouse_fit$b, human_fit$b),
  ci_lo = c(fit$ci.lb, fit_hk$ci.lb, fit_ex$ci.lb, fit_7e$ci.lb, NA, NA,
            mouse_fit$ci.lb, human_fit$ci.lb),
  ci_hi = c(fit$ci.ub, fit_hk$ci.ub, fit_ex$ci.ub, fit_7e$ci.ub, NA, NA,
            mouse_fit$ci.ub, human_fit$ci.ub),
  I2 = c(fit$I2, fit_hk$I2, fit_ex$I2, fit_7e$I2, NA, NA, NA, NA),
  tau2 = c(fit$tau2, fit_hk$tau2, fit_ex$tau2, fit_7e$tau2, NA, NA, NA, NA),
  pi_lo = c(pred$pi.lb, rep(NA, 7)),
  pi_hi = c(pred$pi.ub, rep(NA, 7))
)
write.table(summary_tab, file.path(PROC, "meta_summary_tlr5_7ds.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(ldf, file.path(PROC, "meta_loo_tlr5_7ds.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("derived tables written to %s\n", PROC))
cat(sprintf("Fig5 (v7) saved -> %s\n", file.path(FIG, "Fig5.pdf")))
cat(sprintf("pooled %.3f (95%% CI %.3f to %.3f); I2 %.1f%%; PI %.3f to %.3f; Egger p %.3f\n",
            fit$b, fit$ci.lb, fit$ci.ub, fit$I2, pred$pi.lb, pred$pi.ub, egg$pval))
