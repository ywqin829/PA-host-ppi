# Fig4 panels A-F (latest authoritative data, unified style)
PROC <- "C:/Users/yw_qi/Downloads/pa-host-ppi/data/processed"
ED   <- "C:/Users/yw_qi/Downloads/pa-host-ppi/data/processed/editorial"
FIG  <- "C:/Users/yw_qi/Downloads/pa-host-ppi/figures/style_review"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
read_tab <- function(f) { if (!file.exists(f)) stop("missing: ", f)
  read.delim(f, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) }
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})
theme_clean <- theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold", colour = "black"),
        panel.border = element_blank(),
        axis.text = element_text(size = 12), axis.title = element_text(size = 13),
        legend.text = element_text(size = 11), legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))
theme_set(theme_clean)

it <- read_tab(file.path(PROC, "independent_datasets_tlr5_check.tsv"))
obs <- read.delim(file.path(PROC, "gse274823_obs.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
sc <- read_tab(file.path(PROC, "gse274823_tlr5_celltype.tsv"))
equiv <- read_tab(file.path(ED, "gse287469_flic_equiv.tsv"))
logit <- read_tab(file.path(ED, "scRNA_tlr5_detection_logistic.tsv"))
deseq <- read_tab(file.path(ED, "gse192528_deseq2_allgenes.tsv"))
ttest <- read_tab(file.path(ED, "gse324794_ttest.tsv"))

# A
mdc <- it %>% filter(dataset == "GSE287469_mouseMDSC" & gene %in% c("Tlr5","Myd88")) %>%
  select(cell_type, gene, NT_mean, WT_mean, dFliC_mean, WT_vs_dFliC_p) %>%
  mutate(across(c(NT_mean, WT_mean, dFliC_mean, WT_vs_dFliC_p), as.numeric)) %>%
  pivot_longer(c(NT_mean, WT_mean, dFliC_mean), names_to = "grp", values_to = "expr") %>%
  mutate(grp = factor(grp, levels = c("NT_mean","WT_mean","dFliC_mean"), labels = c("NT","WT","dFliC")))
pvalM <- it %>% filter(dataset == "GSE287469_mouseMDSC" & gene %in% c("Tlr5","Myd88")) %>%
  mutate(across(c(WT_vs_dFliC_p), as.numeric))
ymM <- mdc %>% group_by(gene, cell_type) %>% summarise(ymax = max(expr), .groups = "drop")
brM <- pvalM %>% left_join(ymM, by = c("gene","cell_type")) %>%
  mutate(y = ymax * 1.18, lab = ifelse(WT_vs_dFliC_p < 0.05, "*", "")) %>% filter(lab != "")
A4a <- ggplot(mdc, aes(grp, expr, fill = grp)) +
  geom_col(width = 0.62) + facet_wrap(~ gene + cell_type, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("#D9C4A9","#B03A2E","#E89B6C")) +
  geom_segment(data = brM, aes(x = 1.92, xend = 3.08, y = y, yend = y), inherit.aes = FALSE) +
  geom_text(data = brM, aes(x = 2.5, y = y + 0.03 * ymax, label = lab), inherit.aes = FALSE, size = 4.2) +
  labs(x = NULL, y = "mean counts", title = "A  FliC exposure in mouse MDSC (bracket: WT vs dFliC; * p<0.05)") +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))
ggsave(file.path(FIG, "Fig4A.pdf"), A4a, width = 8.0, height = 3.8, device = "pdf", useDingbats = FALSE)

# B
comp <- obs %>% filter(genotype == "WT") %>% count(condition, celltype) %>%
  group_by(condition) %>% mutate(pct = n / sum(n) * 100)
ct_levels <- c("AlvMac","B","cDC","Epi","Mono","NK","Neu","Other","T","pDC")
comp$celltype <- factor(comp$celltype, levels = ct_levels)
comp$condition <- factor(comp$condition, levels = c("sham","PA"))
comp <- comp %>% arrange(condition, celltype) %>%
  group_by(condition) %>% mutate(cum = cumsum(pct) - 0.5 * pct)
pal4b <- c("#B03A2E","#C65D3A","#D98E73","#E8B48C","#F1D3A6","#C89B7B","#A67B5B","#8E6B4E","#6F543F","#5D4037")
names(pal4b) <- ct_levels
key_txt <- sprintf("sham:  AlvMac %.1f%%  |  cDC %.1f%%\nPA:  Neu %.1f%%  |  AlvMac %.1f%%  |  cDC %.1f%%",
                   comp$pct[comp$condition == "sham" & comp$celltype == "AlvMac"],
                   comp$pct[comp$condition == "sham" & comp$celltype == "cDC"],
                   comp$pct[comp$condition == "PA" & comp$celltype == "Neu"],
                   comp$pct[comp$condition == "PA" & comp$celltype == "AlvMac"],
                   comp$pct[comp$condition == "PA" & comp$celltype == "cDC"])
A4b <- ggplot(comp, aes(condition, pct, fill = celltype)) +
  geom_col(position = "stack", width = 0.62) +
  scale_fill_manual(values = pal4b) +
  labs(x = NULL, y = "% BALF cells (WT)", title = "B  PA remodels BALF composition (scRNA)",
       subtitle = key_txt) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.subtitle = element_text(size = 9, colour = "#3E2723", lineheight = 0.95))
ggsave(file.path(FIG, "Fig4B.pdf"), A4b, width = 6.4, height = 4.6, device = "pdf", useDingbats = FALSE)

# C
scp <- sc %>% filter(celltype %in% c("AlvMac","cDC") & genotype == "WT") %>%
  mutate(across(c(TLR5_mean_PA, TLR5_mean_sham), as.numeric)) %>%
  pivot_longer(c(TLR5_mean_PA, TLR5_mean_sham), names_to = "grp", values_to = "expr") %>%
  mutate(grp = factor(grp, levels = c("TLR5_mean_sham","TLR5_mean_PA"), labels = c("sham","PA")))
A4c <- ggplot(scp, aes(celltype, expr, fill = grp)) +
  geom_col(position = position_dodge(0.72), width = 0.55) +
  geom_text(aes(label = sprintf("%.3f", expr), group = grp),
            position = position_dodge(0.72), vjust = -0.4, size = 3.4) +
  scale_fill_manual(values = c("#D9C4A9","#B03A2E")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "TLR5 mean log-norm expr", title = "C  Lower TLR5 within AlvMac and cDC",
       subtitle = "Descriptive pooled libraries; n=1 per condition") +
  theme(legend.position = "top", legend.title = element_blank(), plot.subtitle = element_text(size = 9.5))
ggsave(file.path(FIG, "Fig4C.pdf"), A4c, width = 6.4, height = 4.2, device = "pdf", useDingbats = FALSE)

# D
inf <- logit[logit$term == "inf", ]
A4d <- ggplot(inf, aes("PA infection", OR)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey60") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.12, colour = "#C65D3A") +
  geom_point(size = 3.8, colour = "#B03A2E") +
  geom_text(aes(label = sprintf("OR = %.2f (95%% CI %.2f-%.2f), p=%.3g", OR, ci_lo, ci_hi, p)),
            vjust = -1.2, size = 3.5) +
  coord_cartesian(ylim = c(0, 1.08)) +
  labs(x = NULL, y = "Odds ratio for Tlr5 detection", title = "D  UMI-adjusted Tlr5 detection",
       subtitle = "PA vs sham; descriptive pooled-library model") +
  theme(plot.subtitle = element_text(size = 9.5))
ggsave(file.path(FIG, "Fig4D.pdf"), A4d, width = 5.8, height = 4.2, device = "pdf", useDingbats = FALSE)

# E
myd_fc <- c("Human lung explants" = deseq[deseq$gene == "MYD88", "lfc_PA_ctrl"][1],
            "Mouse lung (SP-AKI)" = ttest[tolower(ttest$gene) == "myd88", "lfc"][1],
            "Mouse lung (meta-DEG)" = 1.55, "MDSC M-" = 0.67, "MDSC PMN-" = -0.15)
myd <- data.frame(ctx = names(myd_fc), fc = as.numeric(myd_fc))
myd$ctx <- factor(myd$ctx, levels = myd$ctx[order(myd$fc)])
A4e <- ggplot(myd, aes(fc, ctx)) +
  geom_col(fill = ifelse(myd$fc < 0, "#B03A2E", "#E89B6C"), width = 0.62) + geom_vline(xintercept = 0) +
  labs(x = "MYD88 log2FC (PA vs control)", y = NULL, title = "E  MYD88 dynamics are cell-type-dependent (not co-suppressed)")
ggsave(file.path(FIG, "Fig4E.pdf"), A4e, width = 6.4, height = 4.4, device = "pdf", useDingbats = FALSE)

# F
equiv$label <- paste0(equiv$gene, " ", equiv$celltype)
equiv$effect_label <- sprintf("%.2f [%.1f, %.1f]; pTOST=%.2g", equiv$d, equiv$ci_lo, equiv$ci_hi, equiv$tost_p_0.5)
equiv$yn <- as.numeric(factor(equiv$label, levels = rev(unique(equiv$label))))
A4f <- ggplot(equiv, aes(d, yn)) +
  annotate("rect", xmin = -0.5, xmax = 0.5, ymin = 0.4, ymax = 4.6, fill = "#fef9e7", alpha = 0.7) +
  annotate("text", x = 0, y = 4.75, label = "equivalence bound +/-0.5", size = 3.6, colour = "#8A5A3B") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, y = yn, yend = yn),
               linewidth = 0.8, colour = ifelse(equiv$tost_p_0.5 < 0.05, "#2E7D32", "#B03A2E")) +
  geom_point(size = 3.2, colour = ifelse(equiv$tost_p_0.5 < 0.05, "#2E7D32", "#B03A2E")) +
  geom_text(aes(y = yn + 0.22, label = effect_label), hjust = 0.5, size = 2.8) +
  scale_y_continuous(breaks = equiv$yn, labels = equiv$label) +
  coord_cartesian(xlim = c(-6.2, 7.2), ylim = c(0.4, 4.6)) +
  labs(x = "WT vs dFliC log2FC contrast (limma-voom)", y = NULL,
       title = "F  WT-vs-dFliC equivalence analysis",
       subtitle = "TLR5 equivalence not established; TOST bound +/-0.5") +
  theme(plot.subtitle = element_text(size = 9.5))
ggsave(file.path(FIG, "Fig4F.pdf"), A4f, width = 6.4, height = 4.6, device = "pdf", useDingbats = FALSE)
cat("Fig4A-F saved ->", FIG, "\n")
