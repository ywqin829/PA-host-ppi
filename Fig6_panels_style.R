# Fig6 panels A-E (latest authoritative data, unified style; F omitted - user's own)
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

tr <- read_tab(file.path(PROC, "time_resolved_hubs_GSE262433.tsv"))
HVS <- read_tab(file.path(PROC, "hvs_v2_matrix.tsv")); HVS$score <- as.numeric(HVS$HVS_v2)
drug <- read_tab(file.path(PROC, "drug_direction_reasoned.tsv"))
rk <- read_tab(file.path(PROC, "tlr5_rank_under_systems.tsv"))
mcda <- read_tab(file.path(ED, "pts_mcda_summary.tsv"))
drop <- read_tab(file.path(ED, "pts_mcda_drop_spec.tsv"))

# A
seld <- c("TLR5","C3","VWF","RND1","MUC1","CFH","ITGB1","EZR")
dt <- tr %>% filter(gene %in% seld) %>%
  mutate(across(c(log2FC_24H, log2FC_120H, padj_120H), as.numeric)) %>%
  pivot_longer(c(log2FC_24H, log2FC_120H), names_to = "t", values_to = "fc") %>%
  mutate(t = factor(t, levels = c("log2FC_24H","log2FC_120H"), labels = c("24 h","120 h")),
         gene = factor(gene, levels = seld))
trp <- tr %>% filter(gene %in% seld) %>% mutate(across(c(padj_120H), as.numeric))
ymD <- dt %>% group_by(gene) %>% summarise(ymax = max(fc, na.rm = TRUE), ymin = min(fc, na.rm = TRUE), .groups = "drop") %>%
  mutate(span = pmax(ymax - ymin, 1))
brD <- trp %>% left_join(ymD, by = "gene") %>%
  mutate(y = pmax(ymax, 0) + 0.22 * span, lab = ifelse(padj_120H < 0.05, "*", "")) %>% filter(lab != "")
A6a <- ggplot(dt, aes(t, fc, fill = t)) +
  geom_col(position = "dodge", width = 0.62) + facet_wrap(~ gene, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("#E89B6C","#B03A2E")) +
  geom_segment(data = brD, aes(x = 0.92, xend = 2.08, y = y, yend = y), inherit.aes = FALSE) +
  geom_text(data = brD, aes(x = 1.5, y = y + 0.03 * span, label = lab), inherit.aes = FALSE, size = 4.2) +
  geom_hline(yintercept = 0) +
  labs(x = NULL, y = "signed log2FC (PAO1 vs control; >0 = induced)",
       title = "A  Persistent-infection boundary (Calu-3 ALI; * BH q<0.05)") +
  theme(legend.position = "none", plot.margin = margin(5.5, 5.5, 5.5, 30))
ggsave(file.path(FIG, "Fig6A.pdf"), A6a, width = 8.4, height = 4.2, device = "pdf", useDingbats = FALSE)

# B
A6b <- ggplot(HVS %>% arrange(desc(score)) %>% head(20) %>% mutate(gene = factor(gene, levels = rev(gene))),
              aes(score, gene)) +
  geom_col(aes(fill = score), width = 0.7) +
  scale_fill_gradient(low = "#F6E3CE", high = "#B03A2E", guide = "none") +
  labs(x = "HVS v2 score", y = NULL, title = "B  Host Vulnerability Score top-20") +
  theme(axis.text.y = element_text(size = 9.5))
ggsave(file.path(FIG, "Fig6B.pdf"), A6b, width = 6.4, height = 4.6, device = "pdf", useDingbats = FALSE)

# C
hsc <- HVS %>% select(gene, score)
dru <- drug %>% filter(gene %in% c("TLR5","ELANE","IFNG","CAMP","LTF")) %>% left_join(hsc, by = "gene") %>%
  mutate(dir2 = ifelse(grepl("^(AGONIZE|RESTORE)|ENHANCE", direction), "enhance", "inhibit"),
         gene = factor(gene, levels = gene[order(score)]))
A6c <- ggplot(dru, aes(score, gene, fill = dir2)) +
  geom_col(width = 0.62) + geom_text(aes(label = round(score, 1)), hjust = -0.2, size = 3.7) +
  scale_fill_manual(values = c("enhance" = "#2E7D32", "inhibit" = "#B03A2E")) +
  coord_cartesian(xlim = c(0, max(dru$score) * 1.25)) +
  labs(x = "Host Vulnerability Score (v2)", y = NULL, title = "C  Direction-reasoned candidates",
       subtitle = "green = enhance/restore; red = inhibit") +
  theme(legend.position = "none", plot.subtitle = element_text(size = 9.5))
ggsave(file.path(FIG, "Fig6C.pdf"), A6c, width = 6.4, height = 4.2, device = "pdf", useDingbats = FALSE)

# D
rk <- rk %>% mutate(across(c(HVS_v2, PTS), as.numeric)) %>% filter(!is.na(HVS_v2) & !is.na(PTS))
hl <- rk %>% filter(gene %in% c("TLR5","FN1","C3","VWF","ELANE","CAMP","MUC1"))
hl$lab <- ifelse(hl$gene %in% c("TLR5", "FN1"),
                 sprintf("%s (HVS %.1f, PTS %.1f)", hl$gene, hl$HVS_v2, hl$PTS), hl$gene)
A6d <- ggplot(rk, aes(HVS_v2, PTS)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "#7f8c8d") +
  geom_point(colour = "#D9C4A9", size = 1.8, alpha = 0.6) +
  geom_point(data = hl, colour = "#C65D3A", size = 2.8) +
  geom_point(data = hl %>% filter(gene == "TLR5"), colour = "#B03A2E", shape = 8, size = 4.5) +
  ggrepel::geom_text_repel(data = hl, aes(label = lab), size = 3.7, seed = 42, max.overlaps = 20,
                           colour = ifelse(hl$gene == "TLR5", "#B03A2E", "#4E342E"),
                           fontface = ifelse(hl$gene %in% c("TLR5","FN1"), "bold", "plain")) +
  coord_cartesian(xlim = c(0, 90), ylim = c(0, 90)) +
  labs(x = "HVS rank (low = more vulnerable)", y = "PTS rank (low = higher priority)",
       title = "D  Vulnerability is not therapeutic priority")
ggsave(file.path(FIG, "Fig6D.pdf"), A6d, width = 6.4, height = 5.2, device = "pdf", useDingbats = FALSE)

# E
sc <- data.frame(design = c(mcda$design, "drop PA preference"),
                 P_top1 = c(mcda$P_top1, drop$P_top1), P_top5 = c(mcda$P_top5, drop$P_top5),
                 median_rank = c(mcda$median_rank, drop$median_rank))
scm <- sc %>% pivot_longer(c(P_top1, P_top5), names_to = "metric", values_to = "prob")
scm$design <- factor(scm$design, levels = rev(unique(scm$design)))
A6e <- ggplot(scm, aes(prob, design, fill = metric)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  scale_fill_manual(values = c("P_top1" = "#B03A2E", "P_top5" = "#E89B6C"), name = NULL) +
  geom_text(aes(label = sprintf("%.0f%%", 100*prob)), position = position_dodge(0.7), hjust = -0.15, size = 3.6) +
  geom_text(data = sc, aes(x = 1.02, y = design, label = paste0("median rank ", median_rank)),
            inherit.aes = FALSE, hjust = 0, size = 3.6) +
  coord_cartesian(xlim = c(0, 1.55)) +
  labs(x = "Probability that TLR5 reaches the rank (top-1 requires no tie)", y = NULL,
       title = "E  PTS ranking is weight-sensitive",
       subtitle = "10,000 Dirichlet draws; dark = unique top-1; light = top-5") +
  theme(legend.position = "none")
ggsave(file.path(FIG, "Fig6E.pdf"), A6e, width = 6.4, height = 4.2, device = "pdf", useDingbats = FALSE)
cat("Fig6A-E saved ->", FIG, "\n")
