# ============================================================
# Fig3_panels_style.R - Fig3 ?????? PDF??????? + ?????
# ????data/processed + m8_robustness/data + data/processed/editorial
# ???classic ? strip??? 12 / ?? 14???????????? PDF
# ============================================================
PROC <- "C:/Users/yw_qi/Downloads/pa-host-ppi/data/processed"
M8   <- "C:/Users/yw_qi/Downloads/pa-host-ppi/m8_robustness/data"
ED   <- "C:/Users/yw_qi/Downloads/pa-host-ppi/data/processed/editorial"
FIG  <- "C:/Users/yw_qi/Downloads/pa-host-ppi/figures/style_review"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

read_tab <- function(f) {
  if (!file.exists(f)) stop("???????: ", f)
  read.delim(f, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(ggsignif); library(ggrepel)})
theme_clean <- theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold", colour = "black"),
        panel.border = element_blank(),
        axis.text = element_text(size = 12), axis.title = element_text(size = 13),
        legend.text = element_text(size = 11), legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10))
theme_set(theme_clean)

it <- read_tab(file.path(PROC, "independent_datasets_tlr5_check.tsv"))
prr <- read_tab(file.path(PROC, "mouse_prr_specificity.tsv"))
pd <- read_tab(file.path(PROC, "gse192528_per_donor_genes.tsv"))
jk <- read_tab(file.path(M8, "human_gse192528_leave1donor.tsv"))
dualF <- read_tab(file.path(PROC, "dual_evidence_targets.tsv"))
deseq <- read_tab(file.path(ED, "gse192528_deseq2_allgenes.tsv"))

# ---------- A ----------
six <- data.frame(ctx = c("Mouse lung\n(meta-DEG)","Mouse lung\n(SP-AKI)","Human lung\nexplants","Human corneal\nepithelium","Human iPSC\nmacrophages","Mouse MDSC\n(M-MDSC)"),
                  fc = c(-1.69, -0.90, -0.55, -2.08, -2.97, -1.41),
                  lab = c("FDR 0.014 *","p = 0.0009 *","ns (4/5 donors)","p = 0.0005 *","tier 3 (descriptive)","p = 0.25 ns"))
six$fc_abs <- abs(six$fc); six$ctx <- factor(six$ctx, levels = rev(six$ctx))
A3a <- ggplot(six, aes(fc, ctx)) +
  geom_col(aes(fill = fc_abs), width = 0.68) +
  scale_fill_gradient(low = "#F6E3CE", high = "#B03A2E", guide = "none") +
  geom_text(aes(x = -0.05, label = lab), hjust = 1, size = 3.8, colour = "#4E342E") +
  geom_vline(xintercept = 0) +
  labs(x = "TLR5 log2FC (PA vs control)", y = NULL,
       title = "A  TLR5 across six contexts",
       subtitle = "Per-context FDR, nominal p or Tier 3 annotation")
ggsave(file.path(FIG, "Fig3A.pdf"), A3a, width = 6.4, height = 4.6, device = "pdf", useDingbats = FALSE)

# ---------- B ----------
prr$logFC <- as.numeric(prr$logFC_ave); prr$padj <- as.numeric(prr$new_padj)
prr <- prr %>% filter(mouse_gene %in% c("Tlr5","Tlr2","Tlr4","Tlr6","Tlr9","Nod1","Nod2","Nlrc4","Nlrp3","Myd88","Tirap","Irak2","Nfkb1","Nfkb2","Il6","Il1b")) %>%
  arrange(logFC) %>% mutate(mouse_gene = factor(mouse_gene, levels = mouse_gene), star = ifelse(padj < 0.05, "*", ""))
A3b <- ggplot(prr, aes(logFC, mouse_gene)) +
  geom_col(fill = ifelse(prr$mouse_gene == "Tlr5", "#8E2A18",
                  ifelse(prr$logFC < 0, "#C65D3A", "#E89B6C")), width = 0.68) +
  geom_text(aes(label = star), hjust = ifelse(prr$logFC < 0, 1.15, -0.4), size = 4.2) +
  geom_vline(xintercept = 0) +
  labs(x = "logFC (mouse lung infection)", y = NULL,
       title = "B  Mouse receptor/adaptor response",
       subtitle = "TLR5 is the only significant negative gene shown; * FDR<0.05")
ggsave(file.path(FIG, "Fig3B.pdf"), A3b, width = 6.4, height = 4.6, device = "pdf", useDingbats = FALSE)

# ---------- C ----------
pd$log2expr <- as.numeric(pd$log2expr)
pd4 <- pd %>% filter(gene %in% c("TLR5","NFKB2","CXCL8","IL1B")) %>%
  mutate(condition = factor(condition, levels = c("Uninfected","PA")))
A3c <- ggplot(pd4, aes(condition, log2expr)) +
  geom_line(aes(group = donor), colour = "grey70") + geom_point(aes(colour = condition), size = 2.0) +
  scale_colour_manual(values = c("#5D4037","#B03A2E"), guide = "none") +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  geom_signif(comparisons = list(c("Uninfected","PA")), test = "t.test", map_signif_level = TRUE,
              test.args = list(paired = TRUE), textsize = 3.6, tip_length = 0.02) +
  labs(x = NULL, y = "log2(counts+1)", title = "C  Human lung explants: PA vs uninfected (paired, n=5)")
ggsave(file.path(FIG, "Fig3C.pdf"), A3c, width = 7.6, height = 3.8, device = "pdf", useDingbats = FALSE)

# ---------- D ----------
jk <- jk %>% filter(gene %in% c("TLR5","CAMP","MYD88","NFKB1","NFKB2","CXCL8","IL1B")) %>%
  mutate(across(c(ref_log2FC, jackknife_minFC, jackknife_maxFC), as.numeric)) %>%
  arrange(ref_log2FC) %>% mutate(gene = factor(gene, levels = gene))
A3d <- ggplot(jk, aes(ref_log2FC, gene)) +
  geom_segment(aes(x = jackknife_minFC, xend = jackknife_maxFC, y = gene, yend = gene),
               linewidth = 0.7, colour = "#C89B7B") +
  geom_point(aes(colour = ref_log2FC > 0), size = 3.2) +
  scale_colour_manual(values = c("#B03A2E","#E89B6C"), guide = "none") +
  geom_text(aes(x = ifelse(ref_log2FC < 0, jackknife_minFC, jackknife_maxFC),
                label = sprintf("%.2f", ref_log2FC)),
            hjust = ifelse(jk$ref_log2FC < 0, 1.2, -0.2), size = 3.8, colour = "#4E342E") +
  geom_vline(xintercept = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.18, 0.16))) +
  labs(x = "log2FC (PA vs uninfected; point = reference, whisker = 5-donor jackknife range)",
       y = NULL, title = "D  Descriptive paired-expression jackknife",
       subtitle = "Reference TLR5=-0.65; direction remains negative after each exclusion")
ggsave(file.path(FIG, "Fig3D.pdf"), A3d, width = 6.4, height = 4.6, device = "pdf", useDingbats = FALSE)

# ---------- E ----------
tc <- it %>% filter(dataset == "GSE276816_humaniPSCmacro" & gene == "TLR5") %>%
  select(cell_type, t60_vs_0_log2FC, t150_vs_0_log2FC, t450_vs_0_log2FC) %>%
  pivot_longer(-cell_type, names_to = "t", values_to = "fc") %>%
  mutate(h = as.numeric(sub("t(\\d+)_.*", "\\1", t)))
tc$cell_type <- factor(tc$cell_type, levels = c("WT_macrophage","CF_macrophage"), labels = c("CFTR WT","CFTR F508del"))
A3e <- ggplot(tc, aes(h, fc, colour = cell_type, group = cell_type)) +
  geom_line(linewidth = 1) + geom_point(size = 2.6) +
  scale_colour_manual(values = c("#B03A2E","#E89B6C")) +
  geom_hline(yintercept = 0) +
  labs(x = "minutes post PA infection", y = "TLR5 log2FC vs 0 min",
       title = "E  iPSC-macrophage TLR5 time course",
       subtitle = "Three measured time points; Tier 3 descriptive") +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(FIG, "Fig3E.pdf"), A3e, width = 6.4, height = 4.6, device = "pdf", useDingbats = FALSE)

# ---------- F ----------
dfF <- dualF %>% mutate(logFC_ave = as.numeric(logFC_ave)) %>%
  left_join(deseq %>% transmute(gene, log2FC_PA = as.numeric(lfc_PA_ctrl)),
            by = c("gene_clean" = "gene")) %>% filter(!is.na(log2FC_PA))
dfF$is_tlr5 <- dfF$gene_clean == "TLR5"
A3f <- ggplot(dfF, aes(logFC_ave, log2FC_PA)) +
  annotate("rect", xmin = -2.3, xmax = 0, ymin = -3.1, ymax = 0, fill = "#FBEAD8", alpha = 0.7) +
  geom_hline(yintercept = 0, colour = "#7f8c8d") + geom_vline(xintercept = 0, colour = "#7f8c8d") +
  geom_point(aes(size = is_tlr5), colour = ifelse(dfF$is_tlr5, "#B03A2E", "#C89B7B"), alpha = 0.9) +
  scale_size_manual(values = c("TRUE" = 3.6, "FALSE" = 2.6), guide = "none") +
  geom_text_repel(data = dfF %>% filter(!is_tlr5), aes(label = gene_clean), size = 3.6, colour = "#4E342E",
                  seed = 42, max.overlaps = 25, box.padding = 0.35, segment.colour = "#B89B8B") +
  geom_text_repel(data = dfF %>% filter(is_tlr5), aes(label = gene_clean), size = 4.2, colour = "#B03A2E",
                  fontface = "bold", seed = 42, box.padding = 0.4, segment.colour = "#B03A2E") +
  annotate("text", x = -2.25, y = -2.9, label = "convergent attenuation", size = 3.8, colour = "#8E2A18", hjust = 0, fontface = "italic") +
  annotate("text", x = 0.35, y = 4.9, label = "CAMP: species-divergent (up in mouse, down in human)",
           size = 3.6, colour = "#8E2A18", hjust = 0, fontface = "italic") +
  coord_cartesian(xlim = c(-2.3, 5.8), ylim = c(-3.1, 5.4)) +
  labs(x = "Mouse lung meta-DEG log2FC (PA infection)", y = "Human lung explant log2FC (PA vs uninfected)",
       title = "F  Mouse-human dual-evidence comparison",
       subtitle = "13/15 targets had donor-blocked human estimates; only TLR5 was negative in both")
ggsave(file.path(FIG, "Fig3F.pdf"), A3f, width = 6.4, height = 4.8, device = "pdf", useDingbats = FALSE)

# ---------- G ----------
tlr <- deseq[deseq$gene == "TLR5", ]
ag <- data.frame(contrast = c("PA vs control","PA vs IAV","PA vs BCG"),
                 lfc = c(tlr$lfc_PA_ctrl, tlr$lfc_PA_ctrl - tlr$lfc_IAV_ctrl, tlr$lfc_PA_ctrl - tlr$lfc_BCG_ctrl),
                 se = c(tlr$se_PA_ctrl, sqrt(tlr$se_PA_ctrl^2 + tlr$se_IAV_ctrl^2), sqrt(tlr$se_PA_ctrl^2 + tlr$se_BCG_ctrl^2)),
                 p = c(tlr$p_PA_ctrl, NA, NA))
ag$ci_lo <- ag$lfc - 1.96*ag$se; ag$ci_hi <- ag$lfc + 1.96*ag$se
ag$contrast <- factor(ag$contrast, levels = ag$contrast)
A3g <- ggplot(ag, aes(contrast, lfc)) +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, colour = "#C65D3A") +
  geom_point(size = 3.2, colour = "#B03A2E") +
  geom_text(aes(label = sprintf("lfc=%.2f", lfc)), vjust = 1.7, size = 3.6) +
  annotate("text", x = 1, y = 0.85, hjust = 0, size = 3.7,
           label = "TLR5 S-score percentile: 21.4% of expressed genes\n(39.2% among |PA lfc| >= 1)") +
  coord_cartesian(ylim = c(-1.8, 1.1)) +
  labs(x = NULL, y = "TLR5 log2FC (DESeq2, donor-blocked)",
       title = "G  Donor-blocked TLR5 contrasts",
       subtitle = "Human lung explants (GSE192528)")
ggsave(file.path(FIG, "Fig3G.pdf"), A3g, width = 6.4, height = 4.6, device = "pdf", useDingbats = FALSE)

cat("Fig3A-G saved ->", FIG, "\n")
