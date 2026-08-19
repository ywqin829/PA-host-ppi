# ============================================================
# Fig2_panels_style.R - Fig2 ?????? PDF????? + ?????
# ????data/processed + m8_robustness/data???????
# ???classic ? strip??? 12 / ?? 14???????????? PDF
# ============================================================
PROC <- "C:/Users/yw_qi/Downloads/pa-host-ppi/data/processed"
M8   <- "C:/Users/yw_qi/Downloads/pa-host-ppi/m8_robustness/data"
FIG  <- "C:/Users/yw_qi/Downloads/pa-host-ppi/figures/style_review"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

read_tab <- function(f) {
  if (!file.exists(f)) stop("???????: ", f)
  read.delim(f, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})
theme_clean <- theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold", colour = "black"),
        panel.border = element_blank(),
        axis.text = element_text(size = 12), axis.title = element_text(size = 13),
        legend.text = element_text(size = 11), legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10))
theme_set(theme_clean)

hub <- read_tab(file.path(PROC, "hub_convergence_pa_enrichment.tsv"))
hub$n_sp <- as.numeric(hub$n_bacterial_species); hub$is_pa <- hub$is_pa == "True"
hub <- hub[!is.na(hub$n_sp) & hub$n_sp >= 1, ]
hub$decile <- as.integer(hub$decile)
dec <- hub %>% group_by(decile) %>% summarise(frac_pa = mean(is_pa) * 100, .groups = "drop")
glob <- mean(hub$is_pa) * 100
gmap <- read_tab(file.path(PROC, "human_uniprot_gene_map.tsv"))
hubv <- read_tab(file.path(PROC, "hub_intxdb_validation.tsv")) %>%
  left_join(gmap, by = c("human_uniprot" = "uniprot")) %>%
  mutate(gene = ifelse(is.na(gene), human_gene, gene)) %>% distinct(human_uniprot, .keep_all = TRUE)
selh <- c("FN1","CDC42","RAC1","RHOA","PLG","TP53","CFH","C3","MYD88")
loo <- read_tab(file.path(PROC, "m8_leave_one_pathogen_out.tsv"))
deg <- read_tab(file.path(M8, "hub_degree_stratified_auc.tsv"))
tp  <- read_tab(file.path(PROC, "gse192528_three_pathogen_pvalues.tsv"))

# ---------- A ----------
ymA <- max(dec$frac_pa)
A2a <- ggplot(dec, aes(decile, frac_pa)) +
  geom_col(aes(fill = decile), width = 0.7) +
  scale_fill_gradient(low = "#F6E3CE", high = "#B03A2E", guide = "none") +
  geom_hline(yintercept = glob, linetype = 2, colour = "#8E6B4E") +
  annotate("text", x = 1.0, y = glob + 0.30, label = "expected (all hubs)", size = 4.0, colour = "#8E6B4E", hjust = 0) +
  annotate("segment", x = 9, xend = 10, y = ymA * 1.18, yend = ymA * 1.18, colour = "#B03A2E", linewidth = 0.8) +
  annotate("segment", x = 9, xend = 9, y = ymA * 1.12, yend = ymA * 1.18, colour = "#B03A2E", linewidth = 0.8) +
  annotate("segment", x = 10, xend = 10, y = ymA * 1.12, yend = ymA * 1.18, colour = "#B03A2E", linewidth = 0.8) +
  annotate("text", x = 9.45, y = ymA * 1.32, label = "top-2 deciles\nOR 3.47, p=3.8e-6",
           size = 3.3, colour = "#B03A2E", fontface = "bold", lineheight = 0.9) +
  coord_cartesian(ylim = c(0, ymA * 1.7)) +
  labs(x = "Cross-species-degree decile (10 = highest)", y = "% PA targets",
       title = "A  PA-target enrichment by degree decile",
       subtitle = "Top two deciles; stratified-permutation p=1.0e-4")
ggsave(file.path(FIG, "Fig2A.pdf"), A2a, width = 6.4, height = 4.6, device = "pdf", useDingbats = FALSE)

# ---------- B ----------
dh <- hubv %>% filter(gene %in% selh) %>%
  mutate(intact = as.numeric(n_bacterial_species), itx = as.numeric(intxdb_species)) %>%
  arrange(desc(intact)) %>% mutate(gene = factor(gene, levels = gene)) %>%
  tidyr::pivot_longer(c(intact, itx), names_to = "src", values_to = "n")
dh$src <- factor(dh$src, levels = c("intact","itx"), labels = c("IntAct","InTxDB"))
dv <- dh %>% tidyr::pivot_wider(names_from = src, values_from = n) %>% tidyr::drop_na()
rhop <- suppressWarnings(cor.test(dv$IntAct, dv$InTxDB, method = "spearman"))
A2b <- ggplot(dh, aes(n, gene, fill = src)) +
  geom_col(position = position_dodge(0.72), width = 0.6) +
  scale_fill_manual(values = c("#B03A2E","#E89B6C")) +
  coord_cartesian(xlim = c(0, max(dh$n) * 1.35)) +
  labs(x = "co-targeting species", y = NULL, title = "B  Orthogonal InTxDB comparison",
       subtitle = paste0("Shared hubs; Spearman rho=", round(rhop$estimate, 2),
                         ", p=", formatC(rhop$p.value, format = "e", digits = 1))) +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(FIG, "Fig2B.pdf"), A2b, width = 6.4, height = 4.6, device = "pdf", useDingbats = FALSE)

# ---------- C ----------
loo$auc <- as.numeric(loo$auc_roc)
loo2 <- loo %>% filter(experiment != "neg_resample_5seeds") %>%
  mutate(exp2 = sub("leave_", "leave ", experiment), exp2 = ifelse(exp2 == "full", "full", exp2))
loo2$exp2 <- factor(loo2$exp2, levels = c("full","leave Yersinia","leave Bacillus","leave Francisella","leave Salmonella","leave Chlamydia","leave Escherichia"))
A2c <- ggplot(loo2, aes(exp2, auc)) +
  geom_col(fill = ifelse(loo2$exp2 == "full", "#B03A2E", "#E8C9A8"), width = 0.68) +
  geom_hline(yintercept = loo2$auc[loo2$exp2 == "full"], linetype = 2, colour = "#B03A2E") +
  geom_text(aes(label = sprintf("%.3f", auc)), vjust = -0.4, size = 4.2) +
  coord_cartesian(ylim = c(0.68, 0.80)) +
  labs(x = NULL, y = "transfer AUC-ROC", title = "C  Leave-one-pathogen-out transfer",
       subtitle = "AUC range 0.738-0.755") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 10))
ggsave(file.path(FIG, "Fig2C.pdf"), A2c, width = 6.4, height = 4.6, device = "pdf", useDingbats = FALSE)

# ---------- D ----------
deg$auc <- as.numeric(deg$auc); deg$lo <- as.numeric(deg$ci95_low); deg$hi <- as.numeric(deg$ci95_high)
deg$med <- as.numeric(deg$median_nspp)
ddeg <- deg[deg$stratum == "random (all humans)" | grepl("deg>=", deg$stratum), ]
t300 <- deg[deg$stratum == "hub: top-300", ]
A2d <- ggplot(ddeg, aes(med, auc)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "#E89B6C") +
  geom_line(colour = "#C65D3A") + geom_point(size = 2.6, colour = "#C65D3A") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, colour = "#C65D3A") +
  geom_point(data = t300, aes(med, auc), shape = 4, size = 3.4, colour = "#B03A2E") +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "#666") +
  coord_cartesian(ylim = c(0.15, 0.95)) +
  scale_x_continuous(breaks = sort(unique(ddeg$med)), labels = c("0","1","2","3","5","8.5")) +
  labs(x = "Negative-pool median host degree", y = "transfer AUC-ROC",
       title = "D  Hard-negative transfer by host degree",
       subtitle = "Spearman rho=-0.98, p=3.4e-5")
ggsave(file.path(FIG, "Fig2D.pdf"), A2d, width = 6.4, height = 4.6, device = "pdf", useDingbats = FALSE)

# ---------- E ----------
nr <- loo %>% filter(experiment == "neg_resample_5seeds")
vals <- as.numeric(strsplit(gsub("[^0-9.,-]", " ", gsub(".*roc5.: ", "", nr$note)), ",")[[1]])
nrd <- data.frame(seed = seq_along(vals), auc = vals)
A2e <- ggplot(nrd, aes(seed, auc)) +
  geom_point(size = 3, colour = "#D35400") +
  geom_hline(yintercept = mean(vals), linetype = 2, colour = "#D35400") +
  geom_hline(yintercept = 0.7405, linetype = 3, colour = "#8E6B4E") +
  coord_cartesian(ylim = c(0.60, 0.80)) +
  labs(x = "negative-set resampling (seed)", y = "transfer AUC-ROC",
       title = "E  Negative-set resampling",
       subtitle = paste0("mean ", round(mean(vals), 3), "; range 0.63-0.77; full model 0.740"))
ggsave(file.path(FIG, "Fig2E.pdf"), A2e, width = 6.4, height = 4.6, device = "pdf", useDingbats = FALSE)

# ---------- F ----------
tp$log2FC <- as.numeric(tp$log2FC); tp$p <- as.numeric(tp$p)
sel3 <- c("COL1A1","TF","CAMP","FN1","ELANE","MUC1","NLRC4","TLR5")
d3 <- tp %>% filter(gene %in% sel3) %>%
  mutate(FC = 2^log2FC,
         pathogen = factor(pathogen, levels = c("PA","IAV","BCG")),
         star = ifelse(p < 0.05, "*", ""), gene = factor(gene, levels = sel3))
uninf <- data.frame(gene = factor(rep(sel3, each = 1), levels = sel3),
                    pathogen = "Uninfected", log2FC = 0, FC = 1, p = 1, star = "")
d4 <- rbind(d3, uninf)
d4$pathogen <- factor(d4$pathogen, levels = c("Uninfected","PA","IAV","BCG"))
d4s <- d4 %>% filter(star != "") %>% mutate(xn = as.numeric(pathogen), ystar = FC * 1.18 + 0.06)
A2f <- ggplot(d4, aes(pathogen, FC, fill = pathogen)) +
  geom_col(width = 0.68) +
  geom_hline(yintercept = 1, colour = "grey45", linetype = 2) +
  geom_text(data = d4s, aes(x = xn, y = ystar, label = "*"), inherit.aes = FALSE,
            size = 5.5, colour = "black", fontface = "bold") +
  scale_fill_manual(values = c("Uninfected" = "#D9C4A9", "PA" = "#B03A2E", "IAV" = "#D9B48F", "BCG" = "#E89B6C")) +
  scale_x_discrete(labels = c("Uninfected" = "Ctrl", "PA" = "PA", "IAV" = "IAV", "BCG" = "BCG")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
  facet_wrap(~ gene, nrow = 2, scales = "free") +
  labs(x = NULL, y = "Fold change vs uninfected",
       title = "F  Human-lung three-pathogen response",
       subtitle = "n=5 donors; FC=2^log2FC; * nominal paired p<0.05") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(size = 7.5))
ggsave(file.path(FIG, "Fig2F.pdf"), A2f, width = 7.5, height = 6.2, device = "pdf", useDingbats = FALSE)

cat("Fig2A-F saved ->", FIG, "\n")
