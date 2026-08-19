# ============================================================
# Fig1BC_style.R - style preview for Fig1B / Fig1C (separate PDFs)
# ???????? processed ??????????????????
# ???classic ? strip??? 12 / ?? 14??????? PDF
# ============================================================
setwd("C:/Users/yw_qi/Downloads/pa-host-ppi/figure_reproduction")
PROC <- "C:/Users/yw_qi/Downloads/pa-host-ppi/data/processed"
FIG <- "C:/Users/yw_qi/Downloads/pa-host-ppi/figures/style_review"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

read_tab <- function(f) {
  if (!file.exists(f)) stop("???????: ", f)
  read.delim(f, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}
suppressPackageStartupMessages(library(ggplot2))
theme_clean <- theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold", colour = "black"),
        panel.border = element_blank(),
        axis.text = element_text(size = 12), axis.title = element_text(size = 13),
        plot.title = element_text(size = 14, face = "bold"))
theme_set(theme_clean)

# ---------- Fig1B: 132 curated PA-human pairs by source (from data) ----------
pos <- read_tab(file.path(PROC, "pa_human_positives.tsv"))
src_map <- c("frontiers_2025" = "Frontiers 2025",
             "frontiers_2025;literature" = "+ literature",
             "frontiers_2025;intact;literature" = "+ IntAct",
             "literature" = "literature only",
             "phi_base" = "PHI-base")
srcdf <- data.frame(source = names(src_map), n = as.integer(table(pos$source)[names(src_map)]))
srcdf$source <- factor(srcdf$source, levels = srcdf$source)
cat("Fig1B total:", sum(srcdf$n), "| breakdown:", paste(srcdf$source, srcdf$n, collapse = "; "), "\n")
b_src <- ggplot(srcdf, aes(source, n, fill = n)) +
  geom_col(width = 0.62) +
  scale_fill_gradient(low = "#F6E3CE", high = "#B03A2E", guide = "none") +
  geom_text(aes(label = n), vjust = -0.4, size = 3.8, colour = "#3E2723", fontface = "bold") +
  labs(x = NULL, y = "pairs", title = "B  132 curated PA-human pairs by source") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 12))
ggsave(file.path(FIG, "Fig1B.pdf"), b_src, width = 6.2, height = 4.6, device = "pdf", useDingbats = FALSE)

# ---------- Fig1C: 10,491 cross-species pairs by genus (from data) ----------
cs <- read_tab(file.path(PROC, "cross_species_training_positives.tsv"))
genus <- function(s) {
  s <- trimws(as.character(s))
  if (grepl("^E\\.? ?coli", s, ignore.case = TRUE) || grepl("^Escherichia", s, ignore.case = TRUE)) return("E. coli")
  strsplit(s, " ")[[1]][1]
}
cnt <- sort(table(vapply(cs$bacteria_species, genus, character(1))), decreasing = TRUE)
named <- names(cnt)[1:9]
other <- sum(cnt[!(names(cnt) %in% named)])
gendf <- data.frame(genus = c(named, "Other"), n = as.integer(c(cnt[named], other)))
gendf$genus <- factor(gendf$genus, levels = rev(gendf$genus))
cat("Fig1C total:", sum(gendf$n), "| breakdown:", paste(gendf$genus, gendf$n, collapse = "; "), "\n")
b_gen <- ggplot(gendf, aes(n, genus, fill = n)) +
  geom_col(width = 0.7) +
  scale_fill_gradient(low = "#F6E3CE", high = "#B03A2E", guide = "none") +
  geom_text(aes(label = n), hjust = -0.15, size = 3.8, colour = "#3E2723", fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "pairs", y = NULL, title = "C  10,491 cross-species pairs by genus (102 taxa)")
ggsave(file.path(FIG, "Fig1C_latest.pdf"), b_gen, width = 6.5, height = 4.6, device = "pdf", useDingbats = FALSE)
cat("Fig1B / Fig1C saved ->", FIG, "\n")
