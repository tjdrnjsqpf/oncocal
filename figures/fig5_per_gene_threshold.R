# Fig 5 — Per-gene oncogenic threshold heterogeneity (AlphaMissense)
# Input: tables/per_gene_calibration.tsv  Output: fig5_per_gene_threshold.{pdf,png}
# Message: per-gene oncogenic thresholds span widely from 0.07 to 0.99 (a global 0.469 cannot cover this).
# Resolve the repository root from this script's own location so the script runs
# correctly whether invoked as `Rscript figures/<name>.R` from the repo root or as
# `Rscript <name>.R` from inside figures/.
.a <- commandArgs(trailingOnly = FALSE)
.f <- sub("^--file=", "", .a[grep("^--file=", .a)])
if (length(.f)) setwd(normalizePath(file.path(dirname(.f[1]), "..")))
dir.create("figures_R", showWarnings = FALSE)   # output directory, not tracked in git

suppressPackageStartupMessages({library(tidyverse)})

d <- read.delim("tables/per_gene_calibration.tsv", check.names = FALSE) %>%
  filter(!is.na(am_threshold))
global <- 0.469

# (a) Distribution histogram
pa <- ggplot(d, aes(am_threshold)) +
  geom_histogram(bins = 30, fill = "#1A5276", alpha = 0.85) +
  geom_vline(xintercept = global, linetype = "dashed", color = "#C0392B") +
  annotate("text", x = global + 0.02, y = Inf, vjust = 2, hjust = 0,
           label = "global Youden = 0.469", color = "#C0392B", size = 3.2) +
  labs(x = "Per-gene AlphaMissense oncogenic threshold", y = "Genes",
       title = sprintf("Per-gene thresholds span %.2f-%.2f (n=%d genes)",
                       min(d$am_threshold), max(d$am_threshold), nrow(d))) +
  theme_classic(base_size = 11)

ggsave("figures_R/fig5_per_gene_threshold.pdf", pa, width = 5.5, height = 3.3)
ggsave("figures_R/fig5_per_gene_threshold.png", pa, width = 5.5, height = 3.3, dpi = 300)

# (b) Extreme-gene forest — reduced to top/bottom 10 each (=20) to fit on one page
N_EACH <- 10
ext <- d %>% arrange(am_threshold) %>%
  slice(c(1:N_EACH, (n()-N_EACH+1):n())) %>% mutate(gene = fct_inorder(gene))
pal_role <- c(oncogene = "#C0392B", TSG = "#1A5276", both = "#8E44AD",
              unknown = "#95A5A6", fusion_only = "#27AE60")
pb <- ggplot(ext, aes(am_threshold, gene, color = role_class)) +
  geom_segment(aes(x = global, xend = am_threshold, yend = gene), color = "grey80") +
  geom_point(size = 2.4) +
  geom_vline(xintercept = global, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = pal_role, name = "Role") +
  labs(x = "Per-gene AlphaMissense oncogenic threshold", y = NULL,
       title = "Genes requiring very low vs very high thresholds") +
  theme_classic(base_size = 10) + theme(legend.position = "bottom")
ggsave("figures_R/fig5b_per_gene_forest.pdf", pb, width = 5.5, height = 5.0)
ggsave("figures_R/fig5b_per_gene_forest.png", pb, width = 5.5, height = 5.0, dpi = 300)

# (optional) Single-page combined version: a top / b bottom (if patchwork is available)
if (requireNamespace("patchwork", quietly = TRUE)) {
  comb <- patchwork::wrap_plots(pa, pb, ncol = 1, heights = c(1, 1.7)) +
          patchwork::plot_annotation(tag_levels = "a")
  ggsave("figures_R/fig5_per_gene_combined.pdf", comb, width = 5.8, height = 8.2)
  ggsave("figures_R/fig5_per_gene_combined.png", comb, width = 5.8, height = 8.2, dpi = 300)
}
