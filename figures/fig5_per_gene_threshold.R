# Fig 5 — Per-gene oncogenic threshold heterogeneity (AlphaMissense)
# Input: tables/per_gene_calibration.tsv  Output: fig5_per_gene_threshold.{pdf,png}
# 메시지: 유전자별 oncogenic 임계가 0.07~0.99로 넓게 퍼짐(전역 0.469로 못 가림).
suppressPackageStartupMessages({library(tidyverse)})

d <- read.delim("tables/per_gene_calibration.tsv", check.names = FALSE) %>%
  filter(!is.na(am_threshold))
global <- 0.469

# (a) 분포 히스토그램
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

# (b) 극단 유전자 forest — 상/하위 10개씩(=20)으로 축소해 한 페이지에 맞춤
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

# (선택) 한 페이지용 결합본: a 위 / b 아래 (patchwork 있으면)
if (requireNamespace("patchwork", quietly = TRUE)) {
  comb <- patchwork::wrap_plots(pa, pb, ncol = 1, heights = c(1, 1.7)) +
          patchwork::plot_annotation(tag_levels = "a")
  ggsave("figures_R/fig5_per_gene_combined.pdf", comb, width = 5.8, height = 8.2)
  ggsave("figures_R/fig5_per_gene_combined.png", comb, width = 5.8, height = 8.2, dpi = 300)
}
