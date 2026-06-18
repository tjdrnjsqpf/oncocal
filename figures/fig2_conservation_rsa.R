# Fig 2 — Mechanism: missed drivers sit at low-conservation (a) and high-RSA (b) positions
# Input: tables/fig1_input.tsv (phyloP, GERP, AM_status)
#        + RSA: ⚠️ fig1_input.tsv 에 RSA 컬럼 없음 → 서버에서 struct_rsa export 필요(아래 참고).
# Output: fig2_conservation.{pdf,png} (+ fig2_rsa.* RSA TSV 도착 시)
suppressPackageStartupMessages({library(tidyverse)})

d <- read.delim("tables/fig1_input.tsv", check.names = FALSE) %>%
  filter(AM_status %in% c("AM caught","AM missed")) %>%
  mutate(AM_status = factor(AM_status, levels = c("AM missed","AM caught")))

box <- function(df, yvar, ylab, plabel) {
  ggplot(df, aes(AM_status, .data[[yvar]], fill = AM_status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.15, size = 0.6) +
    scale_fill_manual(values = c("AM missed" = "#C0392B","AM caught" = "#1A5276"), guide = "none") +
    # 상단 18% 여백 확보 후 p-value를 그 여백 안에 배치 → 박스와 안 겹침
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    annotate("text", x = 1.5, y = Inf, vjust = 1.6, label = plabel, size = 3.4, fontface = "italic") +
    labs(x = NULL, y = ylab) + theme_classic(base_size = 12)
}

pa <- box(d, "phyloP100way", "phyloP100way conservation", "Mann-Whitney p < 1e-200") +
  ggtitle("Missed oncogenic drivers sit at less-conserved positions")

ggsave("figures_R/fig2_conservation.pdf", pa, width = 4.5, height = 4.5)
ggsave("figures_R/fig2_conservation.png", pa, width = 4.5, height = 4.5, dpi = 300)

# (b) RSA — 서버에서 fig2_rsa_input.tsv (pkey, AM_status, struct_rsa) export 후 아래 활성화:
rsa_path <- "tables/fig2_rsa_input.tsv"
if (file.exists(rsa_path)) {
  r <- read.delim(rsa_path, check.names = FALSE) %>%
    mutate(AM_status = factor(AM_status, levels = c("AM missed","AM caught")))
  pb <- box(r, "struct_rsa", "Relative solvent accessibility (RSA)", "Mann-Whitney p = 2.6e-154") +
    ggtitle("Missed drivers are more solvent-exposed")
  ggsave("figures_R/fig2_rsa.pdf", pb, width = 4.5, height = 4.5)
  ggsave("figures_R/fig2_rsa.png", pb, width = 4.5, height = 4.5, dpi = 300)
} else {
  message("fig2_rsa: tables/fig2_rsa_input.tsv 없음 → 서버에서 struct_rsa export 필요.")
}
