# Fig 2D (후보) — 보존도 × RSA 2D: 두 직교 축이 함께 AM 실패를 예측
# Input: fig1_input.tsv (phyloP) + fig2_rsa_input.tsv (RSA), pkey로 join
# Output: fig2d_conservation_rsa_2d.{pdf,png}
suppressPackageStartupMessages({library(tidyverse)})

a <- read.delim("tables/fig1_input.tsv", check.names = FALSE)[, c("pkey","AM_status","phyloP100way")]
b <- read.delim("tables/fig2_rsa_input.tsv", check.names = FALSE)[, c("pkey","struct_rsa")]
d <- inner_join(a, b, by = "pkey") %>%
  filter(AM_status %in% c("AM missed","AM caught"), !is.na(phyloP100way), !is.na(struct_rsa)) %>%
  mutate(AM_status = factor(AM_status, levels = c("AM missed","AM caught")))

p <- ggplot(d, aes(phyloP100way, struct_rsa, color = AM_status)) +
  # "저보존+표면노출" 코너 음영
  annotate("rect", xmin = -Inf, xmax = 4, ymin = 0.4, ymax = Inf,
           fill = "#C0392B", alpha = 0.06) +
  geom_point(alpha = 0.4, size = 1.3) +
  scale_color_manual(values = c("AM missed" = "#C0392B","AM caught" = "#1A5276"), name = NULL) +
  # (79% 주석은 사용자가 직접 추가) — 음영 코너만 시각 가이드로 유지
  labs(x = "phyloP100way conservation", y = "Relative solvent accessibility (RSA)",
       title = "Missed drivers occupy the low-conservation, exposed corner") +
  theme_classic(base_size = 12) + theme(legend.position = "top")

ggsave("figures_R/fig2d_conservation_rsa_2d.pdf", p, width = 6, height = 5.2)
ggsave("figures_R/fig2d_conservation_rsa_2d.png", p, width = 6, height = 5.2, dpi = 300)
