# Fig 1b (대안) — AlphaMissense raw 점수 vs COSMIC recurrence (oncogenic positives)
# Input: tables/fig1_input.tsv  Output: fig1b_recurrence_scatter.{pdf,png}
# 메시지: 대부분 고점수이지만, 고recurrence인데 0.5 미만으로 놓치는 임상적 소수(JAK2 V617F 등)가 존재.
suppressPackageStartupMessages({library(tidyverse)})
have_repel <- requireNamespace("ggrepel", quietly = TRUE)

d <- read.delim("tables/fig1_input.tsv", check.names = FALSE) %>%
  filter(role_class %in% c("oncogene","both"), !is.na(am_pathogenicity), cosmic_sample_count > 0) %>%
  mutate(status = ifelse(am_pathogenicity < 0.5, "Missed (score < 0.5)", "Caught (score ≥ 0.5)"))

# 라벨: 본문에 명시된 대표 '놓친' driver 고정(항상 표시) + 필요시 보강
named <- c("JAK2:p.V617F","FGFR3:p.Y373C","MPL:p.W515L","STAT3:p.Y640F")
lab <- d %>% filter(pkey %in% named)

p <- ggplot(d, aes(cosmic_sample_count, am_pathogenicity, color = status)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.55, size = 1.8) +
  scale_x_log10(labels = scales::comma,
                breaks = c(1,10,100,1000,10000,40000)) +
  scale_color_manual(values = c("Missed (score < 0.5)" = "#C0392B",
                                "Caught (score ≥ 0.5)" = "#1A5276"), name = NULL) +
  labs(x = "COSMIC recurrence (samples, log scale)",
       y = "AlphaMissense pathogenicity (raw)",
       title = "Recurrent activating drivers are not always scored pathogenic",
       subtitle = "Oncogenic positives; dashed line = pathogenic threshold (0.5)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")   # 범례를 플롯 밖(우측)으로

if (have_repel) {
  txt <- ggrepel::geom_text_repel(data = lab, aes(label = pkey), size = 3.2, color = "#7B241C",
                                  box.padding = 0.6, point.padding = 0.3,
                                  min.segment.length = 0, max.overlaps = Inf,
                                  segment.color = "grey50", seed = 1, show.legend = FALSE)
} else {
  txt <- geom_text(data = lab, aes(label = pkey), size = 3.2, color = "#7B241C",
                   hjust = -0.05, vjust = -0.4, show.legend = FALSE)
}
p <- p + txt

ggsave("figures_R/fig1b_recurrence_scatter.pdf", p, width = 7.5, height = 5)
ggsave("figures_R/fig1b_recurrence_scatter.png", p, width = 7.5, height = 5, dpi = 300)

# 요약 통계 (캡션/본문용)
miss <- d %>% summarise(n = n(), n_missed = sum(am_pathogenicity < 0.5),
                        pct = round(100*mean(am_pathogenicity < 0.5),1))
cat(sprintf("oncogenic positives n=%d, missed(<0.5)=%d (%.1f%%)\n", miss$n, miss$n_missed, miss$pct))
