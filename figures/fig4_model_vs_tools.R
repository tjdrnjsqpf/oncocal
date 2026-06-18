# Fig 4a (main) — Model vs predictors: AUROC & AUPRC, mean±sd over 10 seeds, by role
# Input: tables/model_vs_tools.tsv (mean/sd 컬럼)  Output: fig4_model_vs_tools.{pdf,png}
# 메시지: 모델(스택) > best 단일 도구(robust); 세 모델 변형은 sd 내 동일(mech/struct 무효).
suppressPackageStartupMessages({library(tidyverse)})
PREV <- 0.188   # 양성 비율 = AUPRC 무작위 기준선

disp <- c("MODEL_tools"="OncoCal (tools)", "MODEL_tools+mech"="OncoCal (+mech)",
          "MODEL_tools+mech+struct"="OncoCal (+mech+struct)",
          "best_tool:VARITY_R"="VARITY-R", "AlphaMissense"="AlphaMissense")
pal <- c("OncoCal (tools)"="#C0392B","OncoCal (+mech)"="#D98880","OncoCal (+mech+struct)"="#E59866",
         "VARITY-R"="#1A5276","AlphaMissense"="#5DADE2")

d <- read.delim("tables/model_vs_tools.tsv", check.names = FALSE) %>%
  pivot_longer(-predictor, names_to = c("metric","stratum","stat"), names_sep = "_") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(stratum = recode(stratum, all="All", oncogene="Oncogene (GOF)", TSG="TSG (LOF)"),
         stratum = factor(stratum, levels=c("All","Oncogene (GOF)","TSG (LOF)")),
         metric = factor(metric, levels = c("AUROC","AUPRC")),   # AUROC 위, AUPRC 아래
         predictor = factor(recode(predictor, !!!disp), levels = names(pal)))

hl <- data.frame(metric = "AUPRC", yint = PREV)   # AUPRC facet에만 기준선
SHOW_LABELS <- TRUE    # 막대 위 숫자 표시(에러바 위로 배치). 끄려면 FALSE.

p <- ggplot(d, aes(predictor, mean, fill = predictor)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.25, linewidth = 0.4) +
  geom_hline(data = hl, aes(yintercept = yint), linetype = "dotted", color = "grey40") +
  {if (SHOW_LABELS) geom_text(aes(y = mean + sd, label = sprintf("%.3f", mean)),
                              vjust = -0.6, size = 2.3)} +
  facet_grid(metric ~ stratum) +
  scale_fill_manual(values = pal) +
  coord_cartesian(ylim = c(0.65, 1.0)) +
  labs(x = NULL, y = "Score (mean ± s.d., 10 seeds)",
       title = "Cancer-calibrated model vs predictors (gene-group OOF)",
       caption = "Dotted line in AUPRC panels = positive prevalence (0.188), the random-classifier baseline.") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none",
        plot.caption = element_text(size = 7.5, color = "grey40"))

ggsave("figures_R/fig4_model_vs_tools.pdf", p, width = 8.5, height = 5.5)
ggsave("figures_R/fig4_model_vs_tools.png", p, width = 8.5, height = 5.5, dpi = 300)
