# Fig 3 — Driver under-scoring rate per tool (unsupervised vs supervised)
# Input: tables/driver_miss_rate.tsv  Output: fig3_driver_miss_rate.{pdf,png}
# 메시지: 비지도(unsup) 도구가 고recurrence driver를 더 많이 놓침(14.3% vs 9.1%).
suppressPackageStartupMessages({library(tidyverse)})

# dbNSFP 컬럼명 → 출판용 정식 도구명 (언더바 등 정리). 매핑에 없으면 원명 유지.
disp <- c(
  "phyloP17way_primate"="phyloP (17-way primate)", "phastCons17way_primate"="phastCons (17-way primate)",
  "phyloP100way_vertebrate"="phyloP (100-way vertebrate)", "phastCons100way_vertebrate"="phastCons (100-way vertebrate)",
  "phyloP470way_mammalian"="phyloP (470-way mammalian)", "phastCons470way_mammalian"="phastCons (470-way mammalian)",
  "GERPpp_RS"="GERP++", "GERP_92_mammals"="GERP (92 mammals)",
  "Polyphen2_HDIV"="PolyPhen-2 (HDIV)", "Polyphen2_HVAR"="PolyPhen-2 (HVAR)",
  "Eigen-PC-raw_coding"="Eigen-PC", "Eigen-raw_coding"="Eigen", "CADD_raw"="CADD",
  "fathmm-XF_coding"="FATHMM-XF", "BayesDel_addAF"="BayesDel (+AF)", "BayesDel_noAF"="BayesDel (-AF)",
  "MisFit_S"="MisFit-S", "MisFit_D"="MisFit-D",
  "VARITY_R"="VARITY-R", "VARITY_ER"="VARITY-ER", "VARITY_R_LOO"="VARITY-R (LOO)", "VARITY_ER_LOO"="VARITY-ER (LOO)")

# 계열당 대표 하나만: 아래 중복 변형은 그림에서 제외 (요약 통계는 전체 49개 기준)
drop_tools <- c("VARITY-ER","VARITY-R (LOO)","VARITY-ER (LOO)","PolyPhen-2 (HDIV)",
                "BayesDel (-AF)","Eigen-PC","GERP (92 mammals)",
                "phyloP (17-way primate)","phyloP (470-way mammalian)",
                "phastCons (17-way primate)","phastCons (470-way mammalian)",
                "SIFT4G","MisFit-S")

# type2 = 원칙 기반 2-way 분류(supervised: 임상라벨 학습 / unsupervised: 보존도·진화·언어모델)
d <- read.delim("tables/driver_miss_rate.tsv", check.names = FALSE) %>%
  rename(miss_rate = `miss_rate(<0.5)`) %>%
  mutate(grp = recode(type2, supervised = "Supervised", unsupervised = "Unsupervised"),
         tool = recode(tool, !!!disp)) %>%
  filter(!tool %in% drop_tools) %>%
  mutate(tool = fct_reorder(tool, miss_rate))

p <- ggplot(d, aes(miss_rate, tool, fill = grp)) +
  geom_col() +
  # 색: Fig 1A(빨강=under-scored)와 헷갈리지 않게 1B는 주황/회색으로
  scale_fill_manual(values = c("Unsupervised" = "#E67E22",
                               "Supervised" = "#95A5A6"), name = NULL) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Fraction of high-recurrence drivers under-scored (rankscore < 0.5)",
       y = NULL, title = "Unsupervised tools miss activating drivers more often") +
  theme_classic(base_size = 10) + theme(legend.position = c(0.75, 0.15))

h <- max(5, nrow(d) * 0.17 + 1.2)   # 남은 도구 수에 맞춰 높이 조정
ggsave("figures_R/fig3_driver_miss_rate.pdf", p, width = 6.5, height = h)
ggsave("figures_R/fig3_driver_miss_rate.png", p, width = 6.5, height = h, dpi = 300)
cat(sprintf("표시 도구 수: %d (전체 49개 중 중복 변형 제외)\n", nrow(d)))
