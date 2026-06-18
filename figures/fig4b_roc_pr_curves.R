# Fig 4b — ROC & PR curves: model (tool stack) vs AlphaMissense vs VARITY_R
# Input: tables/oof_predictions.tsv  Output: fig4b_roc_pr_curves.{pdf,png}
# 메시지: 모델(스택) ~0.92 > best single tool ~0.87. full vs tools-only는 거의 동일(기전/구조 무효).
suppressPackageStartupMessages({library(tidyverse)})

d <- read.delim("tables/oof_predictions.tsv", check.names = FALSE)

# rank 기반 ROC/PR 계산 (sklearn 불필요)
roc_df <- function(y, s) {
  ok <- !is.na(s); y <- y[ok]; s <- s[ok]
  o <- order(s, decreasing = TRUE); y <- y[o]
  P <- sum(y == 1); N <- sum(y == 0)
  tibble(tpr = cumsum(y == 1)/P, fpr = cumsum(y == 0)/N)
}
pr_df <- function(y, s) {
  ok <- !is.na(s); y <- y[ok]; s <- s[ok]
  o <- order(s, decreasing = TRUE); y <- y[o]
  tp <- cumsum(y == 1); fp <- cumsum(y == 0); P <- sum(y == 1)
  tibble(recall = tp/P, precision = tp/(tp + fp))
}
auroc <- function(y, s){ ok<-!is.na(s); y<-y[ok]; s<-s[ok]
  r<-rank(s); P<-sum(y==1); N<-sum(y==1)*0+sum(y==0); (sum(r[y==1])-P*(P+1)/2)/(P*N) }

preds <- c(`OncoCal (tools)` = "p_model_tools",
           `OncoCal (+mech+struct)` = "p_model_full",
           AlphaMissense = "AlphaMissense",
           `VARITY-R` = "VARITY_R")

build <- function(fn) imap_dfr(preds, ~ fn(d$y_true, d[[.x]]) %>%
                                 mutate(predictor = .y,
                                        auc = sprintf("%.3f", auroc(d$y_true, d[[.x]]))))
roc <- build(roc_df); pr <- build(pr_df)

# 공통 팔레트 + AUROC 포함 라벨(두 패널 동일 → 범례 1개로 합쳐짐)
# 색은 Fig 3(막대)와 통일: AlphaMissense=연파랑, VARITY-R=진파랑
pal <- c("OncoCal (tools)" = "#C0392B", "OncoCal (+mech+struct)" = "#E59866",
         AlphaMissense = "#5DADE2", "VARITY-R" = "#1A5276")
auc_lab <- setNames(
  paste0(names(preds), " (AUROC ", sapply(preds, function(c) sprintf("%.3f", auroc(d$y_true, d[[c]]))), ")"),
  names(preds))
roc$predictor <- factor(roc$predictor, levels = names(preds))
pr$predictor  <- factor(pr$predictor,  levels = names(preds))

common_col <- scale_color_manual(values = pal, labels = auc_lab, name = NULL)
common_thm <- theme_classic(base_size = 11) +
  theme(legend.text = element_text(size = 8.5))

pROC <- ggplot(roc, aes(fpr, tpr, color = predictor)) +
  geom_abline(linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_line(linewidth = 0.9) + common_col + common_thm +
  coord_equal() + labs(x = "False positive rate", y = "True positive rate", title = "ROC") +
  theme(plot.margin = margin(6, 14, 6, 6))

pPR <- ggplot(pr, aes(recall, precision, color = predictor)) +
  geom_line(linewidth = 0.9) + common_col + common_thm +
  labs(x = "Recall", y = "Precision", title = "Precision-Recall") +
  theme(plot.margin = margin(6, 6, 6, 14))

# 두 패널 결합 — 범례를 아래로 모아 대각선 가림 없음
if (requireNamespace("patchwork", quietly = TRUE)) {
  p <- patchwork::wrap_plots(pROC, pPR, ncol = 2) +
       patchwork::plot_layout(guides = "collect") &
       theme(legend.position = "bottom")
} else { p <- pROC }
ggsave("figures_R/fig4b_roc_pr_curves.pdf", p, width = 11, height = 6)
ggsave("figures_R/fig4b_roc_pr_curves.png", p, width = 11, height = 6, dpi = 300)
