# Fig 4b — ROC & PR curves: model (tool stack) vs AlphaMissense vs VARITY_R
# Input: tables/oof_predictions.tsv  Output: fig4b_roc_pr_curves.{pdf,png}
# Message: model (stack) ~0.92 > best single tool ~0.87. full vs tools-only nearly identical (mechanism/structure add nothing).
suppressPackageStartupMessages({library(tidyverse)})

d <- read.delim("tables/oof_predictions.tsv", check.names = FALSE)

# rank-based ROC/PR computation (no sklearn needed)
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

# Shared palette + labels with AUROC (identical across both panels → merged into one legend)
# Colors consistent with Fig 3 (bars): AlphaMissense=light blue, VARITY-R=dark blue
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

# Combine the two panels — collect the legend at the bottom so the diagonal is not obscured
if (requireNamespace("patchwork", quietly = TRUE)) {
  p <- patchwork::wrap_plots(pROC, pPR, ncol = 2) +
       patchwork::plot_layout(guides = "collect") &
       theme(legend.position = "bottom")
} else { p <- pROC }
ggsave("figures_R/fig4b_roc_pr_curves.pdf", p, width = 11, height = 6)
ggsave("figures_R/fig4b_roc_pr_curves.png", p, width = 11, height = 6, dpi = 300)
