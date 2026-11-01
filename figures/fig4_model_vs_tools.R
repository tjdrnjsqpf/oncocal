# Fig 4a (main) — Model vs predictors: AUROC & AUPRC, mean±sd over 10 seeds, by role
# Input: tables/model_vs_tools.tsv (mean/sd columns)  Output: fig4_model_vs_tools.{pdf,png}
# Message: model (stack) > best single tool (robust); the three model variants are identical within sd (mech/struct add nothing).
# Resolve the repository root from this script's own location so the script runs
# correctly whether invoked as `Rscript figures/<name>.R` from the repo root or as
# `Rscript <name>.R` from inside figures/.
.a <- commandArgs(trailingOnly = FALSE)
.f <- sub("^--file=", "", .a[grep("^--file=", .a)])
if (length(.f)) setwd(normalizePath(file.path(dirname(.f[1]), "..")))
dir.create("figures_R", showWarnings = FALSE)   # output directory, not tracked in git

suppressPackageStartupMessages({library(tidyverse)})
PREV <- 0.188   # positive prevalence = AUPRC random baseline

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
         metric = factor(metric, levels = c("AUROC","AUPRC")),   # AUROC top, AUPRC bottom
         predictor = factor(recode(predictor, !!!disp), levels = names(pal)))

hl <- data.frame(metric = "AUPRC", yint = PREV)   # baseline only on the AUPRC facet
SHOW_LABELS <- TRUE    # show numbers above the bars (placed above the error bars). Set FALSE to disable.

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
