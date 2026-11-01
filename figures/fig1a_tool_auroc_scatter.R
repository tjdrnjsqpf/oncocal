# Fig 1a — Per-tool AUROC: oncogene (GOF) vs TSG (LOF)
# Input: tables/tool_performance.tsv   Output: fig1a_tool_auroc_scatter.{pdf,png}
# Message: 42 of 49 tools fall below the diagonal (GOF AUROC < TSG AUROC).
# Resolve the repository root from this script's own location so the script runs
# correctly whether invoked as `Rscript figures/<name>.R` from the repo root or as
# `Rscript <name>.R` from inside figures/.
.a <- commandArgs(trailingOnly = FALSE)
.f <- sub("^--file=", "", .a[grep("^--file=", .a)])
if (length(.f)) setwd(normalizePath(file.path(dirname(.f[1]), "..")))
dir.create("figures_R", showWarnings = FALSE)   # output directory, not tracked in git

suppressPackageStartupMessages({library(tidyverse)})

tp <- read.delim("tables/tool_performance.tsv", check.names = FALSE)
wide <- tp %>%
  filter(group %in% c("oncogene","TSG")) %>%
  select(tool, group, AUROC) %>%
  pivot_wider(names_from = group, values_from = AUROC) %>%
  mutate(below = oncogene < TSG)

n_below <- sum(wide$below, na.rm = TRUE); n_tot <- sum(!is.na(wide$below))

p <- ggplot(wide, aes(TSG, oncogene, color = below)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#2C3E50"),
                     labels = c(`TRUE` = "GOF under-scored", `FALSE` = "not under-scored"),
                     name = NULL) +
  coord_equal(xlim = c(0.5, 0.92), ylim = c(0.5, 0.92)) +
  labs(x = "AUROC, TSG (loss-of-function)",
       y = "AUROC, oncogene (gain-of-function)",
       title = sprintf("%d of %d tools score GOF drivers worse than LOF", n_below, n_tot)) +
  theme_classic(base_size = 12) + theme(legend.position = c(0.78, 0.12))

ggsave("figures_R/fig1a_tool_auroc_scatter.pdf", p, width = 5, height = 5)
ggsave("figures_R/fig1a_tool_auroc_scatter.png", p, width = 5, height = 5, dpi = 300)
