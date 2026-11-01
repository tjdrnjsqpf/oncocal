# Fig 2D (candidate) — conservation × RSA 2D: two orthogonal axes jointly predict AM failures
# Input: fig1_input.tsv (phyloP) + fig2_rsa_input.tsv (RSA), joined by pkey
# Output: fig2d_conservation_rsa_2d.{pdf,png}
# Resolve the repository root from this script's own location so the script runs
# correctly whether invoked as `Rscript figures/<name>.R` from the repo root or as
# `Rscript <name>.R` from inside figures/.
.a <- commandArgs(trailingOnly = FALSE)
.f <- sub("^--file=", "", .a[grep("^--file=", .a)])
if (length(.f)) setwd(normalizePath(file.path(dirname(.f[1]), "..")))
dir.create("figures_R", showWarnings = FALSE)   # output directory, not tracked in git

suppressPackageStartupMessages({library(tidyverse)})

a <- read.delim("tables/fig1_input.tsv", check.names = FALSE)[, c("pkey","AM_status","phyloP100way")]
b <- read.delim("tables/fig2_rsa_input.tsv", check.names = FALSE)[, c("pkey","struct_rsa")]
d <- inner_join(a, b, by = "pkey") %>%
  filter(AM_status %in% c("AM missed","AM caught"), !is.na(phyloP100way), !is.na(struct_rsa)) %>%
  mutate(AM_status = factor(AM_status, levels = c("AM missed","AM caught")))

p <- ggplot(d, aes(phyloP100way, struct_rsa, color = AM_status)) +
  # Shade the "low-conservation + surface-exposed" corner
  annotate("rect", xmin = -Inf, xmax = 4, ymin = 0.4, ymax = Inf,
           fill = "#C0392B", alpha = 0.06) +
  geom_point(alpha = 0.4, size = 1.3) +
  scale_color_manual(values = c("AM missed" = "#C0392B","AM caught" = "#1A5276"), name = NULL) +
  # (the 79% annotation is added manually by the user) — keep only the shaded corner as a visual guide
  labs(x = "phyloP100way conservation", y = "Relative solvent accessibility (RSA)",
       title = "Missed drivers occupy the low-conservation, exposed corner") +
  theme_classic(base_size = 12) + theme(legend.position = "top")

ggsave("figures_R/fig2d_conservation_rsa_2d.pdf", p, width = 6, height = 5.2)
ggsave("figures_R/fig2d_conservation_rsa_2d.png", p, width = 6, height = 5.2, dpi = 300)
