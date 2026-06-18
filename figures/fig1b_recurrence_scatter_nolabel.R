# Fig 1b (no labels) — AM raw vs COSMIC recurrence, point labels removed (for manual annotation by the user)
# Input: tables/fig1_input.tsv  Output: fig1b_recurrence_scatter_nolabel.{pdf,png,svg}
suppressPackageStartupMessages({library(tidyverse)})

d <- read.delim("tables/fig1_input.tsv", check.names = FALSE) %>%
  filter(role_class %in% c("oncogene","both"), !is.na(am_pathogenicity), cosmic_sample_count > 0) %>%
  mutate(status = ifelse(am_pathogenicity < 0.5, "Missed (score < 0.5)", "Caught (score ≥ 0.5)"))

p <- ggplot(d, aes(cosmic_sample_count, am_pathogenicity, color = status)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.55, size = 1.8) +
  scale_x_log10(labels = scales::comma, breaks = c(1,10,100,1000,10000,40000)) +
  scale_color_manual(values = c("Missed (score < 0.5)" = "#C0392B",
                                "Caught (score ≥ 0.5)" = "#1A5276"), name = NULL) +
  labs(x = "COSMIC recurrence (samples, log scale)",
       y = "AlphaMissense pathogenicity (raw)",
       title = "Recurrent activating drivers are not always scored pathogenic",
       subtitle = "Oncogenic positives; dashed line = pathogenic threshold (0.5)") +
  theme_classic(base_size = 12) + theme(legend.position = "right")

ggsave("figures_R/fig1b_recurrence_scatter_nolabel.pdf", p, width = 7.5, height = 5)  # vector: editable in Illustrator/Inkscape
ggsave("figures_R/fig1b_recurrence_scatter_nolabel.png", p, width = 7.5, height = 5, dpi = 300)
# For SVG: install.packages("svglite") then uncomment the line below
# ggsave("figures_R/fig1b_recurrence_scatter_nolabel.svg", p, width = 7.5, height = 5)
