# Fig 1b — Canonical driver x tool score heatmap, faceted by gene role (GOF vs LOF)
# Input: tables/canonical_driver_panel.tsv  Output: fig1b_canonical_heatmap.{pdf,png}
# Message: low scores (red) are scattered across the oncogene (GOF) block, contrasting with the TSG/both pattern.
suppressPackageStartupMessages({library(tidyverse)})

d0 <- read.delim("tables/canonical_driver_panel.tsv", check.names = FALSE)
tools <- c("AlphaMissense","REVEL","ESM1b","popEVE","MetaRNN")   # 0-1 rankscore family
TOP_ONCO <- 15   # for oncogenes, only the top N by recurrence (figure length limit; full set in Table 1)

# Select representatives per role: top N oncogenes + all both + all TSG
d <- bind_rows(
  d0 %>% filter(role_class == "oncogene") %>% slice_max(cosmic_sample_count, n = TOP_ONCO, with_ties = FALSE),
  d0 %>% filter(role_class != "oncogene")
)

role_lab <- c(oncogene = "Oncogene (GOF)", both = "Oncogene/TSG", TSG = "TSG (LOF)")
d <- d %>%
  mutate(role_f = factor(recode(role_class, !!!role_lab),
                         levels = c("Oncogene (GOF)","Oncogene/TSG","TSG (LOF)")),
         label = sprintf("%s (n=%s)", pkey, format(cosmic_sample_count, big.mark=","))) %>%
  arrange(role_f, cosmic_sample_count) %>%
  mutate(label = factor(label, levels = unique(label)))

long <- d %>% pivot_longer(all_of(tools), names_to = "tool", values_to = "score")

p <- ggplot(long, aes(tool, label, fill = score)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", score)), size = 2.4) +
  scale_fill_gradient2(low = "#C0392B", mid = "#F7F7F7", high = "#1A5276",
                       midpoint = 0.5, limits = c(0,1), name = "rankscore") +
  facet_grid(role_f ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL,
       title = "Canonical cancer drivers: predictor scores by gene role and recurrence",
       caption = sprintf("Top %d oncogenes by recurrence shown; full panel of %d drivers in Table 1.",
                         TOP_ONCO, nrow(d0))) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        strip.placement = "outside",
        plot.caption = element_text(size = 7.5, color = "grey40"))

h <- max(5.5, nrow(d) * 0.26 + 1.6)   # automatic height based on row count
ggsave("figures_R/fig1b_canonical_heatmap.pdf", p, width = 7.5, height = h, limitsize = FALSE)
ggsave("figures_R/fig1b_canonical_heatmap.png", p, width = 7.5, height = h, dpi = 300, limitsize = FALSE)
