# Fig 5 — External validation on independent DMS (ProteinGym v1.3)
# Input: tables/external_validation_dms.tsv  Output: fig6_dms_spearman.{pdf,png}
# Message: model |Spearman| 0.492 ~= AlphaMissense 0.502 > REVEL 0.426; direction correct (generalises).
# Resolve the repository root from this script's own location so the script runs
# correctly whether invoked as `Rscript figures/<name>.R` from the repo root or as
# `Rscript <name>.R` from inside figures/.
.a <- commandArgs(trailingOnly = FALSE)
.f <- sub("^--file=", "", .a[grep("^--file=", .a)])
if (length(.f)) setwd(normalizePath(file.path(dirname(.f[1]), "..")))
dir.create("figures_R", showWarnings = FALSE)   # output directory, not tracked in git

suppressPackageStartupMessages({library(tidyverse)})

# ProteinGym assay ID -> clean label "GENE (Author Year)"
clean_label <- function(assay) {
  g  <- sub("_HUMAN.*", "", assay)
  yr <- str_extract(assay, "(19|20)[0-9]{2}")
  au <- str_match(assay, "_HUMAN_([A-Za-z]+)")[, 2]
  lab <- paste0(g, " (", au, " ", yr, ")")
  ifelse(grepl("_CD_", assay), sub("\\)$", ", CD)", lab), lab)  # disambiguate the duplicate SRC Ahler 2019
}

d <- read.delim("tables/external_validation_dms.tsv", check.names = FALSE) %>%
  filter(!is.na(spearman_MODEL) & spearman_MODEL != "") %>%
  mutate(across(starts_with("spearman_"), as.numeric),
         label = clean_label(assay)) %>%
  pivot_longer(starts_with("spearman_"), names_to = "predictor",
               names_prefix = "spearman_", values_to = "rho") %>%
  mutate(absrho = abs(rho),
         predictor = factor(recode(predictor, MODEL = "OncoCal"),
                            levels = c("OncoCal","AlphaMissense","REVEL")),
         label = fct_reorder(label, absrho, .fun = max))

# Colors: consistent with Fig 3 (OncoCal red, AlphaMissense light blue), REVEL grey
pal <- c("OncoCal" = "#C0392B", "AlphaMissense" = "#5DADE2", "REVEL" = "#7F8C8D")

p <- ggplot(d, aes(absrho, label, fill = predictor)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
  labs(x = "|Spearman| with DMS fitness", y = NULL,
       title = "Generalisation to independent deep mutational scanning") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.justification = "left",
        panel.grid.major.x = element_line(color = "grey92"))

ggsave("figures_R/fig6_dms_spearman.pdf", p, width = 6.5, height = 5.5)
ggsave("figures_R/fig6_dms_spearman.png", p, width = 6.5, height = 5.5, dpi = 300)
