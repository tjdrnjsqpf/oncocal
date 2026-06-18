# Fig 3 — Driver under-scoring rate per tool (unsupervised vs supervised)
# Input: tables/driver_miss_rate.tsv  Output: fig3_driver_miss_rate.{pdf,png}
# Message: unsupervised (unsup) tools miss high-recurrence drivers more often (14.3% vs 9.1%).
suppressPackageStartupMessages({library(tidyverse)})

# dbNSFP column name → official publication tool name (cleaning up underscores etc.). Keep original name if not in the mapping.
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

# Only one representative per family: exclude the duplicate variants below from the figure (summary stats use all 49)
drop_tools <- c("VARITY-ER","VARITY-R (LOO)","VARITY-ER (LOO)","PolyPhen-2 (HDIV)",
                "BayesDel (-AF)","Eigen-PC","GERP (92 mammals)",
                "phyloP (17-way primate)","phyloP (470-way mammalian)",
                "phastCons (17-way primate)","phastCons (470-way mammalian)",
                "SIFT4G","MisFit-S")

# type2 = principle-based 2-way classification (supervised: trained on clinical labels / unsupervised: conservation, evolution, language models)
d <- read.delim("tables/driver_miss_rate.tsv", check.names = FALSE) %>%
  rename(miss_rate = `miss_rate(<0.5)`) %>%
  mutate(grp = recode(type2, supervised = "Supervised", unsupervised = "Unsupervised"),
         tool = recode(tool, !!!disp)) %>%
  filter(!tool %in% drop_tools) %>%
  mutate(tool = fct_reorder(tool, miss_rate))

p <- ggplot(d, aes(miss_rate, tool, fill = grp)) +
  geom_col() +
  # Colors: use orange/grey here so as not to confuse with Fig 1A (red = under-scored)
  scale_fill_manual(values = c("Unsupervised" = "#E67E22",
                               "Supervised" = "#95A5A6"), name = NULL) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Fraction of high-recurrence drivers under-scored (rankscore < 0.5)",
       y = NULL, title = "Unsupervised tools miss activating drivers more often") +
  theme_classic(base_size = 10) + theme(legend.position = c(0.75, 0.15))

h <- max(5, nrow(d) * 0.17 + 1.2)   # adjust height to the number of remaining tools
ggsave("figures_R/fig3_driver_miss_rate.pdf", p, width = 6.5, height = h)
ggsave("figures_R/fig3_driver_miss_rate.png", p, width = 6.5, height = h, dpi = 300)
cat(sprintf("Number of tools shown: %d (duplicate variants excluded from all 49)\n", nrow(d)))
