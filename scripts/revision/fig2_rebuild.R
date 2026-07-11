# Figure 2 regeneration — unified set N=3,922 oncogenic positives (replaces old build 2,415)
# input : fig2_master_3922.tsv  (located alongside in paper_handoff/revision/)
# def   : AM-missed = AlphaMissense raw pathogenicity < 0.5  (1,071 / 3,922 = 27.3%)
# corner: phyloP100way < 4 & RSA > 0.4  ->  81% missed vs 15% elsewhere (S6 recomputed)
# usage : run from the same folder (or adjust paths)

library(ggplot2)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)  # for combined save (optional)

df <- read.delim("fig2_master_3922.tsv", stringsAsFactors = FALSE)
df$AM_status <- factor(df$AM_status, levels = c("AM caught", "AM missed"))
cols <- c("AM caught" = "#4C9F70", "AM missed" = "#C0504D")

N   <- nrow(df)                                  # 3922
nm  <- sum(df$AM_status == "AM missed")          # 1071
pct <- round(100 * nm / N, 1)                    # 27.3

# ---------------- Fig 2A : AM raw pathogenicity vs COSMIC recurrence ----------------
# y=recurrence count is highly skewed, so log10 scale. x=AM raw, vertical line at 0.5 (=missed cut).
# COSMIC recurrence=0 (positives labeled by hotspot/CIViC only, n=57) shown at log-axis floor (0.7).
n0 <- sum(df$cosmic_sample_count == 0)
df$rec_disp <- ifelse(df$cosmic_sample_count == 0, 0.7, df$cosmic_sample_count)
a <- ggplot(df, aes(am_pathogenicity, rec_disp, color = AM_status)) +
  geom_hline(yintercept = 0.85, linetype = "dotted", color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_point(size = 0.7, alpha = 0.35) +
  scale_color_manual(values = cols, name = NULL) +
  scale_y_log10(breaks = c(0.7, 1, 10, 100, 1000, 10000),
                labels = c("0", "1", "10", "100", "1000", "10000")) +
  labs(x = "AlphaMissense pathogenicity (raw)",
       y = "COSMIC recurrence (samples, log scale)",
       title = "A  Drivers underscored by AlphaMissense",
       subtitle = sprintf("n = %d positives; %d (%.1f%%) AM-missed (raw < 0.5)\n%d with 0 COSMIC count shown at axis floor",
                           N, nm, pct, n0)) +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.85, 0.85),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8))

# ---------------- Fig 2B : phyloP vs RSA + low-conservation/exposed corner ----------------
d2 <- df[!is.na(df$phyloP100way) & !is.na(df$struct_rsa), ]   # N = 3691
d2$corner <- d2$phyloP100way < 4 & d2$struct_rsa > 0.4
mi <- round(100 * mean(d2$AM_status[d2$corner]  == "AM missed"), 0)   # 81
mo <- round(100 * mean(d2$AM_status[!d2$corner] == "AM missed"), 0)   # 15

b <- ggplot(d2, aes(phyloP100way, struct_rsa, color = AM_status)) +
  geom_point(size = 0.7, alpha = 0.35) +
  # corner rectangle (phyloP<4 & RSA>0.4)
  annotate("rect", xmin = -Inf, xmax = 4, ymin = 0.4, ymax = Inf,
           fill = NA, color = "grey20", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = -2, y = 1.02, hjust = 0, vjust = 1, size = 3, color = "grey20",
           label = sprintf("low-conservation / exposed corner\n%d%% AM-missed vs %d%% elsewhere",
                           mi, mo)) +
  scale_color_manual(values = cols, name = NULL) +
  labs(x = "phyloP100way conservation", y = "Relative solvent accessibility (RSA)",
       title = "B  AM-missed drivers: low-conservation, exposed",
       subtitle = sprintf("n = %d positives with phyloP and RSA", nrow(d2))) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9))

cat(sprintf("Fig2A: N=%d, AM-missed=%d (%.1f%%)\n", N, nm, pct))
cat(sprintf("Fig2B corner(phyloP<4 & RSA>0.4): %d%% missed vs %d%% elsewhere (N=%d)\n",
            mi, mo, nrow(d2)))

# individual save
ggsave("fig2A_scatter_3922.pdf", a, width = 5.2, height = 4.5)
ggsave("fig2B_corner_3922.pdf",  b, width = 5.6, height = 4.5)
# combined (only if patchwork available; otherwise use the individual PDFs above)
if (has_patchwork) {
  combo <- patchwork::wrap_plots(a, b, widths = c(1, 1.1))
  ggsave("fig2_combined_3922.pdf", combo, width = 10.5, height = 4.5)
  ggsave("fig2_combined_3922.png", combo, width = 10.5, height = 4.5, dpi = 300)
} else {
  ggsave("fig2A_scatter_3922.png", a, width = 5.2, height = 4.5, dpi = 300)
  ggsave("fig2B_corner_3922.png",  b, width = 5.6, height = 4.5, dpi = 300)
  cat("patchwork not installed -> saved individual PDF/PNG (for combined, install.packages('patchwork') locally)\n")
}
