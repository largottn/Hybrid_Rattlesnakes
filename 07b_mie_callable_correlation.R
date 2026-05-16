# ============================================================================
# Step 7b - MIE count vs. callable sites: correlation analysis + scatter plot
# ----------------------------------------------------------------------------
# Reads a per-(sample, chromosome) summary CSV of callable site counts and
# MIE counts, computes Pearson and Spearman correlations and linear
# regressions, and produces a triangle-layout scatter plot (Fig. 3 in
# the accompanying paper).
#
# Inputs:   CSV with at minimum these columns:
#             chr             chromosome label (e.g. "chr1")
#             sample          sample display name (must match HYBRIDS below)
#             callable_sites  number of callable sites on this chromosome
#                             for this sample's trio
#             mie_count       number of MIEs on this chromosome for this
#                             sample's trio
# Outputs:  PNG, PDF, and TIFF of the scatter plot
#           Stats summary printed to stdout
#
# Required R packages: tidyverse, patchwork, scales
#
# Usage:    Set INPUT_FILE, OUT_DIR, and HYBRIDS below, then run with
#           Rscript or interactively.
# ============================================================================

# ---- User-configurable paths and parameters ----
INPUT_FILE <- "./MIE_callable_sites_per_chrom.csv"
OUT_DIR    <- "./figures"

# Display names of the samples to include, in legend order. These must
# match the values in the "sample" column of the input CSV.
HYBRIDS <- c("Hybrid 1", "Hybrid 2", "Hybrid 3")

# Desert palette: one color per hybrid; order matches HYBRIDS.
COLORS  <- c("#696965", "#A05A2C", "#B0A072")  # stone / rust / sand

# Position of the per-panel r-value annotation (in log10 axis units).
ANN_X <- 1.8     # Mb on the (log) x-axis
ANN_Y <- 80      # count on the (log) y-axis

# ---- Packages ----
suppressPackageStartupMessages({
  if (!require("tidyverse")) install.packages("tidyverse")
  if (!require("patchwork")) install.packages("patchwork")
  if (!require("scales"))    install.packages("scales")
  library(tidyverse)
  library(patchwork)
  library(scales)
})

# ---- Read data ----
data <- read_csv(INPUT_FILE, show_col_types = FALSE) %>%
  mutate(chr_label = gsub("chr", "", chr))

required_cols <- c("chr", "sample", "callable_sites", "mie_count")
missing       <- setdiff(required_cols, names(data))
if (length(missing) > 0) {
  stop("Input CSV is missing required column(s): ",
       paste(missing, collapse = ", "))
}

# ---- Pearson correlations ----
cat("\n=== PEARSON CORRELATIONS ===\n")
for (h in HYBRIDS) {
  sub <- data %>% filter(sample == h)
  ct  <- cor.test(sub$callable_sites, sub$mie_count, method = "pearson")
  cat(sprintf("%s: r = %.4f, R-squared = %.4f, p = %.2e\n",
              h, ct$estimate, ct$estimate^2, ct$p.value))
}

# ---- Spearman correlations ----
cat("\n=== SPEARMAN CORRELATIONS ===\n")
for (h in HYBRIDS) {
  sub <- data %>% filter(sample == h)
  ct  <- cor.test(sub$callable_sites, sub$mie_count, method = "spearman")
  cat(sprintf("%s: rho = %.4f, p = %.2e\n", h, ct$estimate, ct$p.value))
}

# ---- Linear regressions ----
cat("\n=== LINEAR REGRESSIONS (mie_count ~ callable_sites) ===\n")
for (h in HYBRIDS) {
  sub  <- data %>% filter(sample == h)
  fit  <- lm(mie_count ~ callable_sites, data = sub)
  s    <- summary(fit)
  slope     <- coef(fit)[["callable_sites"]]
  intercept <- coef(fit)[["(Intercept)"]]
  slope_ci  <- confint(fit)["callable_sites", ]
  int_ci    <- confint(fit)["(Intercept)", ]
  r2        <- s$r.squared
  fstat     <- s$fstatistic
  pval      <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)

  cat(sprintf(
    paste0(
      "\n%s:\n",
      "  slope     = %.3e (95%% CI: %.3e to %.3e)\n",
      "  intercept = %.2f (95%% CI: %.2f to %.2f)\n",
      "  R-squared = %.4f\n",
      "  p         = %.2e\n",
      "  ~%.3f MIEs per 1,000 callable sites\n"
    ),
    h, slope, slope_ci[1], slope_ci[2],
    intercept, int_ci[1], int_ci[2],
    r2, pval, slope * 1000
  ))
}

# ---- Per-chromosome MIE rates ----
cat("\n=== PER-CHROMOSOME MIE RATES ===\n")
data <- data %>%
  mutate(mie_rate = 100 * mie_count / callable_sites)
for (h in HYBRIDS) {
  sub <- data %>% filter(sample == h)
  cat(sprintf("\n%s (mean rate: %.3f%% +/- %.3f%%):\n",
              h, mean(sub$mie_rate), sd(sub$mie_rate)))
  for (i in 1:nrow(sub)) {
    cat(sprintf("  %5s: %.3f%%\n", sub$chr[i], sub$mie_rate[i]))
  }
}

# ---- Scatter plot panels ----
make_panel <- function(data, hybrid, col) {
  sub     <- data %>% filter(sample == hybrid)
  ct      <- cor.test(sub$callable_sites, sub$mie_count, method = "pearson")
  r_label <- sprintf("r = %.3f", ct$estimate)

  ggplot(sub, aes(x = callable_sites / 1e6, y = mie_count)) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.15,
                linewidth = 0.8, color = col, fill = col) +
    geom_point(size = 2.5, alpha = 0.85, color = col) +
    annotate("text", x = ANN_X, y = ANN_Y, label = r_label,
             hjust = 1, size = 4.5, fontface = "italic", color = col) +
    scale_x_log10(
      breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2),
      labels = c("50K", "100K", "200K", "500K", "1M", "2M")
    ) +
    scale_y_log10(
      breaks = c(100, 200, 500, 1000, 2000, 5000),
      labels = comma
    ) +
    labs(title = hybrid, x = "Callable Sites", y = "MIE Count") +
    theme_classic(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", color = "black",
                                      size = 13, hjust = 0.5),
      axis.text        = element_text(color = "black", size = 10),
      axis.title       = element_text(face = "bold", color = "black", size = 12),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.4),
      plot.margin      = margin(10, 15, 10, 15)
    )
}

panels <- mapply(make_panel,
                 hybrid = HYBRIDS,
                 col    = COLORS,
                 MoreArgs = list(data = data),
                 SIMPLIFY = FALSE)

# Triangle layout: two panels on top, one centered below
layout <- c(
  area(1, 1, 1, 2),
  area(1, 3, 1, 4),
  area(2, 2, 2, 3)
)
p_final <- wrap_plots(panels) + plot_layout(design = layout)
print(p_final)

# ---- Save figures ----
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_base <- file.path(OUT_DIR, "MIE_CallableSites_Correlation")

ggsave(paste0(out_base, ".png"),
       plot = p_final, width = 12, height = 10, dpi = 600, bg = "white")
ggsave(paste0(out_base, ".pdf"),
       plot = p_final, width = 12, height = 10, bg = "white")
ggsave(paste0(out_base, ".tif"),
       plot = p_final, width = 12, height = 10, dpi = 600, bg = "white",
       compression = "lzw")

cat("\nSaved scatter plot (PNG/PDF/TIF) to: ", OUT_DIR, "\n", sep = "")
