# ============================================================================
# Step 7a - Circular Manhattan plot of MIE density (350 kb windows)
# ----------------------------------------------------------------------------
# Reads per-trio MIE position TSVs (one file per offspring, output of the
# MIE pipeline step 4) and plots Mendelian inheritance error density as
# three concentric tracks on a circular ideogram of the autosomal genome.
#
# Inputs:   one TSV per offspring named <PREFIX><ID>.mie.positions.tsv,
#           with two unnamed columns: chromosome accession, position
# Outputs:  PNG, TIFF, and PDF of the circular plot
#
# Required R packages: tidyverse, circlize
#
# Usage:    Set INPUT_DIR, OUT_DIR, SAMPLE_MAP, and WINDOW_SIZE below.
#           If using a different reference genome, edit the chromosome
#           accession-to-name map in the "Chromosome map" section.
# ============================================================================

# ---- User-configurable paths and parameters ----
INPUT_DIR   <- "./mie_positions"
OUT_DIR     <- "./figures"
WINDOW_SIZE <- 350000   # bp; window size for binning MIE counts

# Mapping from sample ID (as it appears in the filename, stripped of the
# fixed prefix below) to the display name used in the plot legend.
SAMPLE_MAP <- c(
  "6001" = "Hybrid 1",
  "6002" = "Hybrid 2",
  "6004" = "Hybrid 3"
)
FILE_PREFIX  <- "animal"   # filename prefix stripped before lookup in SAMPLE_MAP
FILE_PATTERN <- paste0("^", FILE_PREFIX, ".*\\.mie\\.positions\\.tsv$")

# Desert palette: one color per sample. Order must match SAMPLE_MAP.
TRACK_COLORS <- c("#696965", "#A05A2C", "#B0A072")  # stone / rust / sand

# Quantile cap for the y-axis (suppresses extreme spikes that would
# otherwise compress the rest of the signal to a flat line).
Y_CAP_QUANTILE <- 0.99

# ---- Packages ----
suppressPackageStartupMessages({
  if (!require("tidyverse")) install.packages("tidyverse")
  if (!require("circlize"))  install.packages("circlize")
  library(tidyverse)
  library(circlize)
})

# ---- Chromosome map ----
# Maps GenBank/RefSeq sequence accessions to short chromosome labels for
# display. Edit this table to match your reference assembly.
# (CroVir_3.0: 17 autosomes, chr8 absent from the assembly.)
chrom_map <- tribble(
  ~accession,   ~new_name,
  "CM012306.1", "chr1",
  "CM012307.1", "chr2",
  "CM012308.1", "chr3",
  "CM012309.1", "chr4",
  "CM012310.1", "chr5",
  "CM012311.1", "chr6",
  "CM012312.1", "chr7",
  "CM012313.1", "chr9",
  "CM012314.1", "chr10",
  "CM012315.1", "chr11",
  "CM012316.1", "chr12",
  "CM012317.1", "chr13",
  "CM012318.1", "chr14",
  "CM012319.1", "chr15",
  "CM012320.1", "chr16",
  "CM012321.1", "chr17",
  "CM012322.1", "chr18"
)

# ---- Track-color helpers ----
sample_names <- unname(SAMPLE_MAP)
fills <- setNames(adjustcolor(TRACK_COLORS, alpha.f = 0.75), sample_names)
bgs   <- setNames(adjustcolor(TRACK_COLORS, alpha.f = 0.12), sample_names)
cols  <- setNames(TRACK_COLORS,                              sample_names)

# ---- Read MIE position files ----
files <- list.files(INPUT_DIR, pattern = FILE_PATTERN, full.names = TRUE)
if (length(files) == 0) {
  stop("No MIE position files found in ", INPUT_DIR,
       " matching pattern ", FILE_PATTERN)
}

# Strip prefix + suffix to recover the sample ID, then map to display name
raw_ids       <- gsub("\\.mie\\.positions\\.tsv$", "", basename(files))
raw_ids       <- gsub(paste0("^", FILE_PREFIX), "", raw_ids)
display_names <- SAMPLE_MAP[raw_ids]

if (any(is.na(display_names))) {
  stop("File ID(s) not in SAMPLE_MAP: ",
       paste(raw_ids[is.na(display_names)], collapse = ", "))
}

read_errors <- function(file, name) {
  read_tsv(file, col_names = c("chr", "pos"), show_col_types = FALSE) %>%
    mutate(sample = name)
}
data_raw <- map2_dfr(files, display_names, read_errors)

# ---- Clean and bin ----
chr_lengths_bp <- data_raw %>%
  inner_join(chrom_map, by = c("chr" = "accession")) %>%
  group_by(new_name) %>%
  summarise(length = max(pos) + WINDOW_SIZE, .groups = "drop") %>%
  rename(chr = new_name)

data_clean <- data_raw %>%
  inner_join(chrom_map, by = c("chr" = "accession")) %>%
  select(-chr) %>%
  rename(chr = new_name)

desired_order <- chrom_map$new_name
data_clean$chr <- factor(data_clean$chr, levels = desired_order)

# Bin per (sample, chr, window)
data_binned <- data_clean %>%
  mutate(window_start = floor(pos / WINDOW_SIZE) * WINDOW_SIZE,
         window_end   = window_start + WINDOW_SIZE) %>%
  group_by(sample, chr, window_start, window_end) %>%
  summarise(count = n(), .groups = "drop")

# Fill in zero-count windows so the tracks are continuous
all_windows <- chr_lengths_bp %>%
  rowwise() %>%
  do(data.frame(chr          = .$chr,
                window_start = seq(0, .$length, by = WINDOW_SIZE))) %>%
  ungroup() %>%
  mutate(window_end = window_start + WINDOW_SIZE) %>%
  crossing(sample = sample_names)

data_binned <- all_windows %>%
  left_join(data_binned, by = c("sample", "chr", "window_start", "window_end")) %>%
  mutate(count = replace_na(count, 0))
data_binned$chr <- factor(data_binned$chr, levels = desired_order)

# ---- Prepare circlize input ----
genome_df <- chr_lengths_bp %>%
  mutate(start = 0) %>%
  select(chr, start, end = length) %>%
  filter(chr %in% desired_order)
genome_df$chr <- factor(genome_df$chr, levels = desired_order)
genome_df     <- genome_df %>% arrange(chr)

# Cap counts at the chosen quantile to keep visual scale informative
y_max <- quantile(data_binned$count, Y_CAP_QUANTILE)

tracks <- lapply(sample_names, function(s) {
  data_binned %>%
    filter(sample == s) %>%
    arrange(chr, window_start) %>%
    mutate(count = pmin(count, y_max))
})
names(tracks) <- sample_names

# ---- Plot ----
draw_plot <- function() {
  circos.clear()
  circos.par(
    "track.height" = 0.22,
    "gap.degree"   = 1,
    "cell.padding" = c(0, 0, 0, 0),
    "start.degree" = 90
  )
  circos.genomicInitialize(
    genome_df,
    plotType = NULL,
    major.by = 50000000,
    tickLabelsStartFromZero = TRUE
  )

  # Chromosome labels
  circos.track(
    ylim = c(0, 1), track.height = 0.05, bg.border = NA,
    panel.fun = function(x, y) {
      label <- gsub("chr", "", CELL_META$sector.index)
      circos.text(CELL_META$xcenter, 0.5, label,
                  cex = 0.8, facing = "bending.inside",
                  niceFacing = TRUE, font = 2)
    }
  )

  # One filled-area track per sample
  for (s in sample_names) {
    circos.genomicTrack(
      tracks[[s]] %>% select(chr, window_start, window_end, count),
      ylim = c(0, y_max),
      track.height = 0.22,
      bg.col    = bgs[s],
      bg.border = cols[s],
      bg.lwd    = 0.5,
      panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, type = "area",
                            col = fills[s], border = cols[s],
                            lwd = 0.6, area = TRUE, baseline = 0)
      }
    )
  }

  legend("center",
         legend    = sample_names,
         fill      = fills[sample_names],
         border    = cols[sample_names],
         cex       = 0.9, bty = "n",
         text.font = 2, text.col = "black")

  circos.clear()
}

# ---- Save figures ----
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_base <- file.path(OUT_DIR, "MIE_Circular_Manhattan_350kb")

png(paste0(out_base, ".png"),
    width = 10, height = 10, units = "in", res = 600)
draw_plot(); dev.off()

tiff(paste0(out_base, ".tif"),
     width = 10, height = 10, units = "in", res = 600, compression = "lzw")
draw_plot(); dev.off()

pdf(paste0(out_base, ".pdf"), width = 10, height = 10)
draw_plot(); dev.off()

cat("Saved circular Manhattan plot (PNG/TIF/PDF) to: ", OUT_DIR, "\n", sep = "")
