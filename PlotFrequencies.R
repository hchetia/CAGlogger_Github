#!/usr/bin/env Rscript

# =============================================================================
# plot_repeat_frequencies.R
#
# Generates publication-quality plots from *_FrequencyOfRepeats.csv files.
# Produces per-sample full-range and zoomed plots, plus combined dashboards.
#
# Usage:
#   Rscript plot_repeat_frequencies.R [input_dir] [output_dir]
#
#   input_dir  : directory containing *_FrequencyOfRepeats.csv files (default: .)
#   output_dir : directory for output plots (default: ./plots)
#
# Dependencies: ggplot2, dplyr, tidyr, scales, patchwork, RColorBrewer
# =============================================================================

# --- User-configurable: inherited allele positions --------------------------
# Set to NULL to disable, or provide a named list per sample:
#   inherited_alleles <- list(
#     "SampleA" = c(17, 45, 117),
#     "SampleB" = c(20, 50)
#   )
# Or set a single vector to apply the same lines to ALL samples:
#   inherited_alleles <- c(17, 45, 117)
#inherited_alleles <- NULL
inherited_alleles <- c(17,45)
# Style for allele annotation lines
allele_line_color  <- "red3"
allele_line_type   <- "dashed"
allele_line_width  <- 0.8
allele_label_size  <- 2.8

# --- Install missing packages quietly --------------------------------------
required_pkgs <- c("ggplot2", "dplyr", "tidyr", "scales", "patchwork",
                   "RColorBrewer", "stringr")

for (pkg in required_pkgs) {
 if (!requireNamespace(pkg, quietly = TRUE)) {
  install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
 }
}

suppressPackageStartupMessages({
 library(ggplot2)
 library(dplyr)
 library(tidyr)
 library(scales)
 library(patchwork)
 library(RColorBrewer)
 library(stringr)
})

# --- Parse command-line arguments -------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_dir  <- ifelse(length(args) >= 1, args[1], ".")
output_dir <- ifelse(length(args) >= 2, args[2], file.path(input_dir, "plots"))

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Discover input files ---------------------------------------------------
csv_files <- list.files(input_dir, pattern = "_FrequencyOfRepeats\\.csv$",
                        full.names = TRUE)

if (length(csv_files) == 0) {
 stop("No *_FrequencyOfRepeats.csv files found in: ", input_dir)
}

cat(sprintf("Found %d file(s) in '%s'\n", length(csv_files), input_dir))

# --- Read and combine all files ---------------------------------------------
all_data <- bind_rows(lapply(csv_files, function(f) {
 df <- read.csv(f, stringsAsFactors = FALSE)
 sample_name <- basename(f) %>%
  str_remove("_FrequencyOfRepeats\\.csv$")
 df$Sample <- sample_name
 return(df)
}))

all_data <- all_data %>%
 mutate(Repeat_Length = as.integer(Repeat_Length),
        Frequency    = as.numeric(Frequency))

# --- Colour palette (warm teal-coral gradient) ------------------------------
palette_pub <- c("#1B7F79", "#E07B54", "#5B5EA6", "#D4A03C", "#8C4B7A",
                 "#3D8EB9", "#C75D4A", "#6BAF6E", "#D17BC1", "#4ECDC4")

n_samples   <- length(unique(all_data$Sample))
sample_cols  <- if (n_samples == 1) {
 setNames("#1B7F79", unique(all_data$Sample))
} else {
 cols <- rep_len(palette_pub, n_samples)
 setNames(cols, sort(unique(all_data$Sample)))
}

# --- Publication theme ------------------------------------------------------
theme_pub <- function(base_size = 12) {
 theme_minimal(base_size = base_size) %+replace%
  theme(
   text              = element_text(color = "grey20"),
   axis.title        = element_text(face = "bold", size = rel(1.05)),
   axis.text         = element_text(color = "grey30", size = rel(0.9)),
   axis.line         = element_line(color = "grey50", linewidth = 0.4),
   axis.ticks        = element_line(color = "grey50", linewidth = 0.3),
   plot.title        = element_text(face = "bold", size = rel(1.3),
                                    hjust = 0.5, margin = margin(b = 2)),
   plot.subtitle     = element_text(size = rel(0.85), hjust = 0.5,
                                    color = "grey45", margin = margin(b = 10)),
   legend.title      = element_text(face = "bold", size = rel(0.9)),
   legend.text       = element_text(size = rel(0.8)),
   legend.position   = "bottom",
   legend.background = element_rect(fill = "transparent", color = NA),
   panel.grid.major.y = element_line(color = "grey88", linewidth = 0.25),
   panel.grid.major.x = element_blank(),
   panel.grid.minor   = element_blank(),
   plot.margin        = margin(14, 14, 10, 10),
   strip.background   = element_rect(fill = "grey96", color = NA),
   strip.text         = element_text(face = "bold", size = rel(0.95))
  )
}

# --- Helper: get allele positions for a sample ------------------------------
get_alleles <- function(sample_name) {
 if (is.null(inherited_alleles)) return(NULL)
 if (is.numeric(inherited_alleles)) return(inherited_alleles)
 if (is.list(inherited_alleles) && sample_name %in% names(inherited_alleles)) {
  return(inherited_alleles[[sample_name]])
 }
 return(NULL)
}

# --- Helper: add allele vlines to a ggplot ----------------------------------
add_allele_lines <- function(p, alleles, ymax = NULL) {
 if (is.null(alleles) || length(alleles) == 0) return(p)
 allele_df <- data.frame(x = alleles)
 p <- p +
  geom_vline(data = allele_df, aes(xintercept = x),
             linetype = allele_line_type,
             color    = allele_line_color,
             linewidth = allele_line_width) +
  geom_label(data = allele_df,
             aes(x = x, y = Inf, label = x),
             vjust    = 1.3,
             size     = allele_label_size,
             color    = allele_line_color,
             fill     = "white",
             fontface = "bold",
             label.size = 0.2,
             label.padding = unit(0.15, "lines"))
 return(p)
}

# --- Helper: detect dominant peaks ------------------------------------------
find_peaks <- function(freq, window = 3, min_prominence = 0.05) {
 n <- length(freq)
 if (n < 2 * window + 1) return(integer(0))
 peaks <- integer(0)
 max_freq <- max(freq, na.rm = TRUE)
 for (i in (window + 1):(n - window)) {
  neighbourhood <- freq[max(1, i - window):min(n, i + window)]
  if (freq[i] == max(neighbourhood) &&
      freq[i] > min_prominence * max_freq) {
   peaks <- c(peaks, i)
  }
 }
 return(peaks)
}

# --- Define zoom ranges -----------------------------------------------------
zoom_ranges <- list(
 "10to120" = c(10, 120),
 "35to120" = c(35, 120)
)

# =============================================================================
# PLOT 1: Per-sample bar plots – full range (linear)
# =============================================================================
cat("Generating per-sample full-range bar plots...\n")

for (samp in unique(all_data$Sample)) {
 df_s <- filter(all_data, Sample == samp)
 col  <- sample_cols[samp]
 alleles <- get_alleles(samp)
 
 peak_idx <- find_peaks(df_s$Frequency)
 df_peaks <- df_s[peak_idx, ] %>%
  arrange(desc(Frequency)) %>%
  slice_head(n = 5)
 
 p <- ggplot(df_s, aes(x = Repeat_Length, y = Frequency)) +
  geom_col(fill = col, color = "black", width = 0.85) +
  geom_text(data = df_peaks,
            aes(label = Repeat_Length),
            vjust = -0.5, size = 3, fontface = "bold", color = "grey30") +
  scale_x_continuous(breaks = pretty_breaks(n = 10),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(title    = samp,
       subtitle = "Distribution of tandem repeat lengths",
       x = "Repeat Length (bp)", y = "Frequency") +
  theme_pub()
 
 p <- add_allele_lines(p, alleles)
 
 fname <- file.path(output_dir, paste0(samp, "_barplot_full.pdf"))
 ggsave(fname, p, width = 10, height = 5, dpi = 300)
 cat(sprintf("  -> %s\n", fname))
}

# =============================================================================
# PLOT 2: Per-sample bar plots – full range (log10)
# =============================================================================
cat("Generating per-sample full-range log-scale bar plots...\n")

for (samp in unique(all_data$Sample)) {
 df_s <- filter(all_data, Sample == samp, Frequency > 0)
 col  <- sample_cols[samp]
 alleles <- get_alleles(samp)
 
 p_log <- ggplot(df_s, aes(x = Repeat_Length, y = Frequency)) +
  geom_col(fill = col, color = NA, width = 0.85) +
  scale_x_continuous(breaks = pretty_breaks(n = 10),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_log10(labels = comma,
                expand = expansion(mult = c(0, 0.08))) +
  annotation_logticks(sides = "l", color = "grey50", size = 0.25) +
  labs(title    = samp,
       subtitle = "Distribution of tandem repeat lengths (log\u2081\u2080 scale)",
       x = "Repeat Length (bp)", y = "Frequency (log\u2081\u2080)") +
  theme_pub()
 
 p_log <- add_allele_lines(p_log, alleles)
 
 fname <- file.path(output_dir, paste0(samp, "_barplot_full_log10.pdf"))
 ggsave(fname, p_log, width = 10, height = 5, dpi = 300)
 cat(sprintf("  -> %s\n", fname))
}

# =============================================================================
# PLOT 3: Zoomed bar plots (10–120 and 35–120, linear)
# =============================================================================
cat("Generating zoomed bar plots...\n")

for (samp in unique(all_data$Sample)) {
 df_s <- filter(all_data, Sample == samp)
 col  <- sample_cols[samp]
 alleles <- get_alleles(samp)
 
 for (rng_name in names(zoom_ranges)) {
  rng <- zoom_ranges[[rng_name]]
  df_z <- filter(df_s, Repeat_Length >= rng[1], Repeat_Length <= rng[2])
  
  peak_idx <- find_peaks(df_z$Frequency)
  df_peaks <- df_z[peak_idx, ] %>%
   arrange(desc(Frequency)) %>%
   slice_head(n = 5)
  
  p <- ggplot(df_z, aes(x = Repeat_Length, y = Frequency)) +
   geom_col(fill = col, color = NA, width = 0.85) +
   geom_text(data = df_peaks,
             aes(label = Repeat_Length),
             vjust = -0.5, size = 3, fontface = "bold", color = "grey30") +
   scale_x_continuous(breaks = pretty_breaks(n = 12),
                      expand = expansion(mult = c(0.01, 0.01))) +
   scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
   labs(title    = samp,
        subtitle = sprintf("Repeat lengths %d\u2013%d bp", rng[1], rng[2]),
        x = "Repeat Length (bp)", y = "Frequency") +
   theme_pub()
  
  vis_alleles <- if (!is.null(alleles)) alleles[alleles >= rng[1] & alleles <= rng[2]] else NULL
  p <- add_allele_lines(p, vis_alleles)
  
  fname <- file.path(output_dir, paste0(samp, "_barplot_", rng_name, ".pdf"))
  ggsave(fname, p, width = 10, height = 5, dpi = 300)
  cat(sprintf("  -> %s\n", fname))
 }
}

# =============================================================================
# PLOT 4: Zoomed bar plots (10–120 and 35–120, log10)
# =============================================================================
cat("Generating zoomed log-scale bar plots...\n")

for (samp in unique(all_data$Sample)) {
 df_s <- filter(all_data, Sample == samp, Frequency > 0)
 col  <- sample_cols[samp]
 alleles <- get_alleles(samp)
 
 for (rng_name in names(zoom_ranges)) {
  rng <- zoom_ranges[[rng_name]]
  df_z <- filter(df_s, Repeat_Length >= rng[1], Repeat_Length <= rng[2])
  
  p <- ggplot(df_z, aes(x = Repeat_Length, y = Frequency)) +
   geom_col(fill = col, color = NA, width = 0.85) +
   scale_x_continuous(breaks = pretty_breaks(n = 12),
                      expand = expansion(mult = c(0.01, 0.01))) +
   scale_y_log10(labels = comma,
                 expand = expansion(mult = c(0, 0.08))) +
   annotation_logticks(sides = "l", color = "grey50", size = 0.25) +
   labs(title    = samp,
        subtitle = sprintf("Repeat lengths %d\u2013%d bp (log\u2081\u2080 scale)",
                           rng[1], rng[2]),
        x = "Repeat Length (bp)", y = "Frequency (log\u2081\u2080)") +
   theme_pub()
  
  vis_alleles <- if (!is.null(alleles)) alleles[alleles >= rng[1] & alleles <= rng[2]] else NULL
  p <- add_allele_lines(p, vis_alleles)
  
  fname <- file.path(output_dir, paste0(samp, "_barplot_", rng_name, "_log10.pdf"))
  ggsave(fname, p, width = 10, height = 5, dpi = 300)
  cat(sprintf("  -> %s\n", fname))
 }
}

# =============================================================================
# PLOT 5: Combined dashboard — 4-panel (linear + log10) x (10–120 + 35–120)
# =============================================================================
cat("Generating combined dashboards...\n")

for (samp in unique(all_data$Sample)) {
 df_s <- filter(all_data, Sample == samp, Frequency > 0)
 col  <- sample_cols[samp]
 alleles <- get_alleles(samp)
 
 make_panel <- function(df, rng, scale_type, panel_label) {
  df_z <- filter(df, Repeat_Length >= rng[1], Repeat_Length <= rng[2])
  
  p <- ggplot(df_z, aes(x = Repeat_Length, y = Frequency)) +
   geom_col(fill = col, color = NA, width = 0.85) +
   scale_x_continuous(breaks = pretty_breaks(n = 10),
                      expand = expansion(mult = c(0.01, 0.01))) +
   labs(title = panel_label) +
   theme_pub(base_size = 10) +
   theme(plot.title = element_text(size = 11, face = "bold", hjust = 0))
  
  if (scale_type == "linear") {
   p <- p +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Frequency")
  } else {
   p <- p +
    scale_y_log10(labels = comma,
                  expand = expansion(mult = c(0, 0.08))) +
    annotation_logticks(sides = "l", color = "grey50", size = 0.2) +
    labs(x = "Repeat Length (bp)", y = "Frequency (log\u2081\u2080)")
  }
  
  vis_alleles <- if (!is.null(alleles)) alleles[alleles >= rng[1] & alleles <= rng[2]] else NULL
  p <- add_allele_lines(p, vis_alleles)
  return(p)
 }
 
 pA <- make_panel(df_s, c(10, 120), "linear", "A.  10\u2013120 bp  \u2014  Linear")
 pB <- make_panel(df_s, c(10, 120), "log",    "B.  10\u2013120 bp  \u2014  Log\u2081\u2080")
 pC <- make_panel(df_s, c(35, 120), "linear", "C.  35\u2013120 bp  \u2014  Linear")
 pD <- make_panel(df_s, c(35, 120), "log",    "D.  35\u2013120 bp  \u2014  Log\u2081\u2080")
 
 p_dash <- (pA | pC) / (pB | pD) +
  plot_annotation(
   title    = samp,
   subtitle = "Tandem repeat length frequency distribution",
   theme    = theme(
    plot.title    = element_text(face = "bold", size = 15, hjust = 0.5,
                                 color = "grey15"),
    plot.subtitle = element_text(size = 10.5, hjust = 0.5, color = "grey45",
                                 margin = margin(b = 6))
   )
  )
 
 fname <- file.path(output_dir, paste0(samp, "_dashboard.pdf"))
 ggsave(fname, p_dash, width = 14, height = 9, dpi = 300)
 cat(sprintf("  -> %s\n", fname))
}

# =============================================================================
# PLOT 6: Multi-sample overlay (if >1 sample)
# =============================================================================
if (n_samples > 1) {
 cat("Generating multi-sample overlay plots...\n")
 
 p_overlay <- ggplot(all_data, aes(x = Repeat_Length, y = Frequency,
                                   color = Sample)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  scale_color_manual(values = sample_cols) +
  scale_x_continuous(breaks = pretty_breaks(n = 10),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Repeat Length Frequency Comparison",
       subtitle = paste(n_samples, "samples"),
       x = "Repeat Length (bp)", y = "Frequency", color = "Sample") +
  theme_pub() +
  theme(legend.position = "right")
 
 fname <- file.path(output_dir, "all_samples_overlay.pdf")
 ggsave(fname, p_overlay, width = 12, height = 6, dpi = 300)
 cat(sprintf("  -> %s\n", fname))
 
 # Faceted panel plot
 p_facet <- ggplot(all_data, aes(x = Repeat_Length, y = Frequency,
                                 fill = Sample)) +
  geom_col(color = NA, width = 0.85, show.legend = FALSE) +
  facet_wrap(~ Sample, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = sample_cols) +
  scale_x_continuous(breaks = pretty_breaks(n = 6),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Repeat Length Frequency Distributions",
       x = "Repeat Length (bp)", y = "Frequency") +
  theme_pub()
 
 n_rows <- ceiling(n_samples / 2)
 fname <- file.path(output_dir, "all_samples_faceted.pdf")
 ggsave(fname, p_facet, width = 12, height = 3 + 2.5 * n_rows, dpi = 300)
 cat(sprintf("  -> %s\n", fname))
 
 # Normalised overlay
 all_data_norm <- all_data %>%
  group_by(Sample) %>%
  mutate(Relative_Frequency = Frequency / sum(Frequency)) %>%
  ungroup()
 
 p_norm <- ggplot(all_data_norm, aes(x = Repeat_Length, y = Relative_Frequency,
                                     color = Sample)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  scale_color_manual(values = sample_cols) +
  scale_x_continuous(breaks = pretty_breaks(n = 10),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Normalised Repeat Length Frequency Comparison",
       subtitle = "Relative frequency per sample",
       x = "Repeat Length (bp)", y = "Relative Frequency (%)",
       color = "Sample") +
  theme_pub() +
  theme(legend.position = "right")
 
 fname <- file.path(output_dir, "all_samples_normalised.pdf")
 ggsave(fname, p_norm, width = 12, height = 6, dpi = 300)
 cat(sprintf("  -> %s\n", fname))
}

# --- Done -------------------------------------------------------------------
cat(sprintf("\nAll plots saved to: %s\n", normalizePath(output_dir)))
cat("Done.\n")