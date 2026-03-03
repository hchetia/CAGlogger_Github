#!/usr/bin/env Rscript

suppressPackageStartupMessages({
 library(optparse)
 library(readr)
 library(dplyr)
 library(stringr)
 library(ggplot2)
 library(ggpubr)
 library(scales)
 library(janitor)
 library(purrr)
})

option_list <- list(
 make_option(c("-i", "--input"), type = "character", default = ".",
             help = "Directory containing *_FrequencyOfRepeats.csv"),
 make_option(c("-o", "--output"), type = "character",
             default = "repeat_frequency_plots",
             help = "Output directory"),
 make_option(c("-p", "--pattern"), type = "character",
             default = "_FrequencyOfRepeats\\.csv$",
             help = "File name pattern"),
 make_option(c("--width"), type = "double", default = 7),
 make_option(c("--height"), type = "double", default = 4.5)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!dir.exists(opt$output)) dir.create(opt$output, recursive = TRUE)

files <- list.files(
 opt$input,
 pattern = opt$pattern,
 recursive = TRUE,
 full.names = TRUE
)

if (length(files) == 0) stop("No matching CSV files found.")

walk(files, function(f) {
 
 df <- read_csv(f, show_col_types = FALSE) |>
  clean_names() |>
  mutate(
   repeat_length = as.integer(repeat_length),
   frequency = as.numeric(frequency)
  ) |>
  arrange(repeat_length)
 
 sample_name <- str_remove(basename(f), opt$pattern)
 total_reads <- sum(df$frequency)
 
 # ===============================
 # Plot 1 : 10–120
 # ===============================
 df1 <- df |> filter(repeat_length >= 10 & repeat_length <= 120)
 ymax1 <- max(df1$frequency)
 
 p1 <- ggplot(df1, aes(repeat_length, frequency)) +
  geom_area(fill = "steelblue", alpha = 0.6) +
  geom_line(color = "black", linewidth = 0.3) +
  geom_vline(xintercept = c(17, 45), linetype = "dashed") +
  geom_text(
   data = df1 |> filter(repeat_length %in% c(17,45)),
   aes(label = paste0("Allele ", repeat_length)),
   vjust = -0.5,
   size = 3.5
  ) +
  scale_y_continuous(labels = comma) +
  coord_cartesian(xlim = c(10,120)) +
  labs(
   title = paste0(sample_name, " — Repeat Length Frequency (10–120)"),
   x = "Repeat Length",
   y = "Frequency",
   caption = paste0("Total reads: ", comma(total_reads))
  ) +
  theme_classic()
 
 p1 <- ggpar(p1, font.title = 12, font.x = 11, font.y = 11)
 
 # ===============================
 # Plot 2 : 36–120
 # ===============================
 df2 <- df |> filter(repeat_length >= 36 & repeat_length <= 120)
 ymax2 <- max(df2$frequency)
 
 p2 <- ggplot(df2, aes(repeat_length, frequency)) +
  geom_area(fill = "darkorange", alpha = 0.6) +
  geom_line(color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 45, linetype = "dashed") +
  geom_text(
   data = df2 |> filter(repeat_length == 45),
   aes(label = "Allele 45"),
   vjust = -0.5,
   size = 3.5
  ) +
  annotate(
   "text",
   x = 38,
   y = ymax2 * 0.85,
   label = "Allele 17 outside range",
   hjust = 0,
   size = 3.5
  ) +
  scale_y_continuous(labels = comma) +
  coord_cartesian(xlim = c(36,120)) +
  labs(
   title = paste0(sample_name, " — Repeat Length Frequency (36–120)"),
   x = "Repeat Length",
   y = "Frequency",
   caption = paste0("Total reads: ", comma(total_reads))
  ) +
  theme_classic()
 
 p2 <- ggpar(p2, font.title = 12, font.x = 11, font.y = 11)
 
 # ===============================
 # Save PDF outputs only
 # ===============================
 ggsave(
  filename = file.path(opt$output,
                       paste0(sample_name, "_RepeatFreq_10-120.pdf")),
  plot = p1,
  width = opt$width,
  height = opt$height,
  device = cairo_pdf
 )
 
 ggsave(
  filename = file.path(opt$output,
                       paste0(sample_name, "_RepeatFreq_36-120.pdf")),
  plot = p2,
  width = opt$width,
  height = opt$height,
  device = cairo_pdf
 )
 
 message("Saved plots for: ", basename(f))
})

message("Done.")