#Authors: Hasnahana Chetia; Ivy Kosater
#V1.0
#Date:01_17_2024

library(ShortRead)
library(ggplot2)
library(tictoc)

#function
extract_cag_structures_with_modified_errors <- function(seq) {
  pattern <- "(CAG|CA[ATCG]|C[ATCG]G|[ATCG]AG)+"
  matches <- gregexpr(pattern, seq, perl = TRUE)
  structures <- regmatches(seq, matches)
  
  has_substitution <- function(triplet) {
    sum(sapply(1:3, function(i) substr(triplet, i, i) != substr("CAG", i, i))) > 0
  }
  validate_triplets <- function(match) {
    triplets <- strsplit(match, "(?<=.{3})", perl = TRUE)[[1]]
    triplets <- sapply(triplets, as.character)
    type_a_error_count <- 0
    type_b_error_count <- 0
    
    for (i in 1:length(triplets)) {
      if (has_substitution(triplets[i])) {
        type_a_error_count <- type_a_error_count + 1
        if (i > 1 && has_substitution(triplets[i - 1])) {
          type_b_error_count <- type_b_error_count + 1
        }
      }
    }
    return(type_a_error_count <= 5 && type_b_error_count <= 1)  # Modified here to allow up to 5 Type A errors
  }
  
  sapply(structures, function(match_group) {
    valid_matches <- Filter(validate_triplets, match_group)
    lengths <- nchar(valid_matches) / 3
    lengths[lengths > 0]
  })
}
#process
process_fastq_file <- function(fastq_path) {
  fastq_data <- readFastq(fastq_path)
  sequences <- sread(fastq_data)
  structure_list <- lapply(as.character(sequences), extract_cag_structures_with_modified_errors)
  structure_tally <- table(unlist(structure_list))
  structure_df <- data.frame(Structure = names(structure_tally), Frequency = as.integer(structure_tally))
  #Total reads
  #total_reads <- sum(structure_df$Frequency)
  #Count reads with CAG repeat size over 35 and over 110
  #reads_over_35 <- sum(structure_df[structure_df$Structure > 35, "Frequency"])
  #reads_over_110 <- sum(structure_df[structure_df$Structure > 110, "Frequency"])
  #Calculate percentage of reads with CAG repeat size over 100
  #percent_over_110 <- (reads_over_110/reads_over_35) * 100
  #Save metrics to a CSV file
  extracted_part <- unlist(strsplit(fastq_path, "_"))[1]
  sample_text <- paste("CAG LOGGER v1.0 Report : SAMPLE ID", extracted_part)
  #metrics_df <- data.frame(
   # File = sample_text,
    #Total_Reads = total_reads,
    #Reads_Over_35 = reads_over_35,
   # Reads_Over_110 = reads_over_110,
   # Percent_Reads_Over_110 = percent_over_110
  #)
  #metrics_csv_filename <- gsub("\\.fastq.gz$", "_cag_repeat_metrics_5typeA.csv", fastq_path)
 # write.csv(metrics_df, metrics_csv_filename, row.names = FALSE)
  csv_filename <- gsub("\\.fastq.gz$", "_cag_structure_counts_5typeA.csv", fastq_path)
  write.csv(structure_df, csv_filename, row.names = FALSE)
  structure_df <- structure_df[structure_df$Frequency > 0, ]
  structure_df$Structure <- as.numeric(as.character(structure_df$Structure))
  structure_df <- structure_df[order(structure_df$Structure), ]
  # Filter out CAG repeats less than 35 in length
  structure_df <- structure_df[structure_df$Structure >= 15, ]
  bar_plot <- ggplot(structure_df, aes(x = Structure, y = Frequency)) +
    geom_bar(stat = "identity", fill="skyblue", colour="black") +
    labs(title = "", x = "(CAG)n Repeat Size", y = "Frequency") + scale_y_continuous(breaks = seq(0, 1000, by=100), limits=c(0,1000)) + 
    scale_x_continuous(breaks = seq(0, 120, by = 2)) +
    theme_classic()
  #range_text <- paste("Range of m(CAG)n is from", min(structure_df$Structure), "to", max(structure_df$Structure), ".")
  #extracted_part <- unlist(strsplit(fastq_path, "_"))[1]
  #sample_text <- paste("SAMPLE ID", extracted_part)
  bar_plot2 <- bar_plot + annotate("text", x = Inf, y = Inf, label = sample_text, hjust = 1.1, vjust = 1.1, size = 6, angle = 0, color = "black")
  png_filename <- gsub("\\.fastq.gz$", "_cag_frequency_plot_5typeA.png", fastq_path)
  ggsave(png_filename, plot = bar_plot2, width = 10, height = 10, device = "png")
}
# Set your working directory if needed
fastq_files <- list.files(pattern = "\\.fastq.gz$")

tic()
for (file in fastq_files) {
  process_fastq_file(file)
}
toc()
