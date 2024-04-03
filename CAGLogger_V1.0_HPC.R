#Authors: Hasnahana Chetia; Ivy Kosater; Kert Matlik
#V1.1
#Date:04_01_2024

library(ShortRead)
library(ggplot2)
library(tictoc)

extract_cag_structures_with_modified_errors <- function(seq) {
  pattern <- "(CAG|CA[ATCG]|C[ATCG]G|[ATCG]AG)+"
  matches <- gregexpr(pattern, seq, perl = TRUE)
  structures <- regmatches(seq, matches)  
  has_substitution <- function(triplet) {
    sum(sapply(1:3, function(i) substr(triplet, i, i) != substr("CAG", i, i))) > 0
  }
  validate_triplets <- function(match) {
  # Split the match into triplets without using lookbehind assertions
  triplets <- sapply(seq(1, nchar(match), by=3), function(i) substr(match, i, min(i+2, nchar(match))))
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
process_fastq_file <- function(fastq_path) {
  fastq_data <- readFastq(fastq_path)
  sequences <- sread(fastq_data)
  structure_list <- lapply(as.character(sequences), extract_cag_structures_with_modified_errors)
  structure_tally <- table(unlist(structure_list))
  structure_df <- data.frame(Structure = names(structure_tally), Frequency = as.integer(structure_tally))
  # Total reads
  total_reads <- sum(structure_df$Frequency)
  # Count reads with CAG repeat size over 35 and over 110
  reads_over_35 <- sum(structure_df[structure_df$Structure > 35, "Frequency"])
  reads_over_110 <- sum(structure_df[structure_df$Structure > 110, "Frequency"])
  
  # Calculate percentage of reads with CAG repeat size over 110
  # Check to avoid division by zero
  if (reads_over_35 > 0) {
    percent_over_110 <- (reads_over_110 / reads_over_35) * 100
  } else {
    percent_over_110 <- 0  # or another default value or handling as per your analysis requirement
  }
  
  # Save metrics to a CSV file
  extracted_part <- unlist(strsplit(fastq_path, "_"))[1]
  sample_text <- paste("CAG LOGGER v1.1 Report : SAMPLE ID", extracted_part)
  metrics_df <- data.frame(
    File = sample_text,
    Total_Reads = total_reads,
    Reads_Over_35 = reads_over_35,
    Reads_Over_110 = reads_over_110,
    Percent_Reads_Over_110 = percent_over_110
  )
  metrics_csv_filename <- gsub("\\.fastq.gz$", "_cagloggerv1.1_RepeatMetrics.csv", fastq_path)
  write.csv(metrics_df, metrics_csv_filename, row.names = FALSE)
  csv_filename <- gsub("\\.fastq.gz$", "_cagloggerv1.1_FrequencyOfRepeats.csv", fastq_path)
  write.csv(structure_df, csv_filename, row.names = FALSE)
  }
#working directory
setwd("/ru-auth/local/home/hchetia/caglogger_test")
fastq_files <- list.files(pattern = "\\.fastq.gz$")
tic()
for (file in fastq_files) {
  process_fastq_file(file)
}
toc()
