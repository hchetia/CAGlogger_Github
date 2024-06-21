#ctg-logger
#version: V1.1
#Authors: Hasnahana Chetia; Ivy Kosater; Kert Matlik
#Date:06-18-2024
library(ShortRead)
library(tictoc)

#working hasna cag function
extract_ctg_structures_with_modified_errors <- function(seq) {
  pattern <- "(CTG|CT[ATCG]|C[ATCG]G|[ATCG]TG)+"
  matches <- gregexpr(pattern, seq, perl = TRUE)
  structures <- regmatches(seq, matches)  
  has_substitution <- function(triplet) {
    sum(sapply(1:3, function(i) substr(triplet, i, i) != substr("CTG", i, i))) > 0
  }
  validate_triplets <- function(match) {
    # Split the match into triplets without using look behind assertions
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
    return(type_a_error_count <= 3 && type_b_error_count <= 1)  # Modified here to allow up to 3 Type A errors and 1 type B errors
  }  
  sapply(structures, function(match_group) {
    valid_matches <- Filter(validate_triplets, match_group)
    lengths <- nchar(valid_matches) / 3
    lengths[lengths > 0]
  })
}
#processing
process_fastq_file <- function(fastq_path) {
    fastq_data <- readFastq(fastq_path)
    #fastq_data <- readFasta("hc_dummy_ctgrepeats.fasta")  
    sequences <- sread(fastq_data)
    structure_list <- lapply(as.character(sequences), extract_ctg_structures_with_modified_errors)
    structure_tally <- table(unlist(structure_list))
    structure_df <- data.frame(Structure = names(structure_tally), Frequency = as.integer(structure_tally))
    structure_df <- subset(structure_df, as.numeric(Structure) >= 5)
    #structure_df
    #structure_df$Structure <- as.numeric(structure_df$Structure) - 3
    structure_df <- subset(structure_df, Structure >= 5)
    # Counting reads with CAG repeat size over 5 and over 110
    reads_over_5 <- sum(structure_df[structure_df$Structure > 5, "Frequency"])
    reads_over_110 <- sum(structure_df[structure_df$Structure > 110, "Frequency"])
    #Metrics Calculation
    #Precheck to avoid division by zero
    if (reads_over_5 > 0) {
      percent_over_110 <- (reads_over_110 / reads_over_5) * 100
    } else {
      percent_over_110 <- NA 
    }
    #Saving metrics file to a CSV file
    extracted_part <- unlist(strsplit(fastq_path, "_"))[1]
    sample_text <- paste("CTG-LOGGER v1.1 Report : SAMPLE ID", extracted_part)
    metrics_df <- data.frame(
      File = sample_text,
      Reads_Over_5 = reads_over_5,
      Reads_Over_110 = reads_over_110,
      Percent_Reads_Over_110 = percent_over_110
    )
    metrics_csv_filename <- gsub("\\.fastq.gz$", "_CTGLogger.v1.1_RepeatMetrics.csv", fastq_path)
    write.csv(metrics_df, metrics_csv_filename, row.names = FALSE)
    csv_filename <- gsub("\\.fastq.gz$", "_CTGLogger.v1.1_FrequencyOfRepeats.csv", fastq_path)
    write.csv(structure_df, csv_filename, row.names = FALSE)
  }
  fastq_files <- list.files(pattern = "\\.fastq.gz$")
tic()
for (file in fastq_files) {
  process_fastq_file(file)
}
toc()
