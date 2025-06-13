library(ShortRead)

#function
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

#fastq_path = full or relative path to fastq file
#lengths = list object containing the structure lengths you wish to pull from fastq file
select_repeat_lengths <- function(fastq_path, lengths){
  seq = fastq_path
  fastq_data <- readFastq(fastq_path)
  sequences <- sread(fastq_data)
  structure_list <- lapply(as.character(sequences), extract_cag_structures_with_modified_errors)
  samples <- c()
  read <- c()
  i <- 1
  
  #Main loop
  while(i < length(structure_list)){
    tmp <- structure_list[i]
    x <- l %in% tmp[[1]]
    y <- TRUE %in% x
    if(y == TRUE){
      samples <- append(samples, tmp)
      read <- append(read, i)
    }
    i <- i+1
    rm(x)
    rm(y)
  }
    
  
  filtered_seq <- sequences[read[1:length(read)]]
  for(i in lengths){
    writeXStringSet(filtered_seq, paste(fastq_path, as.character(i), ".fasta"))
  }
  
}


#Example: using the function to extract reads of length 61 and 62
#select_repeat_lengths("4089_CAG_0109_donor3167_Caudate_DRD2pos_rerun_400_S20_R1_001.fastq.gz",
#                      c(61,62))

