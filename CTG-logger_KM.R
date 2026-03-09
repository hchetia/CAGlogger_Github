extract_ctg_structures_with_modified_errors <- function(seq) {
 pattern <- "(CTG|CT[ATCG]|C[ATCG]G|[ATCG]TG)+"
 matches <- gregexpr(pattern, seq, perl = TRUE)
 structures <- regmatches(seq, matches)  
 
 has_substitution <- function(triplet) {
  sum(sapply(1:3, function(i) substr(triplet, i, i) != substr("CTG", i, i))) > 0
 }
 
 validate_triplets <- function(match) {
  triplets <- sapply(seq(1, nchar(match), by = 3), function(i) substr(match, i, min(i + 2, nchar(match))))
  type_a_error_count <- 0
  type_b_error_count <- 0   
  
  for (i in 1:length(triplets)) {
   if (nchar(triplets[i]) < 3) next  # Skip incomplete triplets
   if (has_substitution(triplets[i])) {
    type_a_error_count <- type_a_error_count + 1
    if (i > 1 && has_substitution(triplets[i - 1])) {
     type_b_error_count <- type_b_error_count + 1
    }
   }
  }
  return(type_a_error_count <= 5 && type_b_error_count <= 1)
 }
 
 sapply(structures, function(match_group) {
  valid_matches <- Filter(validate_triplets, match_group)
  lengths <- nchar(valid_matches) / 3
  lengths[lengths > 0]
 })
}



#gunzip("7992_ATXN3_6114_unsorted_S4_R1_001.fastq.gz")
#fastq_data <- readFastq("7992_ATXN3_6114_unsorted_S4_R1_001.fastq")
fastq_data <- readFastq("atxn3_sim.fastq.gz")
sequences <- sread(fastq_data)

structure_list_TR <- lapply(as.character(sequences), extract_ctg_structures_with_modified_errors)
structure_tally_TR <- table(unlist(structure_list_TR))

structure_df_TR <- data.frame(Structure = names(structure_tally_TR), Frequency = as.integer(structure_tally_TR))
structure_df_TR$Structure<-as.integer(structure_df_TR$Structure)

#whether '-5' or '-6' is subtracted depends on the donor genotype
structure_df_TR$Structure<-structure_df_TR$Structure - 5
df <- tibble(Structure = -5:113)
df3 <- full_join (df, structure_df_TR, by = c ("Structure"))
df3[is.na(df3)] <- 0

structure_df_TR1 <- subset(df3, as.numeric(Structure) >= 1)

structure_df_TR50 <- subset(df3, as.numeric(Structure) >= 50)

write.csv(structure_df_TR1, "oldCTGlogger_results_with_sim.csv", row.names = FALSE)
write.csv(structure_df_TR50, "oldCTGlogger_results_with_sim_over50.csv", row.names = FALSE)