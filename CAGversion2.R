# CAG LOGGER v2.0
# Alignment-free estimation of CAG repeats from targeted amplicon sequencing data
#
# Changes from v1.2 (peer review response):
#
#   1. -3 correction now conditional on boundary_mode (only applied in "none" mode)
#   2. Reverse-complement handling — reads are searched on both strands
#   3. Fixed mode calculation bug (was calling table() on a frequency table)
#   4. Improved regex specificity with configurable max error tolerance
#   5. Basic allele calling via kernel density peak detection for het samples
#   6. FastqStreamer support for memory-efficient processing of large files
#   7. Configurable sample_id parameter (no longer relies on filename convention)
#   8. Version strings unified to v2.0 throughout
#   9. Incomplete triplet edge case handled explicitly
#  10. Repeat-length floor is now a user parameter (min_repeat_length)
#  11. Per-read accounting: each read contributes its single longest match only
#  12. Synthetic validation function included

# ── Libraries ──────────────────────────────────────────────────────────────────

library(ShortRead)
library(tictoc)

VERSION <- "2.0"

# ── Reverse complement ────────────────────────────────────────────────────────
# Suggestion #3: MiSeq amplicon reads come in both orientations.
# We must search both the forward read (CAG) and its reverse complement (CTG).

rc <- function(seq) {
 chartr("ACGTacgt", "TGCAtgca", paste(rev(strsplit(seq, "")[[1]]), collapse = ""))
}

# ── Quality filtering helper ───────────────────────────────────────────────────

passes_quality_filter <- function(phred_quality, q_threshold = 20) {
 qual_matrix <- as(phred_quality, "matrix")
 mean_quals  <- rowMeans(qual_matrix, na.rm = TRUE)
 mean_quals >= q_threshold
}

# ── Repeat-region boundary helper ─────────────────────────────────────────────

apply_boundary <- function(seq,
                           boundary_mode = c("none", "trim", "mask"),
                           left_anchor   = NULL,
                           right_anchor  = NULL) {
 
 boundary_mode <- match.arg(boundary_mode)
 
 if (boundary_mode == "none" || is.null(left_anchor) || is.null(right_anchor)) {
  return(seq)
 }
 
 left_pos  <- regexpr(left_anchor,  seq, fixed = TRUE)
 right_pos <- regexpr(right_anchor, seq, fixed = TRUE)
 
 if (left_pos == -1L || right_pos == -1L) return(seq)
 
 repeat_start <- left_pos  + nchar(left_anchor)
 repeat_end   <- right_pos - 1L
 
 if (repeat_start > repeat_end) return(seq)
 
 if (boundary_mode == "trim") {
  return(substr(seq, repeat_start, repeat_end))
 }
 
 if (boundary_mode == "mask") {
  left_mask  <- paste(rep("N", repeat_start - 1L),       collapse = "")
  right_mask <- paste(rep("N", nchar(seq) - repeat_end), collapse = "")
  return(paste0(left_mask,
                substr(seq, repeat_start, repeat_end),
                right_mask))
 }
}

# ── CAG structure extractor (improved) ────────────────────────────────────────
# Suggestion #2:  Tightened error thresholds; now configurable.
# Suggestion #9:  Incomplete trailing triplets (nchar not divisible by 3)
#                 are explicitly dropped before error counting.
# Suggestion #11: Returns only the LONGEST valid match per read to avoid
#                 inflating counts with multiple sub-matches from one read.

extract_cag_from_read <- function(seq,
                                  max_type_a = 5,
                                  max_type_b = 1) {
 
 pattern  <- "(CAG|CA[ATCG]|C[ATCG]G|[ATCG]AG)+"
 matches  <- gregexpr(pattern, seq, perl = TRUE)
 raw_hits <- regmatches(seq, matches)[[1]]
 
 if (length(raw_hits) == 0) return(NA_integer_)
 
 has_substitution <- function(triplet) {
  sum(vapply(1:3, function(i)
   substr(triplet, i, i) != substr("CAG", i, i), logical(1))) > 0
 }
 
 validate_triplets <- function(match_str) {
  # ── Suggestion #9: drop incomplete trailing triplet ──
  usable_len <- (nchar(match_str) %/% 3) * 3
  if (usable_len < 3) return(NA_integer_)
  match_str <- substr(match_str, 1, usable_len)
  
  triplets <- vapply(seq(1, usable_len, by = 3),
                     function(i) substr(match_str, i, i + 2),
                     character(1))
  
  type_a <- 0L
  type_b <- 0L
  prev_was_error <- FALSE
  
  for (i in seq_along(triplets)) {
   is_error <- has_substitution(triplets[i])
   if (is_error) {
    type_a <- type_a + 1L
    if (prev_was_error) type_b <- type_b + 1L
   }
   prev_was_error <- is_error
  }
  
  if (type_a <= max_type_a && type_b <= max_type_b) {
   return(length(triplets))
  }
  NA_integer_
 }
 
 lengths <- vapply(raw_hits, validate_triplets, integer(1))
 lengths <- lengths[!is.na(lengths)]
 
 # ── Suggestion #11: return single longest match per read ──
 if (length(lengths) == 0) return(NA_integer_)
 max(lengths)
}

# ── Strand-aware wrapper ──────────────────────────────────────────────────────
# Suggestion #3: Search both orientations, keep the longer match.

extract_cag_both_strands <- function(seq, max_type_a = 5, max_type_b = 1) {
 fwd <- extract_cag_from_read(seq,     max_type_a, max_type_b)
 rev <- extract_cag_from_read(rc(seq), max_type_a, max_type_b)
 
 fwd[is.na(fwd)] <- 0L
 rev[is.na(rev)] <- 0L
 
 best <- ifelse(fwd >= rev, fwd, rev)
 best[best == 0L] <- NA_integer_
 best
}

# ── Allele calling (basic) ────────────────────────────────────────────────────
# Suggestion #5: Kernel density peak detection to identify 1 or 2 alleles
# in the repeat-length distribution. Reports peak positions, estimated
# proportion per allele, and a het/hom call.

call_alleles <- function(repeat_lengths, frequencies, bw = 1.5) {
 
 if (length(repeat_lengths) < 3 || sum(frequencies) < 20) {
  return(data.frame(Allele = NA, Peak_Repeat = NA, Proportion = NA,
                    Call = "insufficient_data", stringsAsFactors = FALSE))
 }
 
 # Expand to individual observations for density estimation
 obs <- rep(repeat_lengths, times = frequencies)
 d   <- density(obs, bw = bw, n = 1024)
 
 # Detect peaks (local maxima)
 dy     <- diff(d$y)
 peaks  <- which(dy[-length(dy)] > 0 & dy[-1] < 0) + 1
 peak_x <- d$x[peaks]
 peak_y <- d$y[peaks]
 
 # Keep only peaks with height >= 10% of the tallest peak
 keep   <- peak_y >= 0.10 * max(peak_y)
 peak_x <- peak_x[keep]
 peak_y <- peak_y[keep]
 
 if (length(peak_x) == 0) {
  return(data.frame(Allele = NA, Peak_Repeat = NA, Proportion = NA,
                    Call = "no_peaks", stringsAsFactors = FALSE))
 }
 
 # Sort by height (descending); take top 2 at most
 ord    <- order(peak_y, decreasing = TRUE)
 peak_x <- round(peak_x[ord][1:min(2, length(ord))])
 peak_x <- sort(unique(peak_x))
 
 # Assign each observation to its nearest peak → estimate proportions
 assignments <- vapply(obs, function(o) peak_x[which.min(abs(o - peak_x))], numeric(1))
 props <- table(assignments) / length(obs)
 
 call <- if (length(peak_x) == 1) "homozygous" else "heterozygous"
 
 data.frame(
  Allele      = seq_along(peak_x),
  Peak_Repeat = peak_x,
  Proportion  = as.numeric(props),
  Call        = call,
  stringsAsFactors = FALSE
 )
}

# ── Per-file processor ────────────────────────────────────────────────────────
# Suggestion #6:  Uses FastqStreamer for memory-efficient processing.
# Suggestion #7:  Accepts a configurable sample_id.
# Suggestion #10: min_repeat_length is now a parameter (default 0 to keep all).

process_fastq_file <- function(fastq_path,
                               sample_id        = NULL,
                               q_threshold      = 20,
                               boundary_mode    = c("none", "trim", "mask"),
                               left_anchor      = NULL,
                               right_anchor     = NULL,
                               min_repeat_length = 0,
                               flank_correction = "auto",
                               max_type_a       = 5,
                               max_type_b       = 1,
                               chunk_size       = 1e6,
                               allele_bw        = 1.5) {
 
 boundary_mode <- match.arg(boundary_mode)
 
 # ── Suggestion #7: configurable sample ID ──────────────────────────────────
 if (is.null(sample_id)) {
  # Fallback: strip path + extension
  sample_id <- sub("\\.(fastq|fq)(\\.gz)?$", "",  basename(fastq_path))
 }
 
 # ── Suggestion #1: conditional -3 correction ───────────────────────────────
 # When boundary trimming/masking is active, flanking bleed-in is already
 # removed, so the empirical -3 correction should NOT be applied.
 if (flank_correction == "auto") {
  correction <- if (boundary_mode == "none") 3L else 0L
 } else {
  correction <- as.integer(flank_correction)
 }
 
 # ── Streaming reader (Suggestion #6) ───────────────────────────────────────
 streamer    <- FastqStreamer(fastq_path, n = chunk_size)
 on.exit(close(streamer))
 
 total_reads    <- 0L
 reads_passing  <- 0L
 all_lengths    <- integer(0)  # one entry per read (longest match)
 
 repeat {
  chunk <- yield(streamer)
  if (length(chunk) == 0) break
  
  n <- length(chunk)
  total_reads <- total_reads + n
  
  # ── Quality filter ─────────────────────────────────────────────────────
  keep  <- passes_quality_filter(quality(chunk), q_threshold)
  chunk <- chunk[keep]
  reads_passing <- reads_passing + sum(keep)
  
  if (length(chunk) == 0) next
  
  # ── Boundary processing ────────────────────────────────────────────────
  seq_char <- as.character(sread(chunk))
  
  if (boundary_mode != "none" && !is.null(left_anchor) && !is.null(right_anchor)) {
   seq_char <- vapply(seq_char,
                      apply_boundary,
                      FUN.VALUE = character(1),
                      boundary_mode = boundary_mode,
                      left_anchor   = left_anchor,
                      right_anchor  = right_anchor)
  }
  
  # ── Extract CAG (both strands — Suggestion #3) ─────────────────────────
  chunk_lengths <- vapply(seq_char,
                          extract_cag_both_strands,
                          FUN.VALUE = integer(1),
                          max_type_a = max_type_a,
                          max_type_b = max_type_b)
  
  all_lengths <- c(all_lengths, chunk_lengths)
 }
 
 reads_failed <- total_reads - reads_passing
 
 message(sprintf("[%s] Total reads: %d | Passed Q>=%d: %d | Removed: %d",
                 sample_id, total_reads, q_threshold, reads_passing, reads_failed))
 
 if (reads_passing == 0) {
  warning(sprintf("No reads passed quality filter for %s. Skipping.", sample_id))
  return(invisible(NULL))
 }
 
 # ── Remove NAs (reads with no valid CAG match) ────────────────────────────
 all_lengths <- all_lengths[!is.na(all_lengths)]
 if (length(all_lengths) == 0) {
  warning(sprintf("No CAG structures detected in %s. Skipping.", sample_id))
  return(invisible(NULL))
 }
 
 # ── Apply correction (Suggestion #1) ──────────────────────────────────────
 all_lengths <- all_lengths - correction
 all_lengths <- all_lengths[all_lengths > 0]
 
 # ── Build frequency table ─────────────────────────────────────────────────
 freq_table   <- as.data.frame(table(all_lengths), stringsAsFactors = FALSE)
 names(freq_table) <- c("Repeat_Length", "Frequency")
 freq_table$Repeat_Length <- as.integer(freq_table$Repeat_Length)
 
 # ── Full distribution (before floor) — written for transparency ───────────
 full_freq <- freq_table
 
 # ── Apply user-configurable floor (Suggestion #10) ─────────────────────────
 freq_table <- freq_table[freq_table$Repeat_Length >= min_repeat_length, ]
 
 if (nrow(freq_table) == 0) {
  warning(sprintf("No repeats >= %d after correction in %s.", min_repeat_length, sample_id))
  return(invisible(NULL))
 }
 
 # ── Suggestion #4: Fixed mode calculation ──────────────────────────────────
 mode_repeat <- freq_table$Repeat_Length[which.max(freq_table$Frequency)]
 
 # ── Threshold metrics ──────────────────────────────────────────────────────
 reads_over_35  <- sum(freq_table$Frequency[freq_table$Repeat_Length >= 36])
 reads_over_110 <- sum(freq_table$Frequency[freq_table$Repeat_Length >  110])
 
 percent_over_110 <- if (reads_over_35 > 0) {
  (reads_over_110 / reads_over_35) * 100
 } else {
  NA_real_
 }
 
 # ── Suggestion #5: allele calling ──────────────────────────────────────────
 allele_df <- call_alleles(freq_table$Repeat_Length,
                           freq_table$Frequency,
                           bw = allele_bw)
 
 # ── Write outputs (Suggestion #8: version string v2.0) ─────────────────────
 sample_text <- paste("CAG LOGGER v", VERSION, " Report : SAMPLE ID ", sample_id, sep = "")
 
 metrics_df <- data.frame(
  Sample_ID               = sample_text,
  Total_Reads             = total_reads,
  Reads_Passed_QFilter    = reads_passing,
  Reads_Failed_QFilter    = reads_failed,
  Reads_With_CAG          = length(all_lengths),
  Q_Threshold             = q_threshold,
  Boundary_Mode           = boundary_mode,
  Flank_Correction        = correction,
  Min_Repeat_Length       = min_repeat_length,
  Mode_Repeat_Size        = mode_repeat,
  Reads_Over_35           = reads_over_35,
  Reads_Over_110          = reads_over_110,
  Percent_Reads_Over_110  = percent_over_110,
  Allele_Call             = allele_df$Call[1],
  Allele_1_Peak           = allele_df$Peak_Repeat[1],
  Allele_1_Proportion     = allele_df$Proportion[1],
  Allele_2_Peak           = ifelse(nrow(allele_df) >= 2, allele_df$Peak_Repeat[2], NA),
  Allele_2_Proportion     = ifelse(nrow(allele_df) >= 2, allele_df$Proportion[2], NA)
 )
 
 base_name         <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(fastq_path))
 metrics_csv       <- paste0(base_name, "_CAGLogger.v", VERSION, "_RepeatMetrics.csv")
 freq_csv          <- paste0(base_name, "_CAGLogger.v", VERSION, "_FrequencyOfRepeats.csv")
 full_freq_csv     <- paste0(base_name, "_CAGLogger.v", VERSION, "_FullDistribution.csv")
 allele_csv        <- paste0(base_name, "_CAGLogger.v", VERSION, "_AlleleCalls.csv")
 
 write.csv(metrics_df, metrics_csv,   row.names = FALSE)
 write.csv(freq_table, freq_csv,      row.names = FALSE)
 write.csv(full_freq,  full_freq_csv, row.names = FALSE)
 write.csv(allele_df,  allele_csv,    row.names = FALSE)
 
 message(sprintf("[%s] Mode: %d | Allele call: %s | Reads >35: %d | Reads >110: %d (%.1f%%)",
                 sample_id, mode_repeat, allele_df$Call[1],
                 reads_over_35, reads_over_110,
                 ifelse(is.na(percent_over_110), 0, percent_over_110)))
 
 # Return results invisibly for programmatic use
 invisible(list(
  sample_id  = sample_id,
  metrics    = metrics_df,
  freq_table = freq_table,
  full_freq  = full_freq,
  alleles    = allele_df
 ))
}

# ══════════════════════════════════════════════════════════════════════════════
# Suggestion #12: Synthetic validation suite
# ══════════════════════════════════════════════════════════════════════════════
# Generates FASTQ records with known repeat lengths and runs them through the
# pipeline. Reports observed vs. expected to verify accuracy.

generate_synthetic_fastq <- function(outfile       = "synthetic_validation.fastq",
                                     repeat_lengths = c(20, 36, 42, 70, 120),
                                     reads_per_length = 100,
                                     error_rate     = 0.005,
                                     flanking_left  = "GCGACCCTGGAAAAGCTGA",
                                     flanking_right = "GGCGGCTGAGGAAGCTGAG",
                                     seed           = 42) {
 set.seed(seed)
 bases <- c("A", "T", "C", "G")
 
 lines <- character(0)
 
 read_num <- 0
 for (rl in repeat_lengths) {
  for (i in seq_len(reads_per_length)) {
   read_num <- read_num + 1
   
   # Build perfect repeat
   repeat_seq <- paste(rep("CAG", rl), collapse = "")
   
   # Introduce substitution errors at the given rate
   repeat_chars <- strsplit(repeat_seq, "")[[1]]
   n_errors     <- rbinom(1, length(repeat_chars), error_rate)
   if (n_errors > 0) {
    err_pos <- sample(seq_along(repeat_chars), n_errors)
    for (ep in err_pos) {
     orig <- repeat_chars[ep]
     repeat_chars[ep] <- sample(setdiff(bases, orig), 1)
    }
   }
   repeat_seq <- paste(repeat_chars, collapse = "")
   
   full_seq <- paste0(flanking_left, repeat_seq, flanking_right)
   
   # Randomly reverse-complement ~50% of reads
   if (runif(1) < 0.5) {
    full_seq <- rc(full_seq)
   }
   
   # Quality line: uniform Q30
   qual <- paste(rep("?", nchar(full_seq)), collapse = "")  # ASCII 63 = Q30
   
   lines <- c(lines,
              paste0("@SYNTH_", read_num, "_truelen_", rl),
              full_seq,
              "+",
              qual)
  }
 }
 
 writeLines(lines, outfile)
 message(sprintf("Wrote %d synthetic reads to %s", read_num, outfile))
 invisible(outfile)
}

run_validation <- function(error_rate = 0.005, seed = 42) {
 
 message("\n══════════════════════════════════════════════")
 message("  CAG LOGGER v", VERSION, " — Validation Suite")
 message("══════════════════════════════════════════════\n")
 
 expected_lengths <- c(20, 36, 42, 70, 120)
 reads_per        <- 100
 
 tmpfile <- tempfile(fileext = ".fastq")
 
 generate_synthetic_fastq(
  outfile          = tmpfile,
  repeat_lengths   = expected_lengths,
  reads_per_length = reads_per,
  error_rate       = error_rate,
  seed             = seed
 )
 
 result <- process_fastq_file(
  fastq_path        = tmpfile,
  sample_id         = "VALIDATION",
  q_threshold       = 10,         # synthetic reads are high-quality
  boundary_mode     = "trim",
  left_anchor       = "GCGACCCTGGAAAAGCTGA",
  right_anchor      = "GGCGGCTGAGGAAGCTGAG",
  min_repeat_length = 0
 )
 
 if (is.null(result)) {
  message("VALIDATION FAILED: no results returned.")
  return(invisible(FALSE))
 }
 
 freq <- result$full_freq
 
 message("\n── Validation Results ──────────────────────────")
 message(sprintf("%-12s %-12s %-12s %-10s", "Expected", "Observed_Peak", "Reads", "Status"))
 
 all_pass <- TRUE
 for (el in expected_lengths) {
  # Find the peak within ±2 of expected
  window <- freq[freq$Repeat_Length >= (el - 2) & freq$Repeat_Length <= (el + 2), ]
  if (nrow(window) == 0) {
   message(sprintf("%-12d %-12s %-12s %-10s", el, "MISSING", "0", "FAIL"))
   all_pass <- FALSE
   next
  }
  peak_row <- window[which.max(window$Frequency), ]
  delta    <- abs(peak_row$Repeat_Length - el)
  status   <- ifelse(delta <= 1, "PASS", "WARN")
  if (delta > 2) { status <- "FAIL"; all_pass <- FALSE }
  message(sprintf("%-12d %-12d %-12d %-10s",
                  el, peak_row$Repeat_Length, peak_row$Frequency, status))
 }
 
 message(sprintf("\nOverall: %s", ifelse(all_pass, "ALL PASSED", "ISSUES DETECTED")))
 
 unlink(tmpfile)
 invisible(all_pass)
}

# ── Run ───────────────────────────────────────────────────────────────────────

# To run validation:
#run_validation()

# To process real data:

#fastq_files <- list.files(pattern = "\\.(fastq|fq)(\\.gz)?$")
#
tic()
results <- list()
for (file in fastq_files) {
  results[[file]] <- process_fastq_file(
    fastq_path        = file,
    sample_id         = NULL,         # auto-derived from filename
    q_threshold       = 0, #can change to 20 or 40
    boundary_mode     = NULL,       # "trim" or "mask" with anchors below
    left_anchor       = NULL,         # e.g. "GCGACCCTGGAAAAGCTGA"
    right_anchor      = NULL,         # e.g. "GGCGGCTGAGGAAGCTGAG"
    min_repeat_length = 10,            # set to 35 to replicate v1.2 behaviour
    flank_correction  = "auto",       # "auto", or integer (e.g. 0, 3)
    max_type_a        = 5,
    max_type_b        = 1,
    chunk_size        = 3e6,
    allele_bw         = 1.5
  )
}
toc()