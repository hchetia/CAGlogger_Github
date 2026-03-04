#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# CAG LOGGER v2.1c — R implementation
# Alignment-free estimation of CAG repeats from targeted amplicon MiSeq data
#
# Counting method:
#   Only pure "CAG" triplets within each accepted regex match are counted.
#   Non-CAG triplets (sequencing errors, flank elements) are tolerated for
#   match continuity but are NOT included in the reported repeat length.
#   Type A (total non-CAG) and Type B (consecutive non-CAG) error thresholds
#   gate whether a match is accepted.
#
# Changes from v2.1:
#
#   Fix G — Length-adaptive Type A threshold: fixed max_type_a caused
#           systematic rejection of long expanded-allele reads because
#           sequencing errors scale with tract length.  Type A threshold
#           now scales as:
#
#             effective_max_type_a = max(max_type_a, floor(n_triplets * type_a_rate))
#
#           where max_type_a is the fixed floor (default 5) and type_a_rate
#           is the per-triplet error tolerance (default 0.05, i.e. 5%).
#           Short matches are still gated by the fixed floor; long matches
#           get proportionally more room.
#
#   Diagnostic — threshold_rejected counter: reads whose longest regex
#           match exceeded all error thresholds are now counted and their
#           raw triplet count is recorded for post-hoc review.
#
# Previous changes (v2.1):
#   Fix A — N-aware regex
#   Fix F — CCG tract trimming
#   From v2.0: reverse-complement, fixed mode, configurable params,
#              per-read single longest match, streaming FASTQ, allele calling
#
# Dependencies: ShortRead + Biostrings (required).
# ══════════════════════════════════════════════════════════════════════════════

VERSION <- "2.1c" #typeAerrorsdynamic(for long repeats ~6)

# ── Load required packages ───────────────────────────────────────────────────

if (!requireNamespace("ShortRead", quietly = TRUE)) {
 stop("ShortRead package is required. Install via BiocManager::install('ShortRead')")
}

library(ShortRead)
library(Biostrings)

# ── Reverse complement (both strands searched) ───────────────────────────────

rc <- function(seq) {
 chartr("ACGTacgtN", "TGCAtgcaN",
        paste(rev(strsplit(seq, "")[[1]]), collapse = ""))
}

# ── CCG-like detector (used by Fix F) ────────────────────────────────────────

is_ccg_like <- function(triplet) {
 # TRUE if triplet is CCG or a 1-substitution/N variant, but NOT CAG itself
 if (triplet == "CAG") return(FALSE)
 chars <- strsplit(triplet, "")[[1]]
 ref   <- c("C", "C", "G")
 mismatches <- sum(chars != ref & chars != "N")
 mismatches <= 1L
}

# ── CAG extraction (N-aware, adaptive correction) ───────────────────────────
#
# Fix A: regex includes N at each mismatch position.

CAG_PATTERN <- "(CAG|CA[ATCGN]|C[ATCGN]G|[ATCGN]AG)+"

#' Extract the longest valid CAG match from a single sequence.
#' The reported length is the count of pure "CAG" triplets only.
#' Returns a named list: length (integer or NA), hits_end, rejected_triplets.
#'   rejected_triplets: if the best raw match was rejected by thresholds,
#'                      this is the triplet count of that match (for diagnostics).
#'                      NA if a match was accepted or no match existed.
extract_cag_from_read <- function(seq,
                                  max_type_a     = 5L,
                                  max_type_b     = 1L,
                                  type_a_rate    = 0.05) {
 
 matches <- gregexpr(CAG_PATTERN, seq, perl = TRUE)
 raw_hits <- regmatches(seq, matches)[[1]]
 
 if (length(raw_hits) == 0L) {
  return(list(length = NA_integer_, hits_end = FALSE,
              rejected_triplets = NA_integer_))
 }
 
 starts   <- as.integer(matches[[1]])
 read_len <- nchar(seq)
 
 best_length            <- NA_integer_
 best_hits_end          <- FALSE
 longest_rejected_count <- NA_integer_    # track largest rejected match
 
 for (idx in seq_along(raw_hits)) {
  match_str  <- raw_hits[idx]
  usable_len <- (nchar(match_str) %/% 3L) * 3L
  if (usable_len < 3L) next
  match_str <- substr(match_str, 1L, usable_len)
  
  triplets <- vapply(
   seq(1L, usable_len, by = 3L),
   function(i) substr(match_str, i, i + 2L),
   character(1)
  )
  
  # ── Fix F: Trim absorbed CCG tract ─────────────────────────────────
  consec_ccg <- 0L
  trim_at    <- length(triplets) + 1L
  for (ti in seq_along(triplets)) {
   if (is_ccg_like(triplets[ti])) {
    consec_ccg <- consec_ccg + 1L
    if (consec_ccg >= 2L) {
     trim_at <- ti - consec_ccg + 1L
     break
    }
   } else {
    consec_ccg <- 0L
   }
  }
  if (trim_at <= length(triplets)) {
   triplets   <- triplets[seq_len(trim_at - 1L)]
   usable_len <- length(triplets) * 3L
  }
  if (length(triplets) < 1L) next
  
  # ── Error counting ─────────────────────────────────────────────────
  type_a   <- 0L
  type_b   <- 0L
  prev_err <- FALSE
  
  for (t in triplets) {
   is_err <- (t != "CAG")
   
   if (is_err) {
    type_a <- type_a + 1L
    if (prev_err) type_b <- type_b + 1L
   }
   prev_err <- is_err
  }
  
  # ── Fix G: Length-adaptive Type A threshold ────────────────────────
  #
  # The fixed max_type_a works well for short matches but systematically
  # rejects long expanded-allele reads where the absolute number of
  # sequencing errors is expected to be higher.
  #
  # effective_max_type_a = max(fixed_floor, floor(n_triplets * rate))
  #
  # Examples with defaults (floor=5, rate=0.05):
  #   20 triplets  -> max(5, 1) = 5   (floor dominates)
  #   80 triplets  -> max(5, 4) = 5   (floor dominates)
  #  120 triplets  -> max(5, 6) = 6   (rate dominates)
  #  200 triplets  -> max(5, 10) = 10 (rate dominates)
  
  n_triplets           <- length(triplets)
  effective_max_type_a <- max(max_type_a,
                              as.integer(floor(n_triplets * type_a_rate)))
  
  if (type_a > effective_max_type_a || type_b > max_type_b) {
   # Record the size of the rejected match for diagnostics
   if (is.na(longest_rejected_count) || n_triplets > longest_rejected_count) {
    longest_rejected_count <- n_triplets
   }
   next
  }
  
  # ── Pure CAG count ─────────────────────────────────────────────────
  pure_cag <- sum(triplets == "CAG")
  
  # ── Flank CAG exclusion ────────────────────────────────────────────
  n_tr <- length(triplets)
  if (n_tr >= 2L) {
   last    <- triplets[n_tr]
   penult  <- triplets[n_tr - 1L]
   is_caa  <- penult %in% c("CAA", "CAN", "NAA")
   
   if (is_caa && last == "CAG") {
    pure_cag <- pure_cag - 1L
   } else if (is_ccg_like(last) && n_tr >= 3L) {
    ante   <- triplets[n_tr - 2L]
    middle <- triplets[n_tr - 1L]
    if (ante %in% c("CAA", "CAN", "NAA") && middle == "CAG") {
     pure_cag <- pure_cag - 1L
    }
   }
  }
  
  # ── Does the match reach the end of the read? ──────────────────────
  match_end_pos <- starts[idx] + usable_len - 1L
  hits_read_end <- (read_len - match_end_pos) <= 3L
  
  if (is.na(best_length) || pure_cag > best_length) {
   best_length   <- pure_cag
   best_hits_end <- hits_read_end
  }
 }
 
 # If a match was accepted, clear the rejected tracker
 rejected_out <- if (!is.na(best_length)) NA_integer_ else longest_rejected_count
 
 list(length = best_length, hits_end = best_hits_end,
      rejected_triplets = rejected_out)
}

#' Search both orientations, return best result.
extract_cag_both_strands <- function(seq,
                                     max_type_a  = 5L,
                                     max_type_b  = 1L,
                                     type_a_rate = 0.05) {
 
 fwd <- extract_cag_from_read(seq,     max_type_a, max_type_b, type_a_rate)
 rev <- extract_cag_from_read(rc(seq), max_type_a, max_type_b, type_a_rate)
 
 fwd_len <- ifelse(is.na(fwd$length), 0L, fwd$length)
 rev_len <- ifelse(is.na(rev$length), 0L, rev$length)
 
 if (fwd_len >= rev_len) {
  return(list(
   length            = if (fwd_len > 0L) fwd_len else NA_integer_,
   hits_end          = fwd$hits_end,
   rejected_triplets = if (fwd_len > 0L) NA_integer_ else fwd$rejected_triplets
  ))
 } else {
  return(list(
   length            = if (rev_len > 0L) rev_len else NA_integer_,
   hits_end          = rev$hits_end,
   rejected_triplets = if (rev_len > 0L) NA_integer_ else rev$rejected_triplets
  ))
 }
}

# ── Allele calling (kernel density peak detection) ───────────────────────────

call_alleles <- function(repeat_lengths, frequencies, bw = 1.5) {
 
 total <- sum(frequencies)
 if (length(repeat_lengths) < 3 || total < 20) {
  return(data.frame(Allele = NA, Peak_Repeat = NA, Proportion = NA,
                    Call = "insufficient_data", stringsAsFactors = FALSE))
 }
 
 obs <- rep(repeat_lengths, times = frequencies)
 d   <- density(obs, bw = bw, n = 1024)
 
 # Detect peaks (local maxima)
 dy    <- diff(d$y)
 peaks <- which(dy[-length(dy)] > 0 & dy[-1] < 0) + 1L
 peak_x <- d$x[peaks]
 peak_y <- d$y[peaks]
 
 # Keep peaks >= 10% of tallest
 keep   <- peak_y >= 0.10 * max(peak_y)
 peak_x <- peak_x[keep]
 peak_y <- peak_y[keep]
 
 if (length(peak_x) == 0L) {
  return(data.frame(Allele = NA, Peak_Repeat = NA, Proportion = NA,
                    Call = "no_peaks", stringsAsFactors = FALSE))
 }
 
 # Top 2 by height
 ord    <- order(peak_y, decreasing = TRUE)
 peak_x <- round(peak_x[ord][seq_len(min(2, length(ord)))])
 peak_x <- sort(unique(peak_x))
 
 # Assign observations to nearest peak
 assignments <- vapply(obs, function(o) peak_x[which.min(abs(o - peak_x))],
                       numeric(1))
 props <- table(assignments) / length(obs)
 
 call <- if (length(peak_x) == 1L) "homozygous" else "heterozygous"
 
 data.frame(
  Allele      = seq_along(peak_x),
  Peak_Repeat = peak_x,
  Proportion  = as.numeric(props),
  Call        = call,
  stringsAsFactors = FALSE
 )
}

# ── Streaming FASTQ reader (ShortRead backend) ──────────────────────────────

fastq_streamer <- function(filepath, chunk_size = 50000L, q_threshold = 20) {
 
 sr_stream <- FastqStreamer(filepath, n = chunk_size)
 
 list(
  yield = function() {
   chunk <- yield(sr_stream)
   if (length(chunk) == 0L) return(NULL)
   
   n_total <- length(chunk)
   
   qual_matrix <- as(quality(chunk), "matrix")
   mean_quals  <- rowMeans(qual_matrix, na.rm = TRUE)
   keep        <- mean_quals >= q_threshold
   n_passed    <- sum(keep)
   
   if (n_passed == 0L) {
    return(list(
     data     = data.frame(id = character(0), seq = character(0),
                           qual = character(0), stringsAsFactors = FALSE),
     n_total  = n_total,
     n_passed = 0L
    ))
   }
   
   chunk <- chunk[keep]
   seqs  <- as.character(sread(chunk))
   ids   <- as.character(id(chunk))
   quals <- as.character(chunk@quality@quality)
   
   list(
    data     = data.frame(id = ids, seq = seqs, qual = quals,
                          stringsAsFactors = FALSE),
    n_total  = n_total,
    n_passed = n_passed
   )
  },
  close = function() close(sr_stream)
 )
}

# ── Main processor ───────────────────────────────────────────────────────────

process_fastq_file <- function(fastq_path,
                               sample_id         = NULL,
                               q_threshold       = 20,
                               min_repeat_length = 0L,
                               max_type_a        = 5L,
                               max_type_b        = 1L,
                               type_a_rate       = 0.05,
                               chunk_size        = 50000L,
                               allele_bw         = 1.5,
                               output_dir        = ".") {
 
 # ── Sample ID ──────────────────────────────────────────────────────────────
 if (is.null(sample_id)) {
  sample_id <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(fastq_path))
  sample_id <- sub("\\.gz$", "", sample_id)
 }
 
 # ── Streaming reader ───────────────────────────────────────────────────────
 stream <- fastq_streamer(fastq_path, chunk_size, q_threshold)
 on.exit(stream$close())
 
 total_reads     <- 0L
 reads_passing   <- 0L
 all_lengths     <- integer(0)
 
 # Diagnostic counters
 n_bridged_reads       <- 0L
 censored_reads        <- 0L
 threshold_rejected    <- 0L
 rejected_triplet_sizes <- integer(0)   # raw match sizes of rejected reads
 
 chunk_num <- 0L
 repeat {
  result_chunk <- stream$yield()
  if (is.null(result_chunk)) break
  
  chunk_num   <- chunk_num + 1L
  n           <- result_chunk$n_total
  n_passed    <- result_chunk$n_passed
  total_reads <- total_reads + n
  reads_passing <- reads_passing + n_passed
  chunk       <- result_chunk$data
  
  if (nrow(chunk) == 0L) {
   message(sprintf("  chunk %d: %d reads, 0 passed Q filter", chunk_num, n))
   next
  }
  
  for (i in seq_len(nrow(chunk))) {
   seq_i  <- chunk$seq[i]
   
   # ── CAG extraction — both strands ──────────────────────────────────
   result <- extract_cag_both_strands(seq_i, max_type_a, max_type_b,
                                      type_a_rate)
   
   if (!is.na(result$length)) {
    all_lengths <- c(all_lengths, result$length)
    
    if (grepl("N", seq_i, fixed = TRUE)) {
     n_bridged_reads <- n_bridged_reads + 1L
    }
    if (isTRUE(result$hits_end)) {
     censored_reads <- censored_reads + 1L
    }
   } else if (!is.na(result$rejected_triplets)) {
    # Read had a regex match but it was rejected by thresholds
    threshold_rejected     <- threshold_rejected + 1L
    rejected_triplet_sizes <- c(rejected_triplet_sizes,
                                result$rejected_triplets)
   }
  }
  
  message(sprintf("  chunk %d: %d reads processed (%d passed Q filter)",
                  chunk_num, n, n_passed))
 }
 
 reads_failed <- total_reads - reads_passing
 
 message(sprintf(
  "\n[%s] Total reads: %d | Passed Q>=%d: %d | Removed: %d",
  sample_id, total_reads, q_threshold, reads_passing, reads_failed))
 
 if (reads_passing == 0L) {
  warning(sprintf("No reads passed quality filter for %s.", sample_id))
  return(invisible(NULL))
 }
 
 reads_with_cag <- length(all_lengths)
 message(sprintf(
  "[%s] Reads with CAG match: %d / %d (%.1f%%)",
  sample_id, reads_with_cag, reads_passing,
  100 * reads_with_cag / reads_passing))
 
 message(sprintf(
  "[%s] N-bridged: %d | End-of-read censored: %d",
  sample_id, n_bridged_reads, censored_reads))
 
 # ── Threshold rejection diagnostics ────────────────────────────────────────
 message(sprintf(
  "[%s] Threshold-rejected reads: %d (type_a_rate=%.2f, max_type_a=%d, max_type_b=%d)",
  sample_id, threshold_rejected, type_a_rate, max_type_a, max_type_b))
 
 if (threshold_rejected > 0L) {
  message(sprintf(
   "[%s]   Rejected match sizes (triplets): min=%d, median=%d, max=%d",
   sample_id,
   min(rejected_triplet_sizes),
   as.integer(median(rejected_triplet_sizes)),
   max(rejected_triplet_sizes)))
 }
 
 if (reads_with_cag == 0L) {
  warning(sprintf("No CAG structures detected in %s.", sample_id))
  return(invisible(NULL))
 }
 
 # ── Build frequency table ──────────────────────────────────────────────────
 all_lengths <- all_lengths[all_lengths > 0L]
 freq_tbl    <- as.data.frame(table(all_lengths), stringsAsFactors = FALSE)
 names(freq_tbl) <- c("Repeat_Length", "Frequency")
 freq_tbl$Repeat_Length <- as.integer(freq_tbl$Repeat_Length)
 full_freq <- freq_tbl   # unfiltered copy
 
 # ── Apply user-configurable floor ──────────────────────────────────────────
 freq_tbl <- freq_tbl[freq_tbl$Repeat_Length >= min_repeat_length, ]
 
 if (nrow(freq_tbl) == 0L) {
  warning(sprintf("No repeats >= %d after correction in %s.",
                  min_repeat_length, sample_id))
  return(invisible(NULL))
 }
 
 # ── Mode (fixed: uses Frequency column) ────────────────────────────────────
 mode_repeat <- freq_tbl$Repeat_Length[which.max(freq_tbl$Frequency)]
 
 # ── Threshold metrics ──────────────────────────────────────────────────────
 reads_over_35  <- sum(freq_tbl$Frequency[freq_tbl$Repeat_Length >= 36])
 reads_over_110 <- sum(freq_tbl$Frequency[freq_tbl$Repeat_Length >  110])
 
 pct_over_110 <- if (reads_over_35 > 0) {
  (reads_over_110 / reads_over_35) * 100
 } else {
  NA_real_
 }
 
 # ── Allele calling ────────────────────────────────────────────────────────
 allele_df <- call_alleles(freq_tbl$Repeat_Length,
                           freq_tbl$Frequency,
                           bw = allele_bw)
 
 # ── Write outputs ──────────────────────────────────────────────────────────
 dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
 base_name <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(fastq_path))
 base_name <- sub("\\.gz$", "", base_name)
 
 sample_text <- paste0("CAG LOGGER v", VERSION, " : ", sample_id)
 
 # Metrics CSV
 metrics_df <- data.frame(
  Sample_ID                  = sample_text,
  Total_Reads                = total_reads,
  Reads_Passed_QFilter       = reads_passing,
  Reads_Failed_QFilter       = reads_failed,
  Reads_With_CAG             = reads_with_cag,
  Q_Threshold                = q_threshold,
  N_Bridged_Reads            = n_bridged_reads,
  EndOfRead_Censored_Reads   = censored_reads,
  Threshold_Rejected_Reads   = threshold_rejected,
  Type_A_Rate                = type_a_rate,
  Max_Type_A_Floor           = max_type_a,
  Max_Type_B                 = max_type_b,
  Min_Repeat_Length          = min_repeat_length,
  Mode_Repeat_Size           = mode_repeat,
  Reads_Over_35              = reads_over_35,
  Reads_Over_110             = reads_over_110,
  Percent_Reads_Over_110     = ifelse(is.na(pct_over_110), NA, round(pct_over_110, 2)),
  Allele_Call                = allele_df$Call[1],
  Allele_1_Peak              = allele_df$Peak_Repeat[1],
  Allele_1_Proportion        = allele_df$Proportion[1],
  Allele_2_Peak              = ifelse(nrow(allele_df) >= 2, allele_df$Peak_Repeat[2], NA),
  Allele_2_Proportion        = ifelse(nrow(allele_df) >= 2, allele_df$Proportion[2], NA),
  stringsAsFactors           = FALSE
 )
 
 metrics_csv   <- file.path(output_dir,
                            paste0(base_name, "_CAGLogger.v", VERSION, "_RepeatMetrics.csv"))
 freq_csv      <- file.path(output_dir,
                            paste0(base_name, "_CAGLogger.v", VERSION, "_FrequencyOfRepeats.csv"))
 full_freq_csv <- file.path(output_dir,
                            paste0(base_name, "_CAGLogger.v", VERSION, "_FullDistribution.csv"))
 allele_csv    <- file.path(output_dir,
                            paste0(base_name, "_CAGLogger.v", VERSION, "_AlleleCalls.csv"))
 
 write.csv(metrics_df, metrics_csv,   row.names = FALSE)
 write.csv(freq_tbl,   freq_csv,      row.names = FALSE)
 write.csv(full_freq,  full_freq_csv, row.names = FALSE)
 write.csv(allele_df,  allele_csv,    row.names = FALSE)
 
 # ── Write rejected match sizes if any ──────────────────────────────────────
 if (threshold_rejected > 0L) {
  rej_df <- data.frame(
   Rejected_Match_Triplets = rejected_triplet_sizes,
   stringsAsFactors = FALSE
  )
  rej_csv <- file.path(output_dir,
                       paste0(base_name, "_CAGLogger.v", VERSION, "_ThresholdRejected.csv"))
  write.csv(rej_df, rej_csv, row.names = FALSE)
  message(sprintf("[%s] Rejected match sizes written to: %s",
                  sample_id, rej_csv))
 }
 
 message(sprintf(
  "\n[%s] Mode: %d | Allele call: %s | Peaks: %s",
  sample_id, mode_repeat, allele_df$Call[1],
  paste(allele_df$Peak_Repeat, collapse = ", ")))
 
 if (!is.na(pct_over_110)) {
  message(sprintf(
   "[%s] Reads >35: %d | Reads >110: %d (%.1f%%)",
   sample_id, reads_over_35, reads_over_110, pct_over_110))
 } else {
  message(sprintf("[%s] Reads >35: %d | Reads >110: 0", sample_id, reads_over_35))
 }
 
 message(sprintf("[%s] Output: %s/", sample_id, output_dir))
 
 invisible(list(
  sample_id              = sample_id,
  metrics                = metrics_df,
  freq_table             = freq_tbl,
  full_freq              = full_freq,
  alleles                = allele_df,
  n_bridged_reads        = n_bridged_reads,
  censored_reads         = censored_reads,
  threshold_rejected     = threshold_rejected,
  rejected_triplet_sizes = rejected_triplet_sizes
 ))
}

# ══════════════════════════════════════════════════════════════════════════════
# RUN
# ══════════════════════════════════════════════════════════════════════════════

# Usage:
#   Rscript CAG_LOGGER_v2.1c.R  <input.fastq.gz>  [output_dir]  [type_a_rate]
#
# type_a_rate (optional, default 0.05): per-triplet error tolerance for the
# adaptive Type A threshold.  The effective threshold is:
#   max(max_type_a, floor(n_triplets * type_a_rate))

if (!interactive()) {
 
 args <- commandArgs(trailingOnly = TRUE)
 
 fastq_path  <- if (length(args) >= 1) args[1] else stop("Usage: Rscript CAG_LOGGER_v2.1c.R <input.fastq[.gz]> [output_dir] [type_a_rate]")
 output_dir  <- if (length(args) >= 2) args[2] else "."
 type_a_rate <- if (length(args) >= 3) as.numeric(args[3]) else 0.05
 
 message(strrep("=", 54))
 message("  CAG LOGGER v", VERSION)
 message(strrep("=", 54))
 message("Input:       ", fastq_path)
 message("Output:      ", output_dir)
 message("type_a_rate: ", type_a_rate)
 message("")
 
 t0 <- proc.time()
 
 result <- process_fastq_file(
  fastq_path        = fastq_path,
  sample_id         = NULL,
  q_threshold       = 20,
  min_repeat_length = 0L,
  max_type_a        = 5L,
  max_type_b        = 1L,
  type_a_rate       = type_a_rate,
  chunk_size        = 50000L,
  allele_bw         = 1.5,
  output_dir        = output_dir
 )
 
 elapsed <- (proc.time() - t0)["elapsed"]
 message(sprintf("\nDone in %.1f seconds.", elapsed))
 
 if (!is.null(result)) {
  message("\n── Top 15 Repeat Lengths ──")
  top <- head(result$freq_table[order(-result$freq_table$Frequency), ], 15)
  print(top, row.names = FALSE)
  
  if (result$threshold_rejected > 0L) {
   message(sprintf(
    "\n── %d reads were rejected by error thresholds ──",
    result$threshold_rejected))
   message("   See _ThresholdRejected.csv for raw match sizes.")
   message("   If many rejected matches are long, consider increasing type_a_rate.")
  }
 }
}

# ══════════════════════════════════════════════════════════════════════════════
# BATCH PROCESSING EXAMPLE
# ══════════════════════════════════════════════════════════════════════════════
#
# source("CAG_LOGGER_v2.1c.R")
#
# fastq_files <- list.files(pattern = "\\.(fastq|fq)(\\.gz)?$")
#
library(tictoc)
tic()
results <- list()
for (file in fastq_files) {
  results[[file]] <- process_fastq_file(
    fastq_path        = file,
    sample_id         = NULL,
    q_threshold       = 20,
    min_repeat_length = 10L,
    max_type_a        = 5L,
    max_type_b        = 1L,
    type_a_rate       = 0.05,          # 5% per-triplet error tolerance
    output_dir        = "~/Desktop/CAGlogger_Github/cag2.1results/2.1d_combinedBandC"
  )
}
toc()