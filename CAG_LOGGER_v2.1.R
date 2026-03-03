#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# CAG LOGGER v2.1 — R implementation
# Alignment-free estimation of CAG repeats from targeted amplicon MiSeq data
#
# Changes from v1.2 (peer review response + artifact fixes):
#
#   Fix A — N-aware regex: pattern includes N as valid base at mismatch
#           positions, bridging cycle-specific N-calls (Q=2 dropouts) instead
#           of splitting the match.
#
#   Fix B — Adaptive flank correction: the blanket -3 is replaced by a
#           per-read check.  Correction applied ONLY when the regex match
#           demonstrably absorbed the downstream CAA·CAG·CCG flank.
#
#   Fix C — N-triplet accounting: N-containing triplets tracked separately
#           with a configurable ceiling (max_n_triplets, default 2).
#
#   Fix D — Optional per-base quality masking: bases below a quality floor
#           are replaced with N before the regex runs.  Disabled by default.
#
#   Fix E — End-of-read correction: when a match reaches within 3 bp of read
#           end, trailing non-pure-CAG triplets are subtracted as terminal
#           quality junk, and the read is flagged as right-censored.
#
#   From v2.0 peer review:
#     - Reverse-complement handling (both strands searched)
#     - Fixed mode calculation (uses Frequency column, not table-of-table)
#     - Configurable sample_id, min_repeat_length, error thresholds
#     - Per-read single longest match (no multi-match inflation)
#     - Incomplete triplet handling (trailing <3 bases dropped)
#     - Streaming FASTQ reader for memory efficiency
#     - Allele calling via kernel density peak detection
#     - Version strings unified
#
# Dependencies: base R only (no Bioconductor required).
#   Optional: ShortRead + Biostrings (auto-detected at runtime)
#     — If present: uses FastqStreamer for ~3-5× faster FASTQ I/O and
#       vectorised PhredQuality→matrix conversion for quality filtering.
#     — If absent:  graceful fallback to base R gzfile/readLines.
# ══════════════════════════════════════════════════════════════════════════════

VERSION <- "2.1"

# ── Reverse complement (both strands searched) ───────────────────────────────

rc <- function(seq) {
 chartr("ACGTacgtN", "TGCAtgcaN",
        paste(rev(strsplit(seq, "")[[1]]), collapse = ""))
}

# # ── Quality helpers ──────────────────────────────────────────────────────────
# 
# mean_phred <- function(qual_string) {
#  mean(as.integer(charToRaw(qual_string)) - 33L)
# }

#' Vectorised mean-quality filter for a character vector of quality strings.
#' When ShortRead is available, converts via PhredQuality → integer matrix
#' (compiled C, much faster for large chunks).  Otherwise uses base R.
# mean_phred_vec <- function(qual_strings) {
#  if (HAS_SHORTREAD && length(qual_strings) > 100L) {
#   qual_matrix <- as(Biostrings::PhredQuality(qual_strings), "matrix")
#   rowMeans(qual_matrix, na.rm = TRUE)
#  } else {
#   vapply(qual_strings, mean_phred, numeric(1))
#  }
# }

#' Fix D: Replace individual bases with quality < floor with N.
#' Converts low-confidence miscalls into explicit Ns that the N-aware
#' regex can bridge.
# mask_low_quality_bases <- function(seq, qual, per_base_q_floor = 10L) {
#  if (per_base_q_floor <= 0L) return(seq)
#  seq_chars  <- strsplit(seq,  "")[[1]]
#  qual_ints  <- as.integer(charToRaw(qual)) - 33L
#  seq_chars[qual_ints < per_base_q_floor] <- "N"
#  paste(seq_chars, collapse = "")
# }

# ── Boundary processing ─────────────────────────────────────────────────────

apply_boundary <- function(seq,
                           boundary_mode = "none",
                           left_anchor   = NULL,
                           right_anchor  = NULL) {
 
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
 seq
}

# ── Flank absorption detector (Fix B) ───────────────────────────────────────
#
# The known HTT exon 1 structure downstream of the pure CAG tract is:
#   ...CAG CAG CAA CAG CCG CCA CCG CCG CCG ...
#                ^^^^^^^^^^^
#                These 3 triplets are absorbed by the greedy regex.
#
# Strategy: inspect the last 5 triplets of the match.  If any is CCG-like
# AND any is CAA-like, infer full flank absorption (subtract 3).  If only
# CCG-like is present, infer partial absorption (subtract 2).

is_ccg_like <- function(triplet) {
 # TRUE if triplet is CCG or a 1-substitution/N variant, but NOT CAG itself
 if (triplet == "CAG") return(FALSE)
 chars <- strsplit(triplet, "")[[1]]
 ref   <- c("C", "C", "G")
 mismatches <- sum(chars != ref & chars != "N")
 mismatches <= 1L
}

detect_flank_absorption <- function(match_str) {
 usable <- (nchar(match_str) %/% 3L) * 3L
 if (usable < 9L) return(0L)
 
 n_tail     <- min(5L, usable %/% 3L)
 tail_start <- usable - (n_tail * 3L) + 1L
 tail_triplets <- vapply(
  seq(tail_start, usable, by = 3L),
  function(i) substr(match_str, i, i + 2L),
  character(1)
 )
 
 has_ccg <- any(vapply(tail_triplets, is_ccg_like, logical(1)))
 has_caa <- any(tail_triplets %in% c("CAA", "CAN", "NAA"))
 
 if (has_ccg && has_caa) return(3L)
 if (has_ccg)            return(2L)
 0L
}

# ── CAG extraction (N-aware, adaptive correction) ───────────────────────────
#
# Fix A: regex includes N at each mismatch position.

CAG_PATTERN <- "(CAG|CA[ATCGN]|C[ATCGN]G|[ATCGN]AG)+"

#' Extract the longest valid CAG match from a single sequence.
#' Returns a named list: length (integer or NA), correction, hits_end.
extract_cag_from_read <- function(seq,
                                  max_type_a     = 5L,
                                  max_type_b     = 1L,
                                  max_n_triplets = 2L,
                                  boundary_mode  = "none") {
 
 matches <- gregexpr(CAG_PATTERN, seq, perl = TRUE)
 raw_hits <- regmatches(seq, matches)[[1]]
 
 if (length(raw_hits) == 0L) {
  return(list(length = NA_integer_, correction = 0L, hits_end = FALSE))
 }
 
 # We also need match start positions to detect end-of-read
 starts  <- as.integer(matches[[1]])
 lengths <- attr(matches[[1]], "match.length")
 read_len <- nchar(seq)
 
 best_length     <- NA_integer_
 best_correction <- 0L
 best_hits_end   <- FALSE
 
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
  
  # ── Error counting ─────────────────────────────────────────────────
  type_a   <- 0L
  type_b   <- 0L
  n_count  <- 0L
  prev_err <- FALSE
  
  for (t in triplets) {
   has_n  <- grepl("N", t, fixed = TRUE)
   is_err <- (t != "CAG")
   
   if (has_n)  n_count <- n_count + 1L
   if (is_err) {
    type_a <- type_a + 1L
    if (prev_err) type_b <- type_b + 1L
   }
   prev_err <- is_err
  }
  
  # Fix C: reject if too many N-triplets
  if (n_count > max_n_triplets) next
  if (type_a > max_type_a || type_b > max_type_b) next
  
  n_triplets <- length(triplets)
  
  # ── Does the match reach the end of the read? ──────────────────────
  match_end_pos <- starts[idx] + usable_len - 1L
  hits_read_end <- (read_len - match_end_pos) <= 3L
  
  # ── Adaptive correction ────────────────────────────────────────────
  if (boundary_mode != "none") {
   correction <- 0L
  } else if (hits_read_end) {
   # Fix E: end-of-read truncation.
   # Count trailing non-pure-CAG triplets and subtract them.
   tail_junk <- 0L
   for (t in rev(triplets)) {
    if (t != "CAG") { tail_junk <- tail_junk + 1L }
    else            { break }
   }
   correction <- tail_junk
  } else {
   correction <- detect_flank_absorption(match_str)
  }
  
  corrected <- n_triplets - correction
  
  if (is.na(best_length) || corrected > best_length) {
   best_length     <- corrected
   best_correction <- correction
   best_hits_end   <- hits_read_end
  }
 }
 
 list(length = best_length, correction = best_correction, hits_end = best_hits_end)
}

#' Search both orientations, return best result.
extract_cag_both_strands <- function(seq,
                                     max_type_a     = 5L,
                                     max_type_b     = 1L,
                                     max_n_triplets = 2L,
                                     boundary_mode  = "none") {
 
 fwd <- extract_cag_from_read(seq,     max_type_a, max_type_b,
                              max_n_triplets, boundary_mode)
 rev <- extract_cag_from_read(rc(seq), max_type_a, max_type_b,
                              max_n_triplets, boundary_mode)
 
 fwd_len <- ifelse(is.na(fwd$length), 0L, fwd$length)
 rev_len <- ifelse(is.na(rev$length), 0L, rev$length)
 
 if (fwd_len >= rev_len) {
  return(list(
   length     = if (fwd_len > 0L) fwd_len else NA_integer_,
   correction = fwd$correction,
   hits_end   = fwd$hits_end
  ))
 } else {
  return(list(
   length     = if (rev_len > 0L) rev_len else NA_integer_,
   correction = rev$correction,
   hits_end   = rev$hits_end
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

# ── Streaming FASTQ reader ───────────────────────────────────────────────────
#
# Two backends:
#
#   1. ShortRead::FastqStreamer  (preferred when ShortRead is installed)
#      — uses compiled C for FASTQ parsing / gz decompression; ~3-5× faster
#      — quality scores come as PhredQuality objects → integer matrix in bulk
#
#   2. Base-R readLines fallback (no dependencies)
#
# The factory function fastq_streamer() auto-detects which is available and
# returns an object with the same $yield() / $close() interface regardless.

HAS_SHORTREAD <- requireNamespace("ShortRead", quietly = TRUE)

if (HAS_SHORTREAD) {
 message("[CAG LOGGER] ShortRead detected — using FastqStreamer backend")
} else {
 message("[CAG LOGGER] ShortRead not found — using base-R FASTQ reader")
}

#' @param filepath  Path to .fastq or .fastq.gz
#' @param chunk_size  Number of reads per chunk
#' @return List with \code{$yield()} returning a data.frame(id, seq, qual)
#'         and \code{$close()}.
fastq_streamer <- function(filepath, chunk_size = 50000L) {
 
 if (HAS_SHORTREAD) {
  return(.fastq_streamer_shortread(filepath, chunk_size))
 } else {
  return(.fastq_streamer_base(filepath, chunk_size))
 }
}

# ── ShortRead backend ────────────────────────────────────────────────────────

.fastq_streamer_shortread <- function(filepath, chunk_size) {
 
 sr_stream <- ShortRead::FastqStreamer(filepath, n = chunk_size)
 
 list(
  yield = function() {
   chunk <- ShortRead::yield(sr_stream)
   if (length(chunk) == 0L) return(NULL)
   
   # Extract character vectors for sequences and quality strings
   seqs  <- as.character(ShortRead::sread(chunk))
   ids   <- as.character(ShortRead::id(chunk))
   
   # Quality: access the quality slot directly and convert to character.
   # quality() generic lives in Biostrings, not ShortRead, so we use
   # the @ slot accessor to avoid namespace issues.
   quals <- as.character(chunk@quality@quality)
   
   data.frame(id = ids, seq = seqs, qual = quals,
              stringsAsFactors = FALSE)
  },
  close = function() close(sr_stream)
 )
}

# ── Base-R fallback backend ──────────────────────────────────────────────────

.fastq_streamer_base <- function(filepath, chunk_size) {
 
 if (grepl("\\.gz$", filepath)) {
  con <- gzfile(filepath, open = "rt")
 } else {
  con <- file(filepath, open = "rt")
 }
 exhausted <- FALSE
 
 list(
  yield = function() {
   if (exhausted) return(NULL)
   ids   <- character(0)
   seqs  <- character(0)
   quals <- character(0)
   count <- 0L
   
   while (count < chunk_size) {
    line1 <- readLines(con, n = 1L)
    if (length(line1) == 0L) { exhausted <<- TRUE; break }
    line2 <- readLines(con, n = 1L)
    readLines(con, n = 1L)          # "+"
    line4 <- readLines(con, n = 1L)
    if (length(line4) == 0L) { exhausted <<- TRUE; break }
    
    count <- count + 1L
    ids   <- c(ids,   line1)
    seqs  <- c(seqs,  line2)
    quals <- c(quals, line4)
   }
   
   if (count == 0L) return(NULL)
   data.frame(id = ids, seq = seqs, qual = quals,
              stringsAsFactors = FALSE)
  },
  close = function() close(con)
 )
}

# ── Main processor ───────────────────────────────────────────────────────────

process_fastq_file <- function(fastq_path,
                               sample_id         = NULL,
                               q_threshold       = 20,
                               boundary_mode     = c("none", "trim", "mask"),
                               left_anchor       = NULL,
                               right_anchor      = NULL,
                               min_repeat_length = 0L,
                               max_type_a        = 5L,
                               max_type_b        = 1L,
                               max_n_triplets    = 2L,
                               per_base_q_floor  = 0L,
                               chunk_size        = 50000L,
                               allele_bw         = 1.5,
                               output_dir        = ".") {
 
 boundary_mode <- match.arg(boundary_mode)
 
 # ── Sample ID ──────────────────────────────────────────────────────────────
 if (is.null(sample_id)) {
  sample_id <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(fastq_path))
  sample_id <- sub("\\.gz$", "", sample_id)
 }
 
 # ── Streaming reader ───────────────────────────────────────────────────────
 stream <- fastq_streamer(fastq_path, chunk_size)
 on.exit(stream$close())
 
 total_reads     <- 0L
 reads_passing   <- 0L
 all_lengths     <- integer(0)
 
 # Diagnostic counters
 correction_counts <- c("0" = 0L, "2" = 0L, "3" = 0L)
 n_bridged_reads   <- 0L
 censored_reads    <- 0L
 
 chunk_num <- 0L
 repeat {
  chunk <- stream$yield()
  if (is.null(chunk)) break
  
  chunk_num   <- chunk_num + 1L
  n           <- nrow(chunk)
  total_reads <- total_reads + n
  
  # ── Mean quality filter (vectorised when ShortRead available) ───────
  mean_quals <- mean_phred_vec(chunk$qual)
  keep       <- mean_quals >= q_threshold
  chunk      <- chunk[keep, , drop = FALSE]
  reads_passing <- reads_passing + sum(keep)
  
  if (nrow(chunk) == 0L) {
   message(sprintf("  chunk %d: %d reads, 0 passed Q filter", chunk_num, n))
   next
  }
  
  for (i in seq_len(nrow(chunk))) {
   seq_i  <- chunk$seq[i]
   qual_i <- chunk$qual[i]
   
   # ── Fix D: per-base quality masking (optional) ─────────────────────
   if (per_base_q_floor > 0L) {
    seq_i <- mask_low_quality_bases(seq_i, qual_i, per_base_q_floor)
   }
   
   # ── Boundary processing ────────────────────────────────────────────
   if (boundary_mode != "none" && !is.null(left_anchor) && !is.null(right_anchor)) {
    seq_i <- apply_boundary(seq_i, boundary_mode, left_anchor, right_anchor)
   }
   
   # ── CAG extraction — both strands ──────────────────────────────────
   result <- extract_cag_both_strands(seq_i, max_type_a, max_type_b,
                                      max_n_triplets, boundary_mode)
   
   if (!is.na(result$length)) {
    all_lengths <- c(all_lengths, result$length)
    
    corr_key <- as.character(result$correction)
    if (corr_key %in% names(correction_counts)) {
     correction_counts[corr_key] <- correction_counts[corr_key] + 1L
    } else {
     # Other correction values (e.g., tail_junk counts from Fix E)
     if ("0" %in% names(correction_counts)) {
      # Track non-standard corrections under a general bucket
      correction_counts["0"] <- correction_counts["0"]
     }
     # Or add them
     if (is.na(correction_counts[corr_key])) {
      correction_counts[corr_key] <- 1L
     } else {
      correction_counts[corr_key] <- correction_counts[corr_key] + 1L
     }
    }
    
    if (grepl("N", seq_i, fixed = TRUE)) {
     n_bridged_reads <- n_bridged_reads + 1L
    }
    if (isTRUE(result$hits_end)) {
     censored_reads <- censored_reads + 1L
    }
   }
  }
  
  message(sprintf("  chunk %d: %d reads processed (%d passed Q filter)",
                  chunk_num, n, sum(keep)))
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
 
 # Group correction counts for reporting: 0, 2, 3, and "other"
 corr_0 <- sum(correction_counts[names(correction_counts) == "0"], na.rm = TRUE)
 corr_2 <- sum(correction_counts[names(correction_counts) == "2"], na.rm = TRUE)
 corr_3 <- sum(correction_counts[names(correction_counts) == "3"], na.rm = TRUE)
 
 message(sprintf(
  "[%s] Adaptive correction: corr=0: %d, corr=2: %d, corr=3: %d | N-bridged: %d | End-of-read censored: %d",
  sample_id, corr_0, corr_2, corr_3, n_bridged_reads, censored_reads))
 
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
  Boundary_Mode              = boundary_mode,
  Corrections_0              = corr_0,
  Corrections_2              = corr_2,
  Corrections_3              = corr_3,
  N_Bridged_Reads            = n_bridged_reads,
  EndOfRead_Censored_Reads   = censored_reads,
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
  sample_id         = sample_id,
  metrics           = metrics_df,
  freq_table        = freq_tbl,
  full_freq         = full_freq,
  alleles           = allele_df,
  correction_counts = correction_counts,
  n_bridged_reads   = n_bridged_reads,
  censored_reads    = censored_reads
 ))
}

# ══════════════════════════════════════════════════════════════════════════════
# RUN
# ══════════════════════════════════════════════════════════════════════════════

# This section runs when the script is called from the command line:
#   Rscript CAG_LOGGER_v2.1.R  <input.fastq.gz>  [output_dir]
#
# To use interactively or in a pipeline, source() this file and call
# process_fastq_file() directly.

if (!interactive()) {
 
 args <- commandArgs(trailingOnly = TRUE)
 
 fastq_path <- if (length(args) >= 1) args[1] else stop("Usage: Rscript CAG_LOGGER_v2.1.R <input.fastq[.gz]> [output_dir]")
 output_dir <- if (length(args) >= 2) args[2] else "."
 
 message(strrep("=", 54))
 message("  CAG LOGGER v", VERSION)
 message(strrep("=", 54))
 message("Input:  ", fastq_path)
 message("Output: ", output_dir)
 message("")
 
 t0 <- proc.time()
 
 result <- process_fastq_file(
  fastq_path        = fastq_path,
  sample_id         = NULL,
  q_threshold       = 20,
  boundary_mode     = "none",
  left_anchor       = NULL,
  right_anchor      = NULL,
  min_repeat_length = 0L,
  max_type_a        = 5L,
  max_type_b        = 1L,
  max_n_triplets    = 2L,
  per_base_q_floor  = 0L,      # set to 10 to enable Fix D
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
 }
}

# ══════════════════════════════════════════════════════════════════════════════
# BATCH PROCESSING EXAMPLE
# ══════════════════════════════════════════════════════════════════════════════
#
# source("CAG_LOGGER_v2.1.R")
#
# fastq_files <- list.files(pattern = "\\.(fastq|fq)(\\.gz)?$")
#fastq_files <- list.files(path = "../testfqs", pattern = "\\.(fastq|fq)(\\.gz)?$", full.names = TRUE)
#
library(tictoc)
tic()
results <- list()
for (file in fastq_files) {
  results[[file]] <- process_fastq_file(
    fastq_path        = file,
    sample_id         = NULL,
    q_threshold       = 20,
    boundary_mode     = "none",       # "trim" or "mask" with anchors
    left_anchor       = NULL,         # e.g. "GCGACCCTGGAAAAGCTGA"
    right_anchor      = NULL,         # e.g. "GGCGGCTGAGGAAGCTGAG"
    min_repeat_length = 10,           # set to either OL; or 35 for v1.2 behaviour
    max_type_a        = 5L,
    max_type_b        = 1L,
    max_n_triplets    = 2L,
    per_base_q_floor  = 0L,           # set to 10 to enable Fix D
    chunk_size        = 50000L,
    allele_bw         = 1.5,
    output_dir        = "caglogger2.1_results/"
  )
}
toc()