#' Go_intoASV
#'
#' Compute within-taxon nucleotide diversity (π) from ASV-level 16S sequences
#' in a phyloseq object, optionally trimming noisy V3–V4 termini, clustering,
#' and bootstrapping confidence intervals (CIs). Supports both taxonomy-based
#' and similarity-based grouping.
#'
#' @author Heekuk Park <hp2523@cumc.columbia.edu>
#' Created on 2025-11-10
#'
#' @description
#' For each taxonomic group (or all ASVs), this function:
#' \itemize{
#'   \item extracts ASV sequences (\code{refseq(psIN)} or rownames),
#'   \item optionally groups ASVs by global sequence similarity (\code{global_similarity_cutoff}),
#'   \item clusters intra-group sequences at \code{taxonomy_cluster_cutoff} to remove near-identical ASVs,
#'   \item aligns sequences via \code{DECIPHER::AlignSeqs} or external \code{MAFFT},
#'   \item trims both ends by \code{trim_nt} nt to reduce V3–V4 edge noise,
#'   \item computes pairwise distances (\code{"simple"} = Hamming; or \code{"nucdiv"} = TN93 model),
#'   \item calculates per-sample nucleotide diversity (π), weighted by relative abundance
#'         or entropy, and optionally bootstrapped confidence intervals,
#'   \item saves π and ASV-count matrices, and merges them into \code{sample_data(psIN)}.
#' }
#'
#' @param psIN A \code{phyloseq} object containing ASV-level abundance data and
#'   optionally DNA sequences in \code{refseq(psIN)}.
#' @param project Character; project prefix used to create the output directory
#'   \verb{<project_YYMMDD>/table/pi_tab/}.
#' @param level Character or \code{NULL}; taxonomic rank to group ASVs before computing π
#'   (e.g., \code{"Genus"}, \code{"Family"}). Required in taxonomy mode and
#'   ignored if \code{global_similarity_cutoff} is set. Default \code{NULL}.
#' @param target Character; a specific taxon (e.g., \code{"Lactobacillus"}) or
#'   \code{"all"} to compute across all taxa. Ignored if \code{global_similarity_cutoff} is set.
#'   Default \code{"all"}.
#' @param method Character; one of \code{c("simple","nucdiv")}.
#'   \code{"simple"} uses per-position Hamming distances, while \code{"nucdiv"}
#'   applies substitution models via \code{ape::dist.dna()} (e.g., TN93).
#' @param aligner Character; alignment engine, either \code{"DECIPHER"} (R-based)
#'   or \code{"MAFFT"} (external, faster if installed). Default \code{"DECIPHER"}.
#' @param min_asv Integer; minimum number of ASVs required to compute π. Default \code{3}.
#' @param min_abund Numeric; minimum total abundance per sample to be included. Default \code{10}.
#'
#' @param taxonomy_cluster_cutoff Numeric; intra-group clustering threshold applied
#'   only within each taxon in taxonomy mode to remove highly similar ASVs
#'   (e.g., 99.5\% identity) before diversity computation. Default \code{0.995}.
#'   This acts as a denoising step to prevent redundant ASVs from inflating π.
#'
#' @param global_similarity_cutoff Numeric or \code{NULL}.
#'   If set (e.g., \code{0.97}), taxonomy is ignored and ASVs are globally
#'   clustered based on sequence similarity. Each resulting cluster is treated
#'   as a “species-like” group for π calculation. The output files will include
#'   the similarity tag (e.g., \code{_similarity_0.970_}) and a
#'   \code{cluster_map_similarity_<cutoff>_<date>.csv} with taxonomy annotation.
#'   When \code{global_similarity_cutoff} is active, \code{taxonomy_cluster_cutoff} is skipped.
#'
#' @param trim_nt Integer; number of nucleotides trimmed from both sequence ends
#'   to remove noisy termini of V3–V4 amplicons. Default \code{8}.
#' @param distance_gap Character; whether to \code{"exclude"} or \code{"include"}
#'   gap positions when computing nucleotide distances. Default \code{"exclude"}.
#' @param distance_model Character; substitution model for \code{method="nucdiv"}
#'   (e.g., \code{"raw"}, \code{"JC69"}). Default \code{"raw"}.
#' @param weighting Character; weighting scheme for π:
#'   \code{"abundance"} (relative abundance) or \code{"entropy"} (p·(1–p)).
#'   Default \code{"abundance"}.
#' @param compute_ci Logical; if \code{TRUE}, compute 95\% bootstrap confidence
#'   intervals using \code{n_boot} replicates. Default \code{FALSE}.
#' @param n_boot Integer; number of bootstrap replicates for confidence intervals.
#'   Default \code{200}.
#' @param seed Integer; random seed. Default \code{123}.
#' @param n_cores Integer; number of parallel cores for per-taxon computation. Default \code{4}.
#'
#' @details
#' \strong{Conceptual difference between cutoffs:}
#' \itemize{
#'   \item \code{global_similarity_cutoff} — defines global sequence-based grouping
#'         across all ASVs (taxonomy-free, analogous to ANI dereplication).
#'         Recommended for species-level clustering (e.g., 0.97).
#'   \item \code{taxonomy_cluster_cutoff} — defines local denoising within each taxon
#'         to merge nearly identical ASVs (e.g., 0.995).
#'         Used only when \code{global_similarity_cutoff = NULL}.
#' }
#'
#' π (nucleotide diversity) per sample is computed as:
#' \deqn{π = 2 \sum_{i<j} w_i w_j d_{ij}}
#' where \(w_i\) are normalized weights (abundance or entropy) and
#' \(d_{ij}\) is pairwise sequence distance.
#'
#' When \code{compute_ci = TRUE}, 95\% confidence intervals are obtained via
#' bootstrap resampling of ASVs within each sample.
#'
#' The function saves:
#' \itemize{
#'   \item \code{pi_matrix_<method>_<level-or-mode>_<target>_<date>.csv}
#'   \item \code{asv_count_matrix_<level-or-mode>_<target>_<date>.csv}
#'   \item optional \code{pi_ci_low_matrix_<method>_<level-or-mode>_<target>_<date>.csv},
#'         \code{pi_ci_high_matrix_<method>_<level-or-mode>_<target>_<date>.csv}
#'   \item text log: \code{pi_log_<method>_<level-or-mode>_<target>_<date>.txt}
#'   \item and if \code{global_similarity_cutoff} is set:
#'         \code{cluster_map_similarity_<cutoff>_<date>.csv}
#' }
#'
#' These are also merged into \code{sample_data(psIN)} with new columns
#' prefixed by \code{pi_} and \code{asvN_}.
#'
#' @return A \code{phyloseq} object identical to the input but with new columns
#'   added to \code{sample_data(psIN)} for each computed π and ASV count.
#'   Also writes CSVs and a log file under the project directory.
#'
#' @examples
#' \dontrun{
#' # taxonomy mode
#' ps_pi <- Go_intoASV(
#'   psIN = ps,
#'   project = "Gut16S",
#'   level = "Genus",
#'   target = "Lactobacillus",
#'   method = "nucdiv",
#'   aligner = "DECIPHER",
#'   taxonomy_cluster_cutoff = 0.995,
#'   compute_ci = TRUE,
#'   n_boot = 100
#' )
#'
#' # similarity mode (taxonomy ignored)
#' ps_sim <- Go_intoASV(
#'   psIN = ps,
#'   project = "Gut16S",
#'   global_similarity_cutoff = 0.97,
#'   method = "nucdiv"
#' )
#' }
#'
#' @importFrom phyloseq otu_table taxa_are_rows tax_table refseq sample_data
#' @importFrom DECIPHER AlignSeqs DistanceMatrix
#' @importFrom Biostrings DNAStringSet writeXStringSet readDNAStringSet subseq letterFrequency
#' @importFrom ape dist.dna as.DNAbin
#' @importFrom parallel mclapply
#' @importFrom stats quantile
#' @importFrom utils write.csv
#' @export


Go_intoASV <- function(
    psIN,
    project,
    level = NULL,
    target = "all",
    method = c("simple","nucdiv"),
    aligner = c("DECIPHER","MAFFT"),
    min_asv = 3,
    min_abund = 10,
    taxonomy_cluster_cutoff = 0.995,   # taxonomy-mode intra-taxon subclustering
    global_similarity_cutoff = NULL,   # if set: similarity-based mode (no intra-group subclustering)
    trim_nt = 8,
    distance_gap = c("exclude","include"),
    distance_model = c("raw","JC69"),  # kept for compatibility; nucdiv는 TN93 사용
    weighting = c("abundance","entropy"),
    compute_ci = FALSE,
    n_boot = 200,
    seed = 123,
    n_cores = 4,
    clustering_cutoff = NULL,
    similarity_cutoff = NULL
){

  start_time <- Sys.time()
  set.seed(seed)

  if (!is.null(clustering_cutoff)) {
    if (!is.null(taxonomy_cluster_cutoff) &&
        !isTRUE(all.equal(taxonomy_cluster_cutoff, clustering_cutoff))) {
      stop("Conflicting values supplied for 'taxonomy_cluster_cutoff' and deprecated 'clustering_cutoff'.")
    }
    warning("'clustering_cutoff' is deprecated; use 'taxonomy_cluster_cutoff' instead.",
            call. = FALSE)
    taxonomy_cluster_cutoff <- clustering_cutoff
  }
  if (!is.null(similarity_cutoff)) {
    if (!is.null(global_similarity_cutoff) &&
        !isTRUE(all.equal(global_similarity_cutoff, similarity_cutoff))) {
      stop("Conflicting values supplied for 'global_similarity_cutoff' and deprecated 'similarity_cutoff'.")
    }
    warning("'similarity_cutoff' is deprecated; use 'global_similarity_cutoff' instead.",
            call. = FALSE)
    global_similarity_cutoff <- similarity_cutoff
  }

  method         <- match.arg(method)
  aligner        <- match.arg(aligner)
  distance_gap   <- match.arg(distance_gap)
  distance_model <- match.arg(distance_model)
  weighting      <- match.arg(weighting)
  ###############################################
  # 0. Output directory structure
  ###############################################
  date_tag <- format(Sys.Date(), "%y%m%d")
  dir_base <- sprintf("%s_%s/intoASV/pi_tab", project, date_tag)
  dir.create(dir_base, recursive = TRUE, showWarnings = FALSE)

  # ---------- IO helpers ----------
  sanitize <- function(x) gsub("[^A-Za-z0-9._-]+","_", x)
  safe_target  <- if (identical(target,"all")) "all" else sanitize(trimws(target))
  level_label <- if (is.null(level)) "NULL" else as.character(level)

  # --- 파일명 태그 자동 정의 ---
  if (!is.null(global_similarity_cutoff)) {
    tag_label <- sprintf("similarity_%.3f", global_similarity_cutoff)
    level_tag <- "Similarity"
    target_tag <- tag_label
  } else {
    level_tag <- level
    target_tag <- safe_target
  }

  pi_file  <- sprintf("%s/pi_matrix_%s_%s_%s_%s.csv",
                      dir_base, method, level_tag, target_tag, date_tag)
  asv_file <- sprintf("%s/asv_count_matrix_%s_%s_%s.csv",
                      dir_base, level_tag, target_tag, date_tag)
  log_file <- sprintf("%s/pi_log_%s_%s_%s_%s.txt",
                      dir_base, method, level_tag, target_tag, date_tag)

  cat(sprintf("[Go_intoASV v13.4] %s | method=%s | level=%s | target=%s\n",
              Sys.time(), method, level_label, target),
      file = log_file, append = TRUE)

  # ---------- extract tables ----------
  otu_tab <- as(otu_table(psIN), "matrix")
  if (!taxa_are_rows(psIN)) otu_tab <- t(otu_tab)
  tax_tab <- as(tax_table(psIN), "matrix")

  # ----- 1) try refseq -----
  seqs_tmp <- NULL
  if (!is.null(refseq(psIN, errorIfNULL = FALSE))) {
    ref <- as.character(refseq(psIN))
    if (all(grepl("^[ACGTN]+$", ref))) {
      message("[Go_intoASV] Using refseq as DNA sequences.")
      seqs_tmp <- ref
    }
  }

  # ----- 2) try taxa_names -----
  if (is.null(seqs_tmp)) {
    tx <- taxa_names(psIN)
    if (all(grepl("^[ACGTN]+$", tx))) {
      message("[Go_intoASV] Using taxa_names as DNA sequences.")
      seqs_tmp <- tx
    }
  }

  # ----- 3) try rownames/colnames ONLY if DNA -----
  if (is.null(seqs_tmp)) {
    rn <- rownames(as(otu_table(psIN), "matrix"))
    if (all(grepl("^[ACGTN]+$", rn))) seqs_tmp <- rn
    cn <- colnames(as(otu_table(psIN), "matrix"))
    if (is.null(seqs_tmp) && all(grepl("^[ACGTN]+$", cn))) seqs_tmp <- cn
  }

  # ----- 4) stop if no DNA found -----
  if (is.null(seqs_tmp)) {
    stop("[Go_intoASV] ❌ No valid DNA sequences found in refseq / taxa_names / OTU names.")
  }

  # ----- 5) assign IDs -----
  names(seqs_tmp) <- taxa_names(psIN)
  seqs_all <- seqs_tmp

  message(sprintf("[Go_intoASV] Loaded %d DNA sequences", length(seqs_all)))

  # ---------- helper for clustering (average linkage) ----------
  cluster_from_dm <- function(dm, cutoff) {
    if (is.null(dm) || nrow(dm) < 2L) {
      return(setNames(rep(1L, nrow(dm)), rownames(dm)))
    }
    hc <- hclust(as.dist(dm), method = "average")
    cl <- cutree(hc, h = 1 - cutoff)
    names(cl) <- rownames(dm)
    return(cl)
  }

  # ---------- grouping: taxonomy vs similarity ----------
  if (!is.null(global_similarity_cutoff)) {
    message(sprintf(
      "\n[Go_intoASV] Similarity-based mode activated (cutoff = %.3f). 'level' and 'target' will be ignored.\n",
      global_similarity_cutoff
    ))
    cat(sprintf("[Go_intoASV] Similarity mode: cutoff=%.3f | Ignoring level & target\n",
                global_similarity_cutoff), file = log_file, append = TRUE)

    # 모든 ASV를 대상으로 전역 DistanceMatrix
    seqs_valid <- seqs_all[grepl("^[ACGTN]+$", seqs_all)]
    if (length(seqs_valid) < 2L) {
      cat("ERROR: Not enough valid sequences for similarity-based clustering.\n",
          file = log_file, append = TRUE)
      stop("Not enough valid sequences for similarity-based clustering.")
    }

    dna_all <- DNAStringSet(seqs_valid)
    dm_global <- suppressWarnings(
      DECIPHER::DistanceMatrix(dna_all, includeTerminalGaps = FALSE)
    )

    # ---- base R clustering (average linkage) ----
    cl_vec <- cluster_from_dm(dm_global, global_similarity_cutoff)

    # ---- taxonomy info 추가 ----
    tax_cols <- intersect(c("Phylum","Class","Order","Family","Genus","Species"), colnames(tax_tab))
    tax_df <- as.data.frame(tax_tab[names(cl_vec), tax_cols, drop = FALSE])
    cluster_map <- data.frame(
      ASV = names(cl_vec),
      ClusterID = paste0("Cluster_", cl_vec),
      tax_df,
      stringsAsFactors = FALSE
    )

    write.csv(cluster_map,
              sprintf("%s/cluster_map_similarity_%.3f_%s.csv",
                      dir_base, global_similarity_cutoff, date_tag),
              row.names = FALSE)

    cat(sprintf("[Go_intoASV] Similarity mode: cutoff=%.3f | cluster_map with taxonomy saved\n",
                global_similarity_cutoff),
        file = log_file, append = TRUE)

    # In similarity mode, each global sequence cluster becomes the analysis unit.
    tax_labels <- cluster_map$ClusterID
    names(tax_labels) <- cluster_map$ASV
    taxa_targets <- unique(cluster_map$ClusterID)

  } else {
    # ---------- regular taxonomy-based mode ----------
    if (is.null(level) || length(level) != 1L || !nzchar(trimws(level))) {
      stop("'level' must be provided in taxonomy mode when 'global_similarity_cutoff' is NULL.")
    }
    if (!level %in% colnames(tax_tab)) {
      stop(sprintf("'level' (%s) is not a column in tax_table(psIN).", level))
    }
    tax_raw <- tax_tab[, level, drop = TRUE]
    tax_raw[is.na(tax_raw)] <- "Unclassified"
    tax_labels <- trimws(tax_raw)
    names(tax_labels) <- rownames(tax_tab)
    taxa_pool <- unique(tax_labels)
    taxa_targets <- if (identical(target,"all")) taxa_pool else trimws(target)
    message("Target taxa: ", paste(taxa_targets, collapse = ", "))
  }

  # ---------- utils ----------
  trim_alignment <- function(aln, k) {
    if (k <= 0) return(aln)
    w <- Biostrings::width(aln)[1]
    if (k >= w / 2) return(aln)
    rng <- IRanges::IRanges(start = 1 + k, end = w - k)
    Biostrings::DNAStringSet(
      Biostrings::subseq(aln, start = IRanges::start(rng), end = IRanges::end(rng))
    )
  }
  entropy_weights <- function(p) {
    v <- p * (1 - p); if (sum(v) == 0) return(rep(0, length(p))); v / sum(v)
  }
  safe_tapply <- function(x, i, j, fun) {
    res <- try(tapply(x, list(i, j), fun), silent = TRUE)
    if (inherits(res, "try-error") || is.null(dimnames(res))) {
      res <- matrix(NA, nrow = length(unique(i)), ncol = length(unique(j)),
                    dimnames = list(unique(i), unique(j)))
    }
    res
  }

  # ---------- main loop ----------
  results <- mclapply(taxa_targets, function(target_taxon){
    # 1) 해당 cluster / taxon 에 속하는 ASV index
    idx <- which(tax_labels == target_taxon)
    if (length(idx) < min_asv) {
      cat(sprintf("Skip %s: < min_asv\n", target_taxon),
          file = log_file, append = TRUE)
      return(NULL)
    }

    # 2) abundance subset (ASV x sample → sample x ASV)
    sub_abund <- otu_tab[idx, , drop = FALSE]
    sub_abund <- t(sub_abund)

    # 3) per-sample abundance filter (min_abund)
    keep_samples <- rowSums(sub_abund) >= min_abund
    sub_abund <- sub_abund[keep_samples, , drop = FALSE]
    if (nrow(sub_abund) == 0) {
      cat(sprintf("Skip %s: total abundance < min_abund\n", target_taxon),
          file = log_file, append = TRUE)
      return(NULL)
    }

    # 4) sequences subset
    seqs_sub <- seqs_all[colnames(sub_abund)]
    seqs_sub <- seqs_sub[!is.na(seqs_sub)]
    if (any(!grepl("^[ACGTN]+$", seqs_sub))) {
      cat(sprintf("Skip %s: invalid sequences\n", target_taxon),
          file = log_file, append = TRUE)
      return(NULL)
    }
    seqs <- DNAStringSet(seqs_sub)

    # ----- intra-group clustering -----
    # Similarity mode uses the global cluster map only; no subclustering here.
    # Taxonomy mode keeps the largest near-identical subcluster within each taxon.
    if (is.null(global_similarity_cutoff) && !is.null(taxonomy_cluster_cutoff)) {
      dm <- DECIPHER::DistanceMatrix(seqs, includeTerminalGaps = FALSE)
      cl_vec <- cluster_from_dm(dm, taxonomy_cluster_cutoff)
      keep_cluster <- names(which.max(table(cl_vec)))
      keep_ids <- names(cl_vec[cl_vec == keep_cluster])
      seqs <- seqs[names(seqs) %in% keep_ids]
      sub_abund <- sub_abund[, colnames(sub_abund) %in% keep_ids, drop = FALSE]

      if (length(seqs) < min_asv) {
        cat(sprintf("Skip %s: < min_asv after clustering\n", target_taxon),
            file = log_file, append = TRUE)
        return(NULL)
      }
    } else if (!is.null(global_similarity_cutoff)) {
      cat(sprintf("Similarity mode active — skipping intra-group clustering for %s\n",
                  target_taxon),
          file = log_file, append = TRUE)
    }

    # ----- alignment -----
    aln <- tryCatch({
      if (aligner == "MAFFT" &&
          system("mafft --version", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
        tmp <- tempfile(fileext = ".fasta"); out <- tempfile(fileext = ".fasta")
        writeXStringSet(seqs, filepath = tmp)
        system(sprintf("mafft --auto --quiet %s > %s", tmp, out))
        readDNAStringSet(out)
      } else {
        AlignSeqs(seqs, iterations = 2, refinements = 2, verbose = FALSE)
      }
    }, error = function(e) {
      cat(sprintf("Align failed %s: %s\n", target_taxon, e$message),
          file = log_file, append = TRUE)
      NULL
    })
    if (is.null(aln)) return(NULL)

    # ----- trimming & gap QC -----
    if (trim_nt > 0) aln <- trim_alignment(aln, trim_nt)
    gap_prop <- tryCatch({
      sum(letterFrequency(aln, "-")) / (Biostrings::width(aln)[1] * length(aln))
    }, error = function(e) NA_real_)
    cat(sprintf("%s: Gap proportion %.2f%%\n",
                target_taxon, 100*gap_prop), file = log_file, append = TRUE)

    # ----- distance matrix -----
    if (method == "simple") {
      A <- as.matrix(aln)
      n <- nrow(A)
      dist_mat <- matrix(0, n, n, dimnames = list(rownames(A), rownames(A)))
      for (i in seq_len(n - 1)) {
        Ai <- A[i, ]
        for (j in (i + 1):n) {
          Aj <- A[j, ]
          keep <- if (distance_gap == "exclude") (Ai != '-') & (Aj != '-') else rep(TRUE, length(Ai))
          Lij <- sum(keep)
          if (Lij > 0) {
            dist_mat[i, j] <- sum(Ai[keep] != Aj[keep]) / Lij
          } else {
            dist_mat[i, j] <- NA_real_
          }
          dist_mat[j, i] <- dist_mat[i, j]
        }
      }
      diag(dist_mat) <- 0
    } else {
      aln_chr <- as.matrix(aln)
      aln_dna <- as.DNAbin(aln_chr)
      dist_mat <- ape::dist.dna(
        aln_dna,
        model = "TN93",              # nucdiv: TN93 + gamma
        gamma = TRUE,
        pairwise.deletion = (distance_gap == "exclude"),
        as.matrix = TRUE
      )
    }

    # ----- sync order -----
    ids <- intersect(colnames(sub_abund), colnames(dist_mat))
    if (length(ids) < 2) {
      cat(sprintf("Skip %s: <2 ids after sync\n", target_taxon),
          file = log_file, append = TRUE)
      return(NULL)
    }
    dist_mat  <- dist_mat[ids, ids, drop = FALSE]
    sub_abund <- sub_abund[, ids, drop = FALSE]
    L <- lower.tri(dist_mat, diag = FALSE)

    # ----- π calculator -----
    calc_pi <- function(a){
      tot <- sum(a)
      if (tot <= 0) return(list(pi = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
      p <- as.numeric(a) / tot
      w <- if (identical(weighting,"entropy")) entropy_weights(p) else p

      pi_val <- 2 * sum(outer(w, w)[L] * dist_mat[L], na.rm = TRUE)
      if (!is.finite(pi_val)) pi_val <- 0

      ci_low <- NA_real_; ci_high <- NA_real_
      if (compute_ci && sum(a > 0) >= 3) {
        boot_vals <- replicate(n_boot, {
          idx <- sample(seq_along(p), replace = TRUE)
          pb <- p[idx]; db <- dist_mat[idx, idx, drop = FALSE]
          wb <- if (identical(weighting,"entropy")) entropy_weights(pb) else pb
          2 * sum(outer(wb, wb)[lower.tri(db)] * db[lower.tri(db)], na.rm = TRUE)
        })
        if (any(is.finite(boot_vals))) {
          qs <- quantile(boot_vals, c(0.025, 0.975), na.rm = TRUE)
          ci_low <- as.numeric(qs[1]); ci_high <- as.numeric(qs[2])
        }
      }
      list(pi = as.numeric(pi_val), ci_low = ci_low, ci_high = ci_high)
    }

    # per-sample π 계산
    pi_list   <- lapply(rownames(sub_abund), function(s) calc_pi(sub_abund[s, ]))
    pi_values <- sapply(pi_list, function(x) x$pi)
    ci_low    <- sapply(pi_list, function(x) x$ci_low)
    ci_high   <- sapply(pi_list, function(x) x$ci_high)
    asv_counts <- rowSums(sub_abund > 0)

    df_out <- data.frame(
      Sample    = rownames(sub_abund),
      Taxon     = target_taxon,
      ASV_count = asv_counts,
      Pi        = as.numeric(pi_values),
      stringsAsFactors = FALSE
    )
    if (compute_ci) {
      df_out$Pi_CI_low  <- as.numeric(ci_low)
      df_out$Pi_CI_high <- as.numeric(ci_high)
    }
    df_out
  }, mc.cores = n_cores)

  # ---------- 결과 합치기 ----------
  results <- do.call(rbind, results)
  if (is.null(results) || nrow(results) == 0) {
    message("No valid taxa passed thresholds.")
    cat("No valid taxa passed thresholds\n", file = log_file, append = TRUE)
    return(psIN)
  }

  # Apply target filtering only in taxonomy mode.
  if (is.null(global_similarity_cutoff) && !identical(target,"all")) {
    results <- results[trimws(results$Taxon) == trimws(target), , drop = FALSE]
    if (nrow(results) == 0) {
      message("No rows for requested target.")
      return(psIN)
    }
  }

  uniq_idx <- !duplicated(results[, c("Sample","Taxon"), drop = FALSE])
  results  <- results[uniq_idx, , drop = FALSE]

  pi_tab   <- as.data.frame(safe_tapply(results$Pi,        results$Sample, results$Taxon, identity))
  asv_tab  <- as.data.frame(safe_tapply(results$ASV_count, results$Sample, results$Taxon, identity))

  if (is.null(global_similarity_cutoff) && !identical(target,"all")) {
    keep <- colnames(pi_tab) %in% trimws(target)
    pi_tab  <- pi_tab[,  keep, drop = FALSE]
    asv_tab <- asv_tab[, keep, drop = FALSE]
  }

  pi_tab  <- pi_tab[order(rownames(pi_tab)), , drop = FALSE]
  asv_tab <- asv_tab[order(rownames(asv_tab)), , drop = FALSE]
  colnames(pi_tab) <- paste0("pi_", colnames(pi_tab))

  write.csv(pi_tab,  pi_file,  row.names = TRUE, na = "")
  write.csv(asv_tab, asv_file, row.names = TRUE, na = "")

  # Optionally save CI matrices — now taken from results DF (robust)
  if (compute_ci && all(c("Pi_CI_low","Pi_CI_high") %in% names(results))) {
    ci_low_tab  <- as.data.frame(safe_tapply(results$Pi_CI_low,  results$Sample, results$Taxon, identity))
    ci_high_tab <- as.data.frame(safe_tapply(results$Pi_CI_high, results$Sample, results$Taxon, identity))
    if (is.null(global_similarity_cutoff) && !identical(target,"all")) {
      keep <- colnames(ci_low_tab) %in% trimws(target)
      ci_low_tab  <- ci_low_tab[,  keep, drop = FALSE]
      ci_high_tab <- ci_high_tab[, keep, drop = FALSE]
    }
    colnames(ci_low_tab)  <- paste0("pi_CI_low_",  colnames(ci_low_tab))
    colnames(ci_high_tab) <- paste0("pi_CI_high_", colnames(ci_high_tab))
    write.csv(ci_low_tab,  sprintf("%s/pi_ci_low_matrix_%s_%s_%s_%s.csv",  dir_base, method, level_tag, target_tag, date_tag),
              row.names = TRUE, na = "")
    write.csv(ci_high_tab, sprintf("%s/pi_ci_high_matrix_%s_%s_%s_%s.csv", dir_base, method, level_tag, target_tag, date_tag),
              row.names = TRUE, na = "")
  }

  end_time <- Sys.time()
  summary_text <- sprintf("
--------------------------------------------------
π summary (%s/%s): mean=%.5f | sd=%.5f | range=%.5f–%.5f
Seed=%d | Runtime=%.2f min
π matrix saved: %s
ASV count matrix saved: %s
--------------------------------------------------
", method, weighting,
                          mean(results$Pi, na.rm = TRUE),
                          sd(results$Pi,  na.rm = TRUE),
                          min(results$Pi, na.rm = TRUE),
                          max(results$Pi, na.rm = TRUE),
                          seed,
                          round(as.numeric(difftime(end_time, start_time, units='mins')), 2),
                          pi_file, asv_file)
  cat(summary_text)
  cat(summary_text, file = log_file, append = TRUE)

  # ---------- merge π & ASV_count into sample_data ----------
  if (!is.null(sample_data(psIN, errorIfNULL = FALSE))) {
    sd <- as.data.frame(sample_data(psIN))

    pi_cols  <- pi_tab[rownames(sd), , drop = FALSE]
    asv_cols <- asv_tab[rownames(sd), , drop = FALSE]
    colnames(asv_cols) <- paste0("asvN_", colnames(asv_cols))

    if (ncol(pi_cols) > 0)  pi_cols[]  <- lapply(pi_cols,  function(x) suppressWarnings(as.numeric(x)))
    if (ncol(asv_cols) > 0) asv_cols[] <- lapply(asv_cols, function(x) suppressWarnings(as.numeric(x)))

    make_uniq <- function(df) { colnames(df) <- make.unique(colnames(df)); df }
    pi_cols  <- make_uniq(pi_cols)
    asv_cols <- make_uniq(asv_cols)

    if (ncol(pi_cols) > 0) {
      sd$pi_global <- rowMeans(pi_cols, na.rm = TRUE)
      sd$pi_global[!is.finite(sd$pi_global)] <- NA_real_
    }

    sd_merged <- cbind(sd, pi_cols, asv_cols)
    sample_data(psIN) <- sd_merged

    cat(sprintf("\n[merge] Added to sample_data: %d pi cols, %d asv-count cols\n",
                ncol(pi_cols), ncol(asv_cols)))
    if (ncol(pi_cols) > 0) {
      cat("[Go_intoASV] pi_global added to sample_data (mean across pi_* columns)\n")
    }
  } else {
    warning("sample_data(psIN) is NULL — cannot merge π results into sample metadata.")
  }

  return(psIN)
}
