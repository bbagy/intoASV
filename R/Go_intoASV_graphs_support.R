#' intoASV sequence-cloud pi support analysis
#'
#' Runs present-only pi analysis, sample-level sequence evidence, and
#' integrated summary plots from a \code{Go_intoASV()} result. Companion to
#' \code{Go_intoASV_graphs()}.
#'
#' @param psIN A \code{phyloseq} object returned by \code{Go_intoASV()}.
#' @param cluster_map Cluster map data frame or CSV path (must contain \code{ASV} and \code{ClusterID}).
#' @param project Character; project prefix (should embed similarity cutoff and weighting).
#' @param global_similarity_cutoff Numeric; similarity cutoff (e.g., 0.970). If \code{NULL},
#'   extracted from the project name.
#' @param mainGroup Character; group variable in \code{sample_data} for comparison.
#' @param order Character vector; group levels in comparison order (first = reference).
#' @param covariate Character; additional covariate for linear model adjustment.
#' @param inset_layout Optional list of cluster ID vectors (same as in \code{Go_intoASV_graphs()});
#'   used as \code{target_clusters} when \code{target_clusters} is \code{NULL}.
#' @param target_clusters Optional character vector of cluster IDs for sample-level evidence.
#' @param min_present_per_group Minimum samples with pi > 0 per group for a cluster to be tested.
#' @param min_present_total Minimum total present samples for a cluster to be tested.
#' @param n_perm Number of permutations for pi delta test.
#' @param n_boot Number of bootstrap replicates for pi delta.
#' @param graph_n_perm Number of permutations for graph metric test (0 = skip).
#' @param graph_n_boot Number of bootstrap replicates for graph metrics.
#' @param forest_width,support_width,evidence_width PDF widths in inches.
#' @param support_top_clusters Fallback: top N clusters by minimum permutation p for support plot.
#' @param seed Random seed.
#'
#' @return Invisible list with statistical results, evidence data, plot objects, and output paths.
#'
#' @examples
#' \dontrun{
#' support_res <- Go_intoASV_graphs_support(
#'   psIN        = ps_intoASV,
#'   cluster_map = cluster_map,
#'   project     = project_i,
#'   mainGroup   = "case_control",
#'   order       = c("control", "case"),
#'   covariate   = "visit_y",
#'   inset_layout = sequence_cloud_inset_layout
#' )
#' }
#'
#' @export
Go_intoASV_graphs_support <- function(
    psIN,
    cluster_map,
    project               = "intoASV_graphs_support",
    global_similarity_cutoff = NULL,
    mainGroup             = NULL,
    order                 = NULL,
    covariate             = NULL,
    inset_layout          = NULL,
    target_clusters       = NULL,
    min_present_per_group = 3,
    min_present_total     = 6,
    n_perm                = 499,
    n_boot                = 500,
    graph_n_perm          = 0,
    graph_n_boot          = 200,
    forest_width          = 8,
    support_width         = 11,
    evidence_width        = 12,
    support_top_clusters  = 12,
    seed                  = 123
) {
  if (!inherits(psIN, "phyloseq")) {
    stop("[Go_intoASV_graphs_support] psIN must be a phyloseq object.")
  }

  ###############################################
  # 0. Setup
  ###############################################
  date_tag <- format(Sys.Date(), "%y%m%d")

  .extract_weighting <- function(x) {
    x <- tolower(paste(x, collapse = " "))
    hit <- regmatches(x, regexpr("(abundance|entropy)", x, perl = TRUE))
    if (length(hit) == 0 || hit == "") "weightingNA" else hit
  }

  .extract_cutoff <- function(x) {
    x <- paste(x, collapse = " ")
    hit <- regmatches(x, regexpr("similarity[_-]?([0-9]+\\.?[0-9]*)", x, perl = TRUE))
    if (length(hit) == 0 || hit == "") return(NULL)
    val_str <- sub(".*similarity[_-]?", "", hit)
    val_num <- suppressWarnings(as.numeric(val_str))
    if (!is.finite(val_num)) return(NULL)
    if (val_num > 1) val_num <- val_num / 1000
    val_num
  }

  weighting_tag <- .extract_weighting(project)

  if (is.null(global_similarity_cutoff)) {
    global_similarity_cutoff <- .extract_cutoff(project)
    if (is.null(global_similarity_cutoff)) {
      stop("[Go_intoASV_graphs_support] Could not extract similarity cutoff from project name. ",
           "Provide global_similarity_cutoff explicitly.")
    }
  }
  cutoff_i <- global_similarity_cutoff

  group_order <- if (!is.null(order)) as.character(order) else character(0)
  if (!is.null(mainGroup) && length(group_order) == 0) {
    sd_df <- data.frame(phyloseq::sample_data(psIN), check.names = FALSE)
    group_order <- sort(unique(as.character(sd_df[[mainGroup]])))
    warning("[Go_intoASV_graphs_support] 'order' not provided; using alphabetical group order: ",
            paste(group_order, collapse = ", "))
  }
  group_1 <- if (length(group_order) >= 1) group_order[1] else NA_character_
  group_2 <- if (length(group_order) >= 2) group_order[2] else NA_character_

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  dir_base <- sprintf("%s_%s/intoASV/sequence_cloud_pi", project, date_tag)
  dir.create(dir_base, recursive = TRUE, showWarnings = FALSE)

  ###############################################
  # 1. Resolve cluster_map
  ###############################################
  cluster_map_file <- NA_character_
  if (is.character(cluster_map) && length(cluster_map) == 1) {
    cluster_map_file <- cluster_map
    if (!file.exists(cluster_map_file)) {
      stop("[Go_intoASV_graphs_support] cluster_map file does not exist: ", cluster_map_file)
    }
    cluster_map <- utils::read.csv(cluster_map_file, check.names = FALSE)
  }
  if (!all(c("ASV", "ClusterID") %in% colnames(cluster_map))) {
    stop("[Go_intoASV_graphs_support] cluster_map must contain ASV and ClusterID columns.")
  }
  if (!"Genus"   %in% colnames(cluster_map)) cluster_map$Genus   <- NA_character_
  if (!"Species" %in% colnames(cluster_map)) cluster_map$Species <- NA_character_

  ###############################################
  # 2. Internal helpers
  ###############################################
  .mode_or_na <- function(x) {
    x <- x[!is.na(x) & nzchar(as.character(x))]
    if (length(x) == 0) return(NA_character_)
    names(sort(table(x), decreasing = TRUE))[1]
  }

  .find_pi_file <- function(file_prefix) {
    candidates <- Sys.glob(sprintf(
      "%s_*/intoASV/pi_tab/%s*similarity_%.3f_*.csv",
      project, file_prefix, cutoff_i
    ))
    if (length(candidates) == 0) {
      stop("[Go_intoASV_graphs_support] No file found for prefix '", file_prefix,
           "' with similarity ", sprintf("%.3f", cutoff_i))
    }
    candidates[which.max(file.info(candidates)$mtime)]
  }

  .find_cloud_stats <- function(group_var) {
    if (is.null(group_var)) return(NA_character_)
    candidates <- Sys.glob(sprintf(
      "%s_*/intoASV/sequence_cloud/sequence_cloud_stats_%s_*.csv",
      project, group_var
    ))
    if (length(candidates) == 0) return(NA_character_)
    candidates[which.max(file.info(candidates)$mtime)]
  }

  .get_otu_matrix <- function(ps_obj) {
    otu_tab <- as(phyloseq::otu_table(ps_obj), "matrix")
    if (!phyloseq::taxa_are_rows(ps_obj)) otu_tab <- t(otu_tab)
    otu_tab
  }

  .get_asv_sequences <- function(ps_obj) {
    seqs_tmp <- NULL
    ref <- phyloseq::refseq(ps_obj, errorIfNULL = FALSE)
    if (!is.null(ref)) {
      ref <- as.character(ref)
      if (all(grepl("^[ACGTN]+$", ref))) seqs_tmp <- ref
    }
    if (is.null(seqs_tmp)) {
      tx <- phyloseq::taxa_names(ps_obj)
      if (all(grepl("^[ACGTN]+$", tx))) seqs_tmp <- tx
    }
    if (is.null(seqs_tmp)) {
      rn <- rownames(as(phyloseq::otu_table(ps_obj), "matrix"))
      if (all(grepl("^[ACGTN]+$", rn))) seqs_tmp <- rn
    }
    if (is.null(seqs_tmp)) stop("[Go_intoASV_graphs_support] No valid DNA sequences found.")
    names(seqs_tmp) <- phyloseq::taxa_names(ps_obj)
    seqs_tmp
  }

  ###############################################
  # 3. Load data
  ###############################################
  pi_file          <- .find_pi_file("pi_matrix_nucdiv")
  count_file       <- .find_pi_file("asv_count_matrix")
  cloud_stats_file <- .find_cloud_stats(mainGroup)

  pi_mat    <- utils::read.csv(pi_file,    row.names = 1, check.names = FALSE)
  count_mat <- utils::read.csv(count_file, row.names = 1, check.names = FALSE)
  sample_df <- data.frame(phyloseq::sample_data(psIN), check.names = FALSE)
  sample_df$SampleID <- rownames(sample_df)

  if (!is.null(mainGroup) && !mainGroup %in% colnames(sample_df)) {
    stop("[Go_intoASV_graphs_support] mainGroup '", mainGroup, "' not found in sample_data.")
  }

  common_samples <- Reduce(intersect, list(rownames(pi_mat), rownames(count_mat), rownames(sample_df)))
  pi_mat    <- pi_mat[common_samples, , drop = FALSE]
  count_mat <- count_mat[common_samples, , drop = FALSE]
  sample_df <- sample_df[common_samples, , drop = FALSE]

  cluster_annot <- cluster_map %>%
    dplyr::group_by(ClusterID) %>%
    dplyr::summarise(
      cluster_size = dplyr::n(),
      Genus        = .mode_or_na(Genus),
      Species      = .mode_or_na(Species),
      .groups      = "drop"
    )

  cluster_ids <- intersect(
    sub("^pi_", "", grep("^pi_Cluster_", colnames(pi_mat), value = TRUE)),
    colnames(count_mat)
  )

  pi_long <- pi_mat[, paste0("pi_", cluster_ids), drop = FALSE] %>%
    as.data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column("SampleID") %>%
    tidyr::pivot_longer(-SampleID, names_to = "ClusterID", values_to = "pi") %>%
    dplyr::mutate(ClusterID = sub("^pi_", "", ClusterID))

  count_long <- count_mat[, cluster_ids, drop = FALSE] %>%
    as.data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column("SampleID") %>%
    tidyr::pivot_longer(-SampleID, names_to = "ClusterID", values_to = "asv_count")

  df_pi <- pi_long %>%
    dplyr::left_join(count_long, by = c("SampleID", "ClusterID")) %>%
    dplyr::left_join(sample_df, by = "SampleID") %>%
    dplyr::mutate(present = pi > 0 & asv_count > 1)

  if (!is.null(mainGroup)) {
    df_pi <- df_pi %>%
      dplyr::mutate(!!mainGroup := factor(as.character(.data[[mainGroup]]), levels = group_order)) %>%
      dplyr::filter(!is.na(.data[[mainGroup]]))
  }
  if (!is.null(covariate) && covariate %in% colnames(df_pi)) {
    df_pi[[covariate]] <- factor(df_pi[[covariate]])
  }

  ###############################################
  # 4. Pi present summary & wide comparison
  ###############################################
  pi_present_summary <- df_pi %>%
    dplyr::group_by(
      ClusterID,
      group = if (!is.null(mainGroup)) .data[[mainGroup]] else factor("all")
    ) %>%
    dplyr::summarise(
      n_samples              = dplyr::n(),
      n_present              = sum(present, na.rm = TRUE),
      mean_pi_all            = mean(pi, na.rm = TRUE),
      mean_pi_present        = mean(pi[present], na.rm = TRUE),
      median_pi_present      = stats::median(pi[present], na.rm = TRUE),
      mean_asv_count_present = mean(asv_count[present], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(cluster_annot, by = "ClusterID")

  if (!is.na(cloud_stats_file) && file.exists(cloud_stats_file)) {
    cloud_stats <- utils::read.csv(cloud_stats_file, check.names = FALSE)
    if (!is.null(mainGroup) && mainGroup %in% colnames(cloud_stats)) {
      cloud_stats <- cloud_stats %>%
        dplyr::rename(group = dplyr::all_of(mainGroup)) %>%
        dplyr::select(ClusterID, group,
          graph_n_nodes                = n_nodes,
          graph_modularity             = modularity,
          graph_dispersion             = cloud_dispersion,
          graph_mean_pairwise_distance = mean_pairwise_distance
        )
      pi_present_summary <- pi_present_summary %>%
        dplyr::left_join(cloud_stats, by = c("ClusterID", "group"))
    }
  }

  utils::write.csv(pi_present_summary,
    file.path(dir_base, sprintf("pi_present_summary_%s_%s.csv", weighting_tag, date_tag)),
    row.names = FALSE)

  compare_value_cols <- intersect(
    c("n_present", "mean_pi_all", "mean_pi_present", "median_pi_present",
      "mean_asv_count_present", "graph_n_nodes", "graph_modularity",
      "graph_dispersion", "graph_mean_pairwise_distance"),
    colnames(pi_present_summary)
  )

  pi_present_compare <- pi_present_summary %>%
    dplyr::filter(group %in% group_order[seq_len(min(2, length(group_order)))]) %>%
    dplyr::select(ClusterID, group, Species, Genus, cluster_size,
                  dplyr::all_of(compare_value_cols)) %>%
    tidyr::pivot_wider(names_from = group, values_from = dplyr::all_of(compare_value_cols))

  if (length(group_order) >= 2) {
    for (metric in c("mean_pi_present", "graph_dispersion", "graph_modularity")) {
      c1 <- paste0(metric, "_", group_order[1])
      c2 <- paste0(metric, "_", group_order[2])
      if (all(c(c1, c2) %in% colnames(pi_present_compare))) {
        pi_present_compare[[paste0("delta_", metric)]] <- pi_present_compare[[c2]] - pi_present_compare[[c1]]
      }
    }
  }

  utils::write.csv(pi_present_compare,
    file.path(dir_base, sprintf("pi_present_compare_%s_%s.csv", weighting_tag, date_tag)),
    row.names = FALSE)

  ###############################################
  # 5. Pi present lm
  ###############################################
  cluster_keep_present <- df_pi %>%
    dplyr::group_by(ClusterID) %>%
    dplyr::summarise(
      nonzero_n      = sum(present, na.rm = TRUE),
      nonzero_group1 = if (!is.null(mainGroup)) sum(present & .data[[mainGroup]] == group_1, na.rm = TRUE) else nonzero_n,
      nonzero_group2 = if (!is.null(mainGroup)) sum(present & .data[[mainGroup]] == group_2, na.rm = TRUE) else 0L,
      .groups = "drop"
    ) %>%
    dplyr::filter(
      nonzero_n      >= min_present_total,
      nonzero_group1 >= min_present_per_group,
      nonzero_group2 >= min_present_per_group
    )

  .run_lm_present <- function(d) {
    d2 <- d %>%
      dplyr::filter(present, pi > 0) %>%
      dplyr::mutate(pi_log10 = log10(pi))
    if (nrow(d2) < min_present_total) return(data.frame())
    if (!is.null(mainGroup) && dplyr::n_distinct(d2[[mainGroup]]) < 2) return(data.frame())
    terms_i <- if (!is.null(mainGroup)) mainGroup else character(0)
    if (!is.null(covariate) && covariate %in% colnames(d2) && dplyr::n_distinct(d2[[covariate]]) >= 2) {
      terms_i <- c(terms_i, covariate)
    }
    if (length(terms_i) == 0) return(data.frame())
    fit <- stats::lm(stats::reformulate(terms_i, response = "pi_log10"), data = d2)
    broom::tidy(fit, conf.int = TRUE)
  }

  res_lm_present <- df_pi %>%
    dplyr::filter(ClusterID %in% cluster_keep_present$ClusterID) %>%
    dplyr::group_by(ClusterID) %>%
    dplyr::group_modify(~ .run_lm_present(.x)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(p.adj = stats::p.adjust(p.value, method = "fdr")) %>%
    dplyr::left_join(cluster_keep_present, by = "ClusterID") %>%
    dplyr::left_join(cluster_annot, by = "ClusterID")

  utils::write.csv(res_lm_present,
    file.path(dir_base, sprintf("pi_present_lm_all_terms_%s_%s.csv", weighting_tag, date_tag)),
    row.names = FALSE)

  res_case_present <- data.frame()
  if (!is.null(mainGroup) && length(group_order) >= 2 && nrow(res_lm_present) > 0) {
    target_term <- paste0(mainGroup, group_order[2])
    res_case_present <- res_lm_present %>%
      dplyr::filter(term == target_term) %>%
      dplyr::mutate(
        sig_cat = dplyr::case_when(
          p.adj < 0.1 & estimate > 0 ~ "FDR+, Positive",
          p.adj < 0.1 & estimate < 0 ~ "FDR+, Negative",
          p.adj >= 0.1 & p.value < 0.05 ~ "Nominal p<0.05",
          TRUE ~ "NS"
        ),
        sig_cat = factor(sig_cat, levels = c("FDR+, Positive", "FDR+, Negative", "Nominal p<0.05", "NS")),
        Cluster_label = ifelse(!is.na(Species) & nzchar(Species),
                               paste0(ClusterID, " (", Species, ")"), ClusterID)
      ) %>%
      dplyr::arrange(estimate) %>%
      dplyr::mutate(Cluster_label = forcats::fct_inorder(Cluster_label))

    utils::write.csv(res_case_present,
      file.path(dir_base, sprintf("pi_present_lm_%s_%s_%s.csv", mainGroup, weighting_tag, date_tag)),
      row.names = FALSE)
  }

  ###############################################
  # 6. Pi permutation / bootstrap
  ###############################################
  .pi_delta <- function(d, labels = NULL) {
    if (!is.null(labels)) d[[mainGroup]] <- labels
    if (is.null(mainGroup)) return(NA_real_)
    if (sum(d[[mainGroup]] == group_1, na.rm = TRUE) < min_present_per_group) return(NA_real_)
    if (sum(d[[mainGroup]] == group_2, na.rm = TRUE) < min_present_per_group) return(NA_real_)
    pi_log10 <- log10(d$pi)
    mean(pi_log10[d[[mainGroup]] == group_2], na.rm = TRUE) -
      mean(pi_log10[d[[mainGroup]] == group_1], na.rm = TRUE)
  }

  pi_perm_boot <- data.frame()
  if (!is.null(mainGroup) && nrow(cluster_keep_present) > 0) {
    pi_perm_boot <- dplyr::bind_rows(lapply(cluster_keep_present$ClusterID, function(cl) {
      d_present <- df_pi %>%
        dplyr::filter(ClusterID == cl, present, pi > 0,
                      .data[[mainGroup]] %in% c(group_1, group_2))
      obs <- .pi_delta(d_present)
      if (!is.finite(obs)) return(NULL)
      perm_d <- replicate(n_perm, .pi_delta(d_present, labels = sample(d_present[[mainGroup]])))
      idx_1 <- which(d_present[[mainGroup]] == group_1)
      idx_2 <- which(d_present[[mainGroup]] == group_2)
      boot_d <- replicate(n_boot, {
        .pi_delta(d_present[c(sample(idx_1, length(idx_1), replace = TRUE),
                              sample(idx_2, length(idx_2), replace = TRUE)), , drop = FALSE])
      })
      perm_d <- perm_d[is.finite(perm_d)]
      boot_d <- boot_d[is.finite(boot_d)]
      data.frame(
        ClusterID = cl, metric = "pi_log10_present", delta_observed = obs,
        permutation_n = length(perm_d),
        permutation_p_two_sided = ifelse(length(perm_d) > 0,
          (1 + sum(abs(perm_d) >= abs(obs), na.rm = TRUE)) / (length(perm_d) + 1), NA_real_),
        bootstrap_n = length(boot_d),
        bootstrap_ci_low  = ifelse(length(boot_d) > 0, stats::quantile(boot_d, 0.025, na.rm = TRUE), NA_real_),
        bootstrap_ci_high = ifelse(length(boot_d) > 0, stats::quantile(boot_d, 0.975, na.rm = TRUE), NA_real_),
        bootstrap_same_direction = ifelse(length(boot_d) > 0,
          mean(sign(boot_d) == sign(obs), na.rm = TRUE), NA_real_),
        stringsAsFactors = FALSE
      )
    })) %>%
      dplyr::left_join(cluster_annot, by = "ClusterID")

    utils::write.csv(pi_perm_boot,
      file.path(dir_base, sprintf("pi_present_permutation_bootstrap_%s_%s.csv", weighting_tag, date_tag)),
      row.names = FALSE)
  }

  ###############################################
  # 7. Graph metric permutation / bootstrap
  ###############################################
  graph_perm_boot <- data.frame()
  if (!is.null(mainGroup) && (graph_n_perm > 0 || graph_n_boot > 0)) {
    otu_tab  <- .get_otu_matrix(psIN)
    otu_tab  <- otu_tab[, common_samples, drop = FALSE]
    seqs_all <- .get_asv_sequences(psIN)

    graph_clusters <- if (!is.na(cloud_stats_file) && file.exists(cloud_stats_file)) {
      unique(utils::read.csv(cloud_stats_file, check.names = FALSE)$ClusterID)
    } else {
      unique(cluster_keep_present$ClusterID)
    }
    graph_clusters <- intersect(graph_clusters, unique(cluster_map$ClusterID))

    .build_cache <- function() {
      caches <- lapply(graph_clusters, function(cl) {
        asv_ids <- intersect(cluster_map$ASV[cluster_map$ClusterID == cl], rownames(otu_tab))
        asv_ids <- intersect(asv_ids, names(seqs_all))
        if (length(asv_ids) < 2) return(NULL)
        seqs <- Biostrings::DNAStringSet(seqs_all[asv_ids])
        dm <- as.matrix(DECIPHER::DistanceMatrix(seqs, includeTerminalGaps = FALSE))
        list(ClusterID = cl, dm = dm[asv_ids, asv_ids, drop = FALSE],
             otu = otu_tab[asv_ids, , drop = FALSE])
      })
      names(caches) <- graph_clusters
      caches[!vapply(caches, is.null, logical(1))]
    }

    .calc_graph_metrics <- function(cache_i, sample_ids) {
      sample_ids <- sample_ids[sample_ids %in% colnames(cache_i$otu)]
      if (length(sample_ids) == 0) return(NULL)
      abund <- rowSums(cache_i$otu[, sample_ids, drop = FALSE])
      ids <- intersect(names(abund)[abund > 0], rownames(cache_i$dm))
      if (length(ids) < 2) return(NULL)
      dm <- cache_i$dm[ids, ids, drop = FALSE]
      full_g <- igraph::graph_from_adjacency_matrix(dm, mode = "undirected", weighted = TRUE, diag = FALSE)
      mst_edges <- as.data.frame(igraph::as_edgelist(igraph::mst(full_g, weights = igraph::E(full_g)$weight)),
                                 stringsAsFactors = FALSE)
      if (ncol(mst_edges) == 2) {
        colnames(mst_edges) <- c("from", "to"); mst_edges$edge_type <- "mst"
      } else {
        mst_edges <- data.frame(from = character(), to = character(), edge_type = character())
      }
      knn_edges <- do.call(rbind, lapply(rownames(dm), function(i) {
        x <- dm[i, ]; x[i] <- Inf
        j <- names(sort(x))[seq_len(min(3, length(x) - 1))]
        data.frame(from = i, to = j, edge_type = "knn", stringsAsFactors = FALSE)
      }))
      edge_df <- rbind(mst_edges, knn_edges)
      if (nrow(edge_df) == 0) return(NULL)
      pair <- ifelse(edge_df$from < edge_df$to,
                     paste(edge_df$from, edge_df$to, sep = "__"),
                     paste(edge_df$to, edge_df$from, sep = "__"))
      edge_df <- edge_df[!duplicated(pair), , drop = FALSE]
      edge_df$inv_distance <- 1 / (mapply(function(a, b) dm[a, b], edge_df$from, edge_df$to) + 1e-6)
      nodes <- data.frame(name = ids, total_abundance = as.numeric(abund[ids]), stringsAsFactors = FALSE)
      g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = nodes)
      igraph::E(g)$inv_distance <- edge_df$inv_distance
      set.seed(seed)
      comm <- tryCatch(igraph::cluster_louvain(g, weights = igraph::E(g)$inv_distance), error = function(e) NULL)
      set.seed(seed)
      layout <- igraph::layout_with_fr(g, weights = igraph::E(g)$inv_distance, niter = 1000)
      radial <- sqrt((layout[, 1] - mean(layout[, 1]))^2 + (layout[, 2] - mean(layout[, 2]))^2)
      data.frame(
        n_nodes    = igraph::vcount(g),
        modularity = if (!is.null(comm)) igraph::modularity(comm) else NA_real_,
        dispersion = mean(radial, na.rm = TRUE),
        mean_pairwise_distance = mean(dm[lower.tri(dm)], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }

    .graph_delta_table <- function(s1, s2, caches) {
      dplyr::bind_rows(lapply(names(caches), function(cl) {
        m1 <- .calc_graph_metrics(caches[[cl]], s1)
        m2 <- .calc_graph_metrics(caches[[cl]], s2)
        if (is.null(m1) || is.null(m2)) return(NULL)
        data.frame(ClusterID = cl,
          delta_graph_dispersion             = m2$dispersion - m1$dispersion,
          delta_graph_modularity             = m2$modularity - m1$modularity,
          delta_graph_mean_pairwise_distance = m2$mean_pairwise_distance - m1$mean_pairwise_distance,
          n_nodes_group1 = m1$n_nodes, n_nodes_group2 = m2$n_nodes,
          stringsAsFactors = FALSE)
      }))
    }

    group_samples <- rownames(sample_df)[sample_df[[mainGroup]] %in% c(group_1, group_2)]
    group_labels  <- as.character(sample_df[group_samples, mainGroup])
    s1 <- group_samples[group_labels == group_1]
    s2 <- group_samples[group_labels == group_2]

    graph_cache <- .build_cache()
    graph_obs   <- .graph_delta_table(s1, s2, graph_cache)

    if (nrow(graph_obs) > 0) {
      empty_graph_df <- data.frame(ClusterID = character(), delta_graph_dispersion = numeric(),
        delta_graph_modularity = numeric(), delta_graph_mean_pairwise_distance = numeric(),
        iteration = integer())

      graph_perm <- if (graph_n_perm > 0) {
        dplyr::bind_rows(lapply(seq_len(graph_n_perm), function(i) {
          pl <- sample(group_labels)
          out <- .graph_delta_table(group_samples[pl == group_1], group_samples[pl == group_2], graph_cache)
          if (nrow(out) == 0) return(NULL)
          out$iteration <- i; out
        }))
      } else { empty_graph_df }

      graph_boot <- if (graph_n_boot > 0) {
        dplyr::bind_rows(lapply(seq_len(graph_n_boot), function(i) {
          out <- .graph_delta_table(sample(s1, length(s1), replace = TRUE),
                                    sample(s2, length(s2), replace = TRUE), graph_cache)
          if (nrow(out) == 0) return(NULL)
          out$iteration <- i; out
        }))
      } else { empty_graph_df }

      graph_metrics <- c("graph_dispersion", "graph_modularity", "graph_mean_pairwise_distance")
      graph_perm_boot <- dplyr::bind_rows(lapply(seq_len(nrow(graph_obs)), function(i) {
        cl <- graph_obs$ClusterID[i]
        dplyr::bind_rows(lapply(graph_metrics, function(metric_i) {
          dc   <- paste0("delta_", metric_i)
          obs  <- graph_obs[[dc]][i]
          pd   <- graph_perm[[dc]][graph_perm$ClusterID == cl]; pd <- pd[is.finite(pd)]
          bd   <- graph_boot[[dc]][graph_boot$ClusterID == cl]; bd <- bd[is.finite(bd)]
          data.frame(ClusterID = cl, metric = metric_i, delta_observed = obs,
            permutation_n = length(pd),
            permutation_p_two_sided = ifelse(length(pd) > 0,
              (1 + sum(abs(pd) >= abs(obs), na.rm = TRUE)) / (length(pd) + 1), NA_real_),
            bootstrap_n = length(bd),
            bootstrap_ci_low  = ifelse(length(bd) > 0, stats::quantile(bd, 0.025, na.rm = TRUE), NA_real_),
            bootstrap_ci_high = ifelse(length(bd) > 0, stats::quantile(bd, 0.975, na.rm = TRUE), NA_real_),
            bootstrap_same_direction = ifelse(length(bd) > 0,
              mean(sign(bd) == sign(obs), na.rm = TRUE), NA_real_),
            stringsAsFactors = FALSE)
        }))
      })) %>%
        dplyr::left_join(cluster_annot, by = "ClusterID")

      utils::write.csv(graph_perm_boot,
        file.path(dir_base, sprintf("graph_permutation_bootstrap_%s_%s.csv", weighting_tag, date_tag)),
        row.names = FALSE)
    }
  }

  ###############################################
  # 8. Sample evidence analysis
  ###############################################
  target_cl <- if (!is.null(target_clusters)) {
    as.character(target_clusters)
  } else if (!is.null(inset_layout)) {
    unique(unlist(inset_layout, use.names = FALSE))
  } else {
    unique(cluster_map$ClusterID)
  }
  target_cl <- intersect(as.character(target_cl), unique(cluster_map$ClusterID))

  sample_metric_df       <- data.frame()
  prevalence_compare     <- data.frame()
  sample_metric_perm_boot <- data.frame()

  if (length(target_cl) > 0 && !is.null(mainGroup)) {
    otu_se   <- .get_otu_matrix(psIN)
    seqs_se  <- .get_asv_sequences(psIN)
    sids     <- intersect(colnames(otu_se), rownames(sample_df))
    otu_se   <- otu_se[, sids, drop = FALSE]
    sdf_se   <- sample_df[sids, , drop = FALSE]
    sdf_se[[mainGroup]] <- factor(as.character(sdf_se[[mainGroup]]), levels = group_order)
    sdf_se   <- sdf_se[!is.na(sdf_se[[mainGroup]]), , drop = FALSE]

    .calc_cluster_sample_metrics <- function(cluster_i) {
      asv_ids <- intersect(cluster_map$ASV[cluster_map$ClusterID == cluster_i], rownames(otu_se))
      asv_ids <- intersect(asv_ids, names(seqs_se))
      if (length(asv_ids) == 0) return(NULL)
      dm <- NULL
      if (length(asv_ids) >= 2) {
        seqs <- Biostrings::DNAStringSet(seqs_se[asv_ids])
        dm <- as.matrix(DECIPHER::DistanceMatrix(seqs, includeTerminalGaps = FALSE))
        dm <- dm[asv_ids, asv_ids, drop = FALSE]
      }
      do.call(rbind, lapply(rownames(sdf_se), function(sid) {
        abund <- otu_se[asv_ids, sid]; abund <- abund[abund > 0]
        detected        <- length(abund) > 0
        diverse_present <- length(abund) > 1
        mpd <- NA_real_; awpd <- NA_real_; mnnd <- NA_real_; wce <- NA_real_; wcev <- NA_real_
        if (detected) {
          rel <- as.numeric(abund) / sum(abund)
          wce  <- -sum(rel * log(rel), na.rm = TRUE)
          wcev <- ifelse(length(abund) > 1, wce / log(length(abund)), NA_real_)
        }
        if (diverse_present && !is.null(dm)) {
          ids   <- names(abund)
          dm_s  <- dm[ids, ids, drop = FALSE]
          ld    <- dm_s[lower.tri(dm_s)]
          mpd   <- mean(ld, na.rm = TRUE)
          w     <- outer(as.numeric(abund), as.numeric(abund), "*")
          lw    <- w[lower.tri(w)]
          awpd  <- if (sum(lw, na.rm = TRUE) > 0) stats::weighted.mean(ld, lw, na.rm = TRUE) else NA_real_
          diag(dm_s) <- Inf
          mnnd  <- mean(apply(dm_s, 1, min, na.rm = TRUE), na.rm = TRUE)
        }
        data.frame(SampleID = sid, ClusterID = cluster_i,
          detected = detected, diverse_present = diverse_present,
          total_abundance = sum(abund), asv_richness = length(abund),
          mean_pairwise_distance = mpd, abundance_weighted_pairwise_distance = awpd,
          mean_nearest_neighbor_distance = mnnd, within_cluster_entropy = wce,
          within_cluster_evenness = wcev, stringsAsFactors = FALSE)
      }))
    }

    sample_metric_df <- dplyr::bind_rows(lapply(target_cl, .calc_cluster_sample_metrics)) %>%
      dplyr::left_join(sdf_se, by = "SampleID") %>%
      dplyr::left_join(cluster_annot, by = "ClusterID")

    utils::write.csv(sample_metric_df,
      file.path(dir_base, sprintf("sample_metrics_%s_%s.csv", weighting_tag, date_tag)),
      row.names = FALSE)

    .fisher_flag <- function(d, flag_col) {
      gr   <- as.character(d[[mainGroup]])
      flag <- as.logical(d[[flag_col]])
      d2   <- data.frame(group = gr, flag = flag)[gr %in% group_order[seq_len(2)], ]
      if (nrow(d2) == 0) return(NA_real_)
      x1 <- sum(d2$flag[d2$group == group_order[1]], na.rm = TRUE)
      n1 <- sum(d2$group == group_order[1])
      x2 <- sum(d2$flag[d2$group == group_order[2]], na.rm = TRUE)
      n2 <- sum(d2$group == group_order[2])
      if (n1 == 0 || n2 == 0) return(NA_real_)
      tryCatch(
        stats::fisher.test(matrix(c(x1, n1 - x1, x2, n2 - x2), nrow = 2, byrow = TRUE))$p.value,
        error = function(e) NA_real_
      )
    }

    prevalence_by_group <- sample_metric_df %>%
      dplyr::group_by(ClusterID, group = .data[[mainGroup]]) %>%
      dplyr::summarise(
        n_samples = dplyr::n(), n_detected = sum(detected, na.rm = TRUE),
        n_diverse_present = sum(diverse_present, na.rm = TRUE),
        detected_prevalence = n_detected / n_samples,
        diverse_prevalence  = n_diverse_present / n_samples,
        .groups = "drop"
      ) %>%
      dplyr::left_join(cluster_annot, by = "ClusterID")

    prevalence_tests <- dplyr::bind_rows(lapply(split(sample_metric_df, sample_metric_df$ClusterID), function(d) {
      data.frame(ClusterID = unique(d$ClusterID)[1],
        detected_fisher_p = .fisher_flag(d, "detected"),
        diverse_fisher_p  = .fisher_flag(d, "diverse_present"),
        stringsAsFactors = FALSE)
    }))

    prevalence_compare <- prevalence_tests %>%
      dplyr::left_join(
        prevalence_by_group %>%
          dplyr::select(ClusterID, group, n_samples, n_detected, n_diverse_present,
                        detected_prevalence, diverse_prevalence) %>%
          tidyr::pivot_wider(names_from = group,
            values_from = c(n_samples, n_detected, n_diverse_present, detected_prevalence, diverse_prevalence)),
        by = "ClusterID"
      ) %>%
      dplyr::left_join(cluster_annot, by = "ClusterID")

    if (length(group_order) >= 2) {
      for (metric in c("detected_prevalence", "diverse_prevalence")) {
        c1 <- paste0(metric, "_", group_order[1])
        c2 <- paste0(metric, "_", group_order[2])
        if (all(c(c1, c2) %in% colnames(prevalence_compare))) {
          prevalence_compare[[paste0("delta_", metric)]] <- prevalence_compare[[c2]] - prevalence_compare[[c1]]
        }
      }
    }
    prevalence_compare$detected_fdr <- stats::p.adjust(prevalence_compare$detected_fisher_p, method = "fdr")
    prevalence_compare$diverse_fdr  <- stats::p.adjust(prevalence_compare$diverse_fisher_p,  method = "fdr")

    utils::write.csv(prevalence_compare,
      file.path(dir_base, sprintf("prevalence_compare_%s_%s.csv", weighting_tag, date_tag)),
      row.names = FALSE)

    metric_long <- sample_metric_df %>%
      dplyr::select(SampleID, ClusterID, dplyr::all_of(mainGroup),
                    dplyr::any_of(covariate), Species, Genus,
                    mean_pairwise_distance, abundance_weighted_pairwise_distance,
                    mean_nearest_neighbor_distance, within_cluster_entropy, within_cluster_evenness) %>%
      tidyr::pivot_longer(
        cols = c(mean_pairwise_distance, abundance_weighted_pairwise_distance,
                 mean_nearest_neighbor_distance, within_cluster_entropy, within_cluster_evenness),
        names_to = "metric", values_to = "value"
      ) %>%
      dplyr::filter(is.finite(value)) %>%
      dplyr::mutate(
        value_model = dplyr::case_when(
          grepl("distance", metric) ~ log10(value + 1e-6),
          TRUE ~ value
        )
      )

    .metric_delta <- function(d, labels = NULL) {
      d2 <- d[d[[mainGroup]] %in% group_order[seq_len(2)], , drop = FALSE]
      if (!is.null(labels)) d2[[mainGroup]] <- labels
      if (sum(d2[[mainGroup]] == group_order[1], na.rm = TRUE) < min_present_per_group) return(NA_real_)
      if (sum(d2[[mainGroup]] == group_order[2], na.rm = TRUE) < min_present_per_group) return(NA_real_)
      mean(d2$value_model[d2[[mainGroup]] == group_order[2]], na.rm = TRUE) -
        mean(d2$value_model[d2[[mainGroup]] == group_order[1]], na.rm = TRUE)
    }

    if (nrow(metric_long) > 0) {
      sample_metric_perm_boot <- metric_long %>%
        dplyr::group_by(ClusterID, metric) %>%
        dplyr::group_modify(~ {
          d <- .x
          obs <- .metric_delta(d)
          if (!is.finite(obs)) return(data.frame())
          pd <- replicate(n_perm, .metric_delta(d, labels = sample(d[[mainGroup]])))
          idx_1 <- which(d[[mainGroup]] == group_order[1])
          idx_2 <- which(d[[mainGroup]] == group_order[2])
          bd <- replicate(n_boot, {
            .metric_delta(d[c(sample(idx_1, length(idx_1), replace = TRUE),
                              sample(idx_2, length(idx_2), replace = TRUE)), ])
          })
          pd <- pd[is.finite(pd)]; bd <- bd[is.finite(bd)]
          data.frame(delta_observed = obs,
            permutation_n = length(pd),
            permutation_p_two_sided = ifelse(length(pd) > 0,
              (1 + sum(abs(pd) >= abs(obs), na.rm = TRUE)) / (length(pd) + 1), NA_real_),
            bootstrap_n = length(bd),
            bootstrap_ci_low  = ifelse(length(bd) > 0, stats::quantile(bd, 0.025, na.rm = TRUE), NA_real_),
            bootstrap_ci_high = ifelse(length(bd) > 0, stats::quantile(bd, 0.975, na.rm = TRUE), NA_real_),
            bootstrap_same_direction = ifelse(length(bd) > 0,
              mean(sign(bd) == sign(obs), na.rm = TRUE), NA_real_),
            n_group1 = sum(d[[mainGroup]] == group_order[1], na.rm = TRUE),
            n_group2 = sum(d[[mainGroup]] == group_order[2], na.rm = TRUE),
            stringsAsFactors = FALSE)
        }) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(metric) %>%
        dplyr::mutate(permutation_fdr_by_metric = stats::p.adjust(permutation_p_two_sided, method = "fdr")) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(cluster_annot, by = "ClusterID")

      utils::write.csv(sample_metric_perm_boot,
        file.path(dir_base, sprintf("sample_metric_permutation_bootstrap_%s_%s.csv", weighting_tag, date_tag)),
        row.names = FALSE)
    }
  }

  ###############################################
  # 9. Plots
  ###############################################
  plot_objs  <- list()
  plot_files <- list()

  # 9a. Forest plot
  if (nrow(res_case_present) > 0) {
    plot_colors <- c("FDR+, Positive" = "red", "FDR+, Negative" = "blue",
                     "Nominal p<0.05" = "grey40", "NS" = "grey80")
    p_forest <- ggplot2::ggplot(res_case_present,
        ggplot2::aes(x = estimate, y = Cluster_label, color = sig_cat)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = conf.low, xmax = conf.high),
                             width = 0.2, alpha = 0.7, orientation = "y") +
      ggplot2::scale_color_manual(values = plot_colors, drop = FALSE) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        axis.text.y          = ggplot2::element_text(size = 10, face = "italic"),
        plot.title           = ggplot2::element_text(size = 12, hjust = 0),
        plot.subtitle        = ggplot2::element_text(size = 10, hjust = 0),
        plot.title.position  = "plot"
      ) +
      ggplot2::labs(
        title    = "Sequence-cloud-compatible cluster pi association",
        subtitle = sprintf(
          "Present-only lm(log10(pi) ~ %s%s)\nNegative in %s, Positive in %s",
          mainGroup,
          ifelse(!is.null(covariate) && nzchar(covariate), paste0(" + ", covariate), ""),
          group_order[1], group_order[2]
        ),
        x = paste0("Estimate for ", mainGroup, group_order[2]),
        y = "Cluster", color = "Significance"
      )
    forest_pdf <- file.path(dir_base, sprintf("pi_forest_%s_%s.pdf", weighting_tag, date_tag))
    ggplot2::ggsave(forest_pdf, p_forest, width = forest_width,
      height = max(4, 0.3 * nrow(res_case_present) + 1))
    plot_objs$forest  <- p_forest
    plot_files$forest <- forest_pdf
    message("[Go_intoASV_graphs_support] Forest plot: ", forest_pdf)
  }

  # 9b. Support summary plot
  support_plot_df <- dplyr::bind_rows(
    if (nrow(pi_perm_boot)    > 0) pi_perm_boot    else data.frame(),
    if (nrow(graph_perm_boot) > 0) graph_perm_boot else data.frame()
  )

  if (nrow(support_plot_df) > 0 && "Species" %in% colnames(support_plot_df)) {
    support_plot_df <- support_plot_df %>%
      dplyr::filter(metric %in% c("pi_log10_present", "graph_dispersion",
                                  "graph_modularity", "graph_mean_pairwise_distance")) %>%
      dplyr::mutate(
        metric_label = dplyr::recode(metric,
          pi_log10_present             = "present-only log10(pi)",
          graph_dispersion             = "graph dispersion",
          graph_modularity             = "graph modularity",
          graph_mean_pairwise_distance = "graph mean pairwise distance"
        ),
        Cluster_label = ifelse(!is.na(Species) & nzchar(Species),
                               paste0(ClusterID, " (", Species, ")"), ClusterID),
        support = dplyr::case_when(
          metric == "pi_log10_present" & permutation_p_two_sided < 0.05 &
            bootstrap_same_direction >= 0.90 & delta_observed > 0 ~ "p<0.05 + stable positive",
          metric == "pi_log10_present" & permutation_p_two_sided < 0.05 &
            bootstrap_same_direction >= 0.90 & delta_observed < 0 ~ "p<0.05 + stable negative",
          bootstrap_same_direction >= 0.90 & delta_observed > 0 ~ "stable positive",
          bootstrap_same_direction >= 0.90 & delta_observed < 0 ~ "stable negative",
          TRUE ~ "weak"
        ),
        support = factor(support, levels = c("p<0.05 + stable positive", "p<0.05 + stable negative",
                                             "stable positive", "stable negative", "weak"))
      )

    support_cl <- if (!is.null(inset_layout)) {
      unique(unlist(inset_layout, use.names = FALSE))
    } else if (!is.null(target_clusters)) {
      as.character(target_clusters)
    } else {
      support_plot_df %>%
        dplyr::filter(is.finite(permutation_p_two_sided)) %>%
        dplyr::group_by(ClusterID) %>%
        dplyr::summarise(min_p = min(permutation_p_two_sided, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(min_p) %>%
        dplyr::slice_head(n = support_top_clusters) %>%
        dplyr::pull(ClusterID)
    }
    support_cl <- support_cl[support_cl %in% unique(support_plot_df$ClusterID)]

    sdf_f <- support_plot_df %>%
      dplyr::filter(ClusterID %in% support_cl) %>%
      dplyr::mutate(Cluster_label = forcats::fct_reorder(Cluster_label, delta_observed, .fun = mean))

    if (nrow(sdf_f) > 0) {
      p_support <- ggplot2::ggplot(sdf_f,
          ggplot2::aes(x = delta_observed, y = Cluster_label, color = support)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = bootstrap_ci_low, xmax = bootstrap_ci_high),
                               width = 0.15, alpha = 0.55, orientation = "y") +
        ggplot2::geom_point(ggplot2::aes(fill = support), shape = 21, size = 2.7, stroke = 0.75) +
        ggplot2::facet_wrap(~ metric_label, scales = "free_x", ncol = 2) +
        ggplot2::scale_color_manual(values = c(
          "p<0.05 + stable positive" = "red",  "p<0.05 + stable negative" = "blue",
          "stable positive"          = "red",  "stable negative"          = "blue",
          "weak"                     = "grey70"), drop = FALSE) +
        ggplot2::scale_fill_manual(values = c(
          "p<0.05 + stable positive" = "red",   "p<0.05 + stable negative" = "blue",
          "stable positive"          = "white", "stable negative"          = "white",
          "weak"                     = "grey70"), drop = FALSE) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(
          axis.text.y      = ggplot2::element_text(size = 7, face = "italic"),
          legend.position  = "bottom",
          strip.text       = ggplot2::element_text(face = "bold")
        ) +
        ggplot2::labs(
          title    = paste0("Sequence-cloud support summary: ", weighting_tag),
          subtitle = "Delta = case - control; error bars = bootstrap 95% CI\nFilled = pi permutation p<0.05 + stable; open = stable direction/topology",
          x = "Delta", y = "Cluster", color = "Support", fill = "Support"
        )
      support_pdf <- file.path(dir_base, sprintf("support_summary_%s_%s.pdf", weighting_tag, date_tag))
      ggplot2::ggsave(support_pdf, p_support, width = support_width,
        height = max(6, 0.23 * length(unique(sdf_f$Cluster_label)) + 3))
      plot_objs$support  <- p_support
      plot_files$support <- support_pdf
      message("[Go_intoASV_graphs_support] Support plot: ", support_pdf)
    }
  }

  # 9c. Integrated evidence
  evidence_list <- list()
  if (nrow(prevalence_compare) > 0) {
    for (spec in list(
      list("detected_prevalence", "cluster prevalence",  "detected_fisher_p", "detected_fdr"),
      list("diverse_prevalence",  "ASV>1 prevalence",    "diverse_fisher_p",  "diverse_fdr")
    )) {
      d_col <- paste0("delta_", spec[[1]])
      need  <- c("ClusterID", "Species", d_col, spec[[3]], spec[[4]])
      if (all(need %in% colnames(prevalence_compare))) {
        evidence_list[[spec[[1]]]] <- prevalence_compare %>%
          dplyr::transmute(ClusterID, Species,
            metric = spec[[1]], metric_label = spec[[2]],
            delta_observed = .data[[d_col]],
            permutation_p_two_sided = .data[[spec[[3]]]],
            fdr = .data[[spec[[4]]]],
            bootstrap_ci_low = NA_real_, bootstrap_ci_high = NA_real_, bootstrap_same_direction = NA_real_)
      }
    }
  }
  if (nrow(pi_perm_boot) > 0) {
    evidence_list$pi <- pi_perm_boot %>%
      dplyr::transmute(ClusterID, Species,
        metric = "pi_log10_present", metric_label = "present-only log10(pi)",
        delta_observed, permutation_p_two_sided, fdr = NA_real_,
        bootstrap_ci_low, bootstrap_ci_high, bootstrap_same_direction)
  }
  if (nrow(sample_metric_perm_boot) > 0) {
    metric_labels_map <- c(
      mean_pairwise_distance                = "sample mean pairwise distance",
      abundance_weighted_pairwise_distance  = "sample abundance-weighted distance",
      mean_nearest_neighbor_distance        = "sample nearest-neighbor distance",
      within_cluster_entropy                = "within-cluster ASV entropy",
      within_cluster_evenness               = "within-cluster ASV evenness"
    )
    evidence_list$sample_metrics <- sample_metric_perm_boot %>%
      dplyr::filter(metric %in% names(metric_labels_map)) %>%
      dplyr::mutate(metric_label = unname(metric_labels_map[metric])) %>%
      dplyr::transmute(ClusterID, Species, metric, metric_label, delta_observed,
        permutation_p_two_sided, fdr = permutation_fdr_by_metric,
        bootstrap_ci_low, bootstrap_ci_high, bootstrap_same_direction)
  }

  evidence_df <- dplyr::bind_rows(evidence_list)

  if (nrow(evidence_df) > 0) {
    ev_cl <- if (!is.null(inset_layout)) unique(unlist(inset_layout, use.names = FALSE)) else {
      if (!is.null(target_clusters)) as.character(target_clusters) else character(0)
    }
    ev_cl <- ev_cl[!is.na(ev_cl) & nzchar(ev_cl)]
    if (length(ev_cl) > 0) evidence_df <- evidence_df[evidence_df$ClusterID %in% ev_cl, , drop = FALSE]
  }

  if (nrow(evidence_df) > 0) {
    evidence_df <- evidence_df %>%
      dplyr::mutate(
        Cluster_label = ifelse(!is.na(Species) & nzchar(Species),
                               paste0(ClusterID, " (", Species, ")"), ClusterID),
        support = dplyr::case_when(
          permutation_p_two_sided < 0.05 &
            (is.na(bootstrap_same_direction) | bootstrap_same_direction >= 0.90) &
            delta_observed > 0 ~ "p<0.05 + stable positive",
          permutation_p_two_sided < 0.05 &
            (is.na(bootstrap_same_direction) | bootstrap_same_direction >= 0.90) &
            delta_observed < 0 ~ "p<0.05 + stable negative",
          bootstrap_same_direction >= 0.90 & delta_observed > 0 ~ "stable positive",
          bootstrap_same_direction >= 0.90 & delta_observed < 0 ~ "stable negative",
          TRUE ~ "weak"
        ),
        support = factor(support, levels = c("p<0.05 + stable positive", "p<0.05 + stable negative",
                                             "stable positive", "stable negative", "weak")),
        Cluster_label = forcats::fct_reorder(Cluster_label, delta_observed, .fun = mean)
      )

    utils::write.csv(evidence_df,
      file.path(dir_base, sprintf("integrated_evidence_%s_%s.csv", weighting_tag, date_tag)),
      row.names = FALSE)

    p_evidence <- ggplot2::ggplot(evidence_df,
        ggplot2::aes(x = delta_observed, y = Cluster_label, color = support)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = bootstrap_ci_low, xmax = bootstrap_ci_high),
                             width = 0.15, alpha = 0.55, orientation = "y", na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(fill = support), shape = 21, size = 2.7, stroke = 0.75) +
      ggplot2::facet_wrap(~ metric_label, scales = "free_x", ncol = 2) +
      ggplot2::scale_color_manual(values = c(
        "p<0.05 + stable positive" = "red",  "p<0.05 + stable negative" = "blue",
        "stable positive"          = "red",  "stable negative"          = "blue",
        "weak"                     = "grey70"), drop = FALSE) +
      ggplot2::scale_fill_manual(values = c(
        "p<0.05 + stable positive" = "red",   "p<0.05 + stable negative" = "blue",
        "stable positive"          = "white", "stable negative"          = "white",
        "weak"                     = "grey70"), drop = FALSE) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        axis.text.y     = ggplot2::element_text(size = 7, face = "italic"),
        legend.position = "bottom",
        strip.text      = ggplot2::element_text(face = "bold")
      ) +
      ggplot2::labs(
        title    = paste0("Sequence-cloud integrated evidence: ", weighting_tag),
        subtitle = "Sample-level evidence: prevalence, present-only pi, sequence-distance metrics\nDelta = case - control; error bars = bootstrap 95% CI where available",
        x = "Delta", y = "Cluster", color = "Support", fill = "Support"
      )

    evidence_pdf <- file.path(dir_base, sprintf("integrated_evidence_%s_%s.pdf", weighting_tag, date_tag))
    ggplot2::ggsave(evidence_pdf, p_evidence, width = evidence_width,
      height = max(7, 0.25 * length(unique(evidence_df$Cluster_label)) + 5))
    plot_objs$evidence  <- p_evidence
    plot_files$evidence <- evidence_pdf
    message("[Go_intoASV_graphs_support] Integrated evidence plot: ", evidence_pdf)
  }

  message("[Go_intoASV_graphs_support] Done. Output: ", dir_base)

  ###############################################
  # 10. Return
  ###############################################
  invisible(list(
    pi_summary                          = pi_present_summary,
    pi_compare                          = pi_present_compare,
    pi_lm                               = res_lm_present,
    pi_case_control                     = res_case_present,
    pi_permutation_bootstrap            = pi_perm_boot,
    graph_permutation_bootstrap         = graph_perm_boot,
    sample_metrics                      = sample_metric_df,
    prevalence                          = prevalence_compare,
    sample_metric_permutation_bootstrap = sample_metric_perm_boot,
    integrated_evidence                 = evidence_df,
    plot                                = plot_objs,
    plot_file                           = plot_files,
    dir_base                            = dir_base
  ))
}
