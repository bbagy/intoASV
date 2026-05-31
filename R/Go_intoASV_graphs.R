#' intoASV sequence-cloud graph visualization
#'
#' Build sparse ASV sequence-cloud graphs from a \code{Go_intoASV()} result and
#' draw a copangraph-style overview of intra-cluster diversification topology.
#'
#' @param psIN A \code{phyloseq} object returned by \code{Go_intoASV()}.
#' @param project Character; project prefix for output folders.
#' @param cluster_map Cluster map data frame or CSV path. Must contain
#'   \code{ASV} and \code{ClusterID}; \code{Species} is used when available.
#' @param clusters Optional vector of cluster IDs to plot. If \code{NULL},
#'   clusters are selected from \code{target_species_patterns} or by size.
#' @param target_species_patterns Optional species names or regular expressions
#'   used to select biologically relevant clusters by species annotation. Simple
#'   names such as \code{"Lactobacillus iners"} are internally converted to
#'   flexible patterns that tolerate spaces, dots, underscores, semicolons, and
#'   hyphens between genus and species.
#' @param clusters_per_species Integer; maximum clusters retained per matched
#'   species pattern.
#' @param max_clusters Integer; maximum clusters retained when no species
#'   pattern is supplied.
#' @param inset_layout Optional list assigning inset clusters to zoom columns.
#'   Example: \code{list(zoom1 = c("Cluster_1","Cluster_3"),
#'   zoom2 = c("Cluster_2","Cluster_4"))}.
#' @param mainGroup Optional sample metadata variable used to draw up to four
#'   subgroup sequence-cloud panels in one 2 x 2 comparison image. If
#'   \code{NULL}, the standard single overview plot is drawn.
#' @param order Optional character vector controlling the plotting order of
#'   \code{mainGroup} levels. Values not present in the data are ignored.
#' @param name Optional output-name suffix, similar to Gotools naming style.
#' @param width,height Numeric PDF size in inches.
#'
#' @return A list containing selected clusters, graph objects, sequence-cloud
#'   objects, graph statistics, the ggplot object, and output file paths.
#'
#' @details
#' This is not a phylogenetic tree, assembly graph, or co-occurrence network.
#' Each node is an ASV. Edges are sparse sequence-neighborhood links based on a
#' minimum spanning tree plus k-nearest neighbors. The plot is designed to expose
#' intra-taxonomic sequence diversification as a sequence-space ecological cloud.
#'
#' @examples
#' \dontrun{
#' ps_intoASV <- Go_intoASV(psIN = ps, project = "Study", global_similarity_cutoff = 0.97)
#' graph_res <- Go_intoASV_graphs(psIN = ps_intoASV, cluster_map = cluster_map_i, project = "Study")
#' }
#'
#' @export
Go_intoASV_graphs <- function(
    psIN,
    cluster_map,
    project = "intoASV_graphs",
    clusters = NULL,
    target_species_patterns = NULL,
    clusters_per_species = 2,
    max_clusters = 20,
    inset_layout = NULL,
    mainGroup = NULL,
    order = NULL,
    name = NULL,
    width = 12,
    height = 8
) {

  if (!inherits(psIN, "phyloseq")) {
    stop("[Go_intoASV_graphs] psIN must be a phyloseq object.")
  }
  settings <- list(
    graph_distance_cutoff = 0.030,
    knn = 3,
    axis_limit = 1.08,
    circle_radius = 1.12,
    show_circle = FALSE,
    cluster_center_radius = 0.74,
    cluster_radius_range = c(0.48, 1.00),
    cluster_distance_transform = "log1p",
    cluster_distance_alpha = 25,
    cluster_scale_range = c(0.30, 0.70),
    node_size_range = c(1.50, 7.80),
    edge_width = 0.18,
    node_stroke = 0.04,
    edge_alpha_range = c(0.12, 0.55),
    inset_node_size_range = c(1.60, 6.50),
    inset_edge_width_range = c(0.06, 0.42),
    seed = 123
  )

  species_palette <- NULL
  inset_clusters <- NULL
  graph_distance_cutoff <- settings$graph_distance_cutoff
  knn <- settings$knn
  axis_limit <- settings$axis_limit
  circle_radius <- settings$circle_radius
  show_circle <- settings$show_circle
  cluster_center_radius <- settings$cluster_center_radius
  cluster_radius_range <- settings$cluster_radius_range
  cluster_distance_transform <- settings$cluster_distance_transform
  cluster_distance_alpha <- settings$cluster_distance_alpha
  cluster_scale_range <- settings$cluster_scale_range
  node_size_range <- settings$node_size_range
  edge_width <- settings$edge_width
  node_stroke <- settings$node_stroke
  edge_alpha_range <- settings$edge_alpha_range
  inset_node_size_range <- settings$inset_node_size_range
  inset_edge_width_range <- settings$inset_edge_width_range
  seed <- settings$seed

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  ###############################################
  # 0. Output directory structure
  ###############################################
  date_tag <- format(Sys.Date(), "%y%m%d")
  dir_base <- sprintf("%s_%s/intoASV/sequence_cloud", project, date_tag)
  dir.create(dir_base, recursive = TRUE, showWarnings = FALSE)

  .rescale <- function(x, to, from = range(x, na.rm = TRUE)) {
    if (length(x) == 0) return(numeric())
    if (!all(is.finite(from)) || abs(diff(from)) < .Machine$double.eps) {
      return(rep(mean(to), length(x)))
    }
    (x - from[1]) / diff(from) * diff(to) + to[1]
  }

  .transform_distance <- function(d) {
    if (identical(cluster_distance_transform, "log1p")) {
      return(log1p(cluster_distance_alpha * d) / log1p(cluster_distance_alpha))
    }
    if (identical(cluster_distance_transform, "sqrt")) {
      return(sqrt(d))
    }
    d
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
    if (is.null(seqs_tmp)) {
      stop("[Go_intoASV_graphs] No valid DNA sequences found in refseq, taxa_names, or OTU rownames.")
    }
    names(seqs_tmp) <- phyloseq::taxa_names(ps_obj)
    seqs_tmp
  }

  .get_otu_matrix <- function(ps_obj) {
    otu_tab <- as(phyloseq::otu_table(ps_obj), "matrix")
    if (!phyloseq::taxa_are_rows(ps_obj)) otu_tab <- t(otu_tab)
    otu_tab
  }

  .get_tax_table <- function(ps_obj) {
    as.data.frame(as(phyloseq::tax_table(ps_obj), "matrix"), stringsAsFactors = FALSE)
  }

  .safe_file_token <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x[1])) return(NA_character_)
    x <- trimws(as.character(x[1]))
    if (!nzchar(x)) return(NA_character_)
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    if (!nzchar(x)) NA_character_ else x
  }

  .extract_weighting <- function(x) {
    x <- tolower(paste(x, collapse = " "))
    hit <- regmatches(x, regexpr("(abundance|entropy)", x, perl = TRUE))
    if (length(hit) == 0 || hit == "") "weightingNA" else hit
  }

  .extract_similarity <- function(x) {
    x <- paste(x, collapse = " ")
    hit <- regmatches(x, regexpr("similarity[_-]?([0-9]{3,4}|[0-9]+\\.[0-9]+)", x, perl = TRUE))
    if (length(hit) == 0 || hit == "") {
      return("similarityNA")
    }
    val <- sub(".*similarity[_-]?", "", hit)
    if (grepl("\\.", val)) val <- gsub("\\.", "", sprintf("%.3f", as.numeric(val)))
    paste0("similarity", val)
  }

  .escape_regex <- function(x) {
    gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
  }

  .looks_like_regex <- function(x) {
    grepl("\\\\|\\[|\\]|\\(|\\)|\\||\\+|\\*|\\^|\\$", x)
  }

  .clean_species_label <- function(x) {
    x <- gsub("\\\\s\\+", " ", x)
    x <- gsub("\\\\\\.", ".", x)
    x <- gsub("\\$", "", x)
    x <- gsub("\\s+", " ", x)
    trimws(x)
  }

  .species_query_to_regex <- function(x) {
    x <- trimws(x)
    if (!nzchar(x)) return(NA_character_)
    if (.looks_like_regex(x)) return(x)

    tokens <- unlist(strsplit(x, "[[:space:]_.;|:/-]+"))
    tokens <- tokens[nzchar(tokens)]
    if (length(tokens) == 0) return(NA_character_)

    token_patterns <- vapply(tokens, function(tok) {
      if (tolower(tok) == "sp") return("sp\\.?")
      .escape_regex(tok)
    }, character(1))
    sep <- "[[:space:]_.;|:/-]+"
    paste0("(^|[^[:alnum:]])", paste(token_patterns, collapse = sep), "($|[^[:alnum:]])")
  }

  .prepare_species_queries <- function(x) {
    if (is.null(x)) return(NULL)
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(trimws(x))]
    if (length(x) == 0) return(NULL)
    data.frame(
      input = x,
      label = vapply(x, .clean_species_label, character(1)),
      pattern = vapply(x, .species_query_to_regex, character(1)),
      stringsAsFactors = FALSE
    )
  }

  target_species_queries <- .prepare_species_queries(target_species_patterns)

  cluster_map_file <- NA_character_
  if (is.character(cluster_map) && length(cluster_map) == 1) {
    cluster_map_file <- cluster_map
    if (!file.exists(cluster_map_file)) {
      stop("[Go_intoASV_graphs] cluster_map file does not exist: ", cluster_map_file)
    }
    cluster_map <- utils::read.csv(cluster_map_file, check.names = FALSE)
  }

  if (!all(c("ASV", "ClusterID") %in% colnames(cluster_map))) {
    stop("[Go_intoASV_graphs] cluster_map must contain ASV and ClusterID columns.")
  }
  if (!"Species" %in% colnames(cluster_map)) cluster_map$Species <- NA_character_

  .select_target_clusters <- function(cm) {
    cm$n_key <- paste(cm$ClusterID, cm$Species, sep = "||")
    count_tab <- sort(table(cm$n_key), decreasing = TRUE)
    target_df <- do.call(rbind, strsplit(names(count_tab), "\\|\\|", fixed = FALSE))
    target_df <- data.frame(
      ClusterID = target_df[, 1],
      Species = target_df[, 2],
      n = as.integer(count_tab),
      stringsAsFactors = FALSE
    )

    if (!is.null(clusters)) {
      return(target_df[target_df$ClusterID %in% clusters, , drop = FALSE])
    }

    if (!is.null(target_species_queries)) {
      target_df$target_species <- vapply(target_df$Species, function(x) {
        if (is.na(x)) return(NA_character_)
        hit <- which(vapply(target_species_queries$pattern, grepl, logical(1), x = x, ignore.case = TRUE))
        ifelse(length(hit) == 0, NA_character_, target_species_queries$label[hit[1]])
      }, character(1))
      target_df <- target_df[!is.na(target_df$target_species), , drop = FALSE]
      target_df <- target_df[order(target_df$target_species, -target_df$n), , drop = FALSE]
      keep <- unlist(tapply(seq_len(nrow(target_df)), target_df$target_species, head, n = clusters_per_species))
      return(target_df[keep, , drop = FALSE])
    }

    head(target_df[order(-target_df$n), , drop = FALSE], max_clusters)
  }

  target_clusters <- .select_target_clusters(cluster_map)
  if (nrow(target_clusters) == 0) {
    stop("[Go_intoASV_graphs] No target clusters selected.")
  }

  .build_intoasv_graph <- function(cluster_id, ps_obj = psIN) {
    otu_tab <- .get_otu_matrix(ps_obj)
    tax_tab <- .get_tax_table(ps_obj)
    seqs_all <- .get_asv_sequences(ps_obj)
    cm_sub <- cluster_map[cluster_map$ClusterID == cluster_id, , drop = FALSE]
    asv_ids <- intersect(cm_sub$ASV, rownames(otu_tab))
    asv_ids <- intersect(asv_ids, names(seqs_all))
    if (length(asv_ids) < 2) return(NULL)

    seqs <- Biostrings::DNAStringSet(seqs_all[asv_ids])
    dm <- as.matrix(DECIPHER::DistanceMatrix(seqs, includeTerminalGaps = FALSE))
    dm <- dm[asv_ids, asv_ids, drop = FALSE]

    edge_keep <- which(upper.tri(dm) & dm <= graph_distance_cutoff, arr.ind = TRUE)
    edge_df <- data.frame()
    if (nrow(edge_keep) > 0) {
      edge_df <- data.frame(
        from = rownames(dm)[edge_keep[, 1]],
        to = colnames(dm)[edge_keep[, 2]],
        distance = dm[edge_keep],
        similarity = 1 - dm[edge_keep],
        edge_type = "threshold",
        stringsAsFactors = FALSE
      )
    }

    if (!is.null(knn) && knn > 0) {
      knn_edges <- lapply(rownames(dm), function(i) {
        x <- dm[i, ]
        x[i] <- Inf
        j <- names(sort(x, decreasing = FALSE))[seq_len(min(knn, length(x) - 1))]
        data.frame(
          from = i, to = j, distance = as.numeric(x[j]),
          similarity = 1 - as.numeric(x[j]), edge_type = "knn",
          stringsAsFactors = FALSE
        )
      })
      edge_df <- rbind(edge_df, do.call(rbind, knn_edges))
      pair <- ifelse(edge_df$from < edge_df$to,
                     paste(edge_df$from, edge_df$to, sep = "__"),
                     paste(edge_df$to, edge_df$from, sep = "__"))
      edge_df <- edge_df[order(pair, edge_df$distance), , drop = FALSE]
      edge_df <- edge_df[!duplicated(pair[order(pair, edge_df$distance)]), , drop = FALSE]
    }

    node_abund <- rowSums(otu_tab[asv_ids, , drop = FALSE])
    node_df <- data.frame(
      name = asv_ids,
      sequence = seqs_all[asv_ids],
      total_abundance = as.numeric(node_abund),
      prevalence = rowSums(otu_tab[asv_ids, , drop = FALSE] > 0),
      stringsAsFactors = FALSE
    )

    tax_cols <- intersect(c("Phylum", "Class", "Order", "Family", "Genus", "Species"), colnames(tax_tab))
    if (length(tax_cols) > 0) {
      node_df <- cbind(node_df, tax_tab[asv_ids, tax_cols, drop = FALSE])
    }

    graph <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = node_df)
    list(graph = graph, nodes = node_df, edges = edge_df, distance = dm, cluster_id = cluster_id)
  }

  .build_sequence_cloud_graph <- function(graph_obj) {
    dm <- graph_obj$distance
    nodes <- graph_obj$nodes
    ids <- intersect(rownames(dm), nodes$name)
    if (length(ids) < 2) return(NULL)
    dm <- dm[ids, ids, drop = FALSE]
    nodes <- nodes[match(ids, nodes$name), , drop = FALSE]

    full_g <- igraph::graph_from_adjacency_matrix(dm, mode = "undirected", weighted = TRUE, diag = FALSE)
    mst_g <- igraph::mst(full_g, weights = igraph::E(full_g)$weight)
    mst_edges <- as.data.frame(igraph::as_edgelist(mst_g), stringsAsFactors = FALSE)
    if (ncol(mst_edges) == 2) {
      colnames(mst_edges) <- c("from", "to")
      mst_edges$edge_type <- "mst"
    } else {
      mst_edges <- data.frame(from = character(), to = character(), edge_type = character())
    }

    knn_edges <- do.call(rbind, lapply(rownames(dm), function(i) {
      x <- dm[i, ]
      x[i] <- Inf
      j <- names(sort(x, decreasing = FALSE))[seq_len(min(knn, length(x) - 1))]
      data.frame(from = i, to = j, edge_type = "knn", stringsAsFactors = FALSE)
    }))

    edge_df <- rbind(mst_edges, knn_edges)
    pair <- ifelse(edge_df$from < edge_df$to,
                   paste(edge_df$from, edge_df$to, sep = "__"),
                   paste(edge_df$to, edge_df$from, sep = "__"))
    edge_df <- edge_df[!duplicated(pair), , drop = FALSE]
    edge_df$distance <- mapply(function(a, b) dm[a, b], edge_df$from, edge_df$to)
    edge_df$similarity <- 1 - edge_df$distance
    edge_df$inv_distance <- 1 / (edge_df$distance + 1e-6)

    cloud_g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = nodes)
    if (igraph::gsize(cloud_g) > 0) {
      igraph::E(cloud_g)$distance <- edge_df$distance
      igraph::E(cloud_g)$similarity <- edge_df$similarity
      igraph::E(cloud_g)$inv_distance <- edge_df$inv_distance
      igraph::E(cloud_g)$edge_type <- edge_df$edge_type
    }

    set.seed(seed)
    comm <- tryCatch(
      igraph::cluster_louvain(cloud_g, weights = igraph::E(cloud_g)$inv_distance),
      error = function(e) NULL
    )
    if (!is.null(comm)) {
      igraph::V(cloud_g)$cloud_module <- as.factor(igraph::membership(comm))
    } else {
      igraph::V(cloud_g)$cloud_module <- factor(rep(1, igraph::vcount(cloud_g)))
    }
    nodes$cloud_module <- as.character(igraph::V(cloud_g)$cloud_module[match(nodes$name, igraph::V(cloud_g)$name)])

    set.seed(seed)
    layout <- igraph::layout_with_fr(cloud_g, weights = igraph::E(cloud_g)$inv_distance, niter = 1000)
    colnames(layout) <- c("x", "y")
    layout <- as.data.frame(layout)
    layout$name <- igraph::V(cloud_g)$name

    list(
      graph = cloud_g,
      nodes = nodes,
      edges = edge_df,
      distance = dm,
      layout = layout,
      cluster_id = graph_obj$cluster_id
    )
  }

  .sequence_cloud_stats <- function(cloud_obj) {
    g <- cloud_obj$graph
    dm <- cloud_obj$distance
    layout <- cloud_obj$layout
    lower_dm <- dm[lower.tri(dm)]
    cx <- mean(layout$x)
    cy <- mean(layout$y)
    radial <- sqrt((layout$x - cx)^2 + (layout$y - cy)^2)
    abund <- igraph::V(g)$total_abundance
    abund[!is.finite(abund)] <- 0
    if (sum(abund) <= 0) abund <- rep(1, length(radial))
    set.seed(seed)
    comm <- tryCatch(igraph::cluster_louvain(g, weights = igraph::E(g)$inv_distance), error = function(e) NULL)
    modularity <- if (!is.null(comm)) igraph::modularity(comm) else NA_real_
    species <- if ("Species" %in% colnames(cloud_obj$nodes)) {
      names(sort(table(cloud_obj$nodes$Species), decreasing = TRUE))[1]
    } else {
      NA_character_
    }
    genus <- if ("Genus" %in% colnames(cloud_obj$nodes)) {
      names(sort(table(cloud_obj$nodes$Genus), decreasing = TRUE))[1]
    } else {
      sub("\\s+.*$", "", species)
    }
    data.frame(
      ClusterID = cloud_obj$cluster_id,
      Species = species,
      Genus = genus,
      n_nodes = igraph::vcount(g),
      n_edges = igraph::gsize(g),
      components = igraph::components(g)$no,
      modularity = modularity,
      diameter = suppressWarnings(igraph::diameter(g, directed = FALSE, weights = NA)),
      clustering_coefficient = igraph::transitivity(g, type = "global", isolates = "zero"),
      degree_centralization = igraph::centr_degree(g, normalized = TRUE)$centralization,
      cloud_dispersion = mean(radial, na.rm = TRUE),
      abundance_weighted_dispersion = stats::weighted.mean(radial, abund, na.rm = TRUE),
      mean_pairwise_distance = mean(lower_dm, na.rm = TRUE),
      max_pairwise_distance = max(lower_dm, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  .build_cloud_set <- function(ps_obj) {
    graph_list_i <- list()
    for (cluster_i in unique(target_clusters$ClusterID)) {
      graph_i <- .build_intoasv_graph(cluster_i, ps_obj = ps_obj)
      if (!is.null(graph_i)) graph_list_i[[cluster_i]] <- graph_i
    }
    if (length(graph_list_i) == 0) return(NULL)

    sequence_cloud_list_i <- list()
    stats_list_i <- list()
    for (cluster_i in names(graph_list_i)) {
      cloud_i <- .build_sequence_cloud_graph(graph_list_i[[cluster_i]])
      if (is.null(cloud_i)) next
      sequence_cloud_list_i[[cluster_i]] <- cloud_i
      stats_list_i[[cluster_i]] <- .sequence_cloud_stats(cloud_i)
    }
    if (length(sequence_cloud_list_i) == 0) return(NULL)
    list(
      graph_list = graph_list_i,
      sequence_cloud_list = sequence_cloud_list_i,
      stats_df = do.call(rbind, stats_list_i)
    )
  }

  cloud_set <- .build_cloud_set(psIN)
  if (is.null(cloud_set)) stop("[Go_intoASV_graphs] No sequence-cloud graph could be built.")
  graph_list <- cloud_set$graph_list
  sequence_cloud_list <- cloud_set$sequence_cloud_list
  stats_df <- cloud_set$stats_df

  utils::write.csv(target_clusters, file.path(dir_base, sprintf("target_clusters_%s.csv", date_tag)), row.names = FALSE)
  utils::write.csv(stats_df, file.path(dir_base, sprintf("sequence_cloud_stats_%s.csv", date_tag)), row.names = FALSE)

  .build_genus_palette <- function(stats_df_i) {
    pal_df <- unique(stats_df_i[, c("Species", "Genus"), drop = FALSE])
    pal_df <- pal_df[!is.na(pal_df$Species) & nzchar(pal_df$Species), , drop = FALSE]
    pal_df$Genus[is.na(pal_df$Genus) | !nzchar(pal_df$Genus)] <- sub("\\s+.*$", "", pal_df$Species[is.na(pal_df$Genus) | !nzchar(pal_df$Genus)])
    pal_df <- pal_df[order(pal_df$Genus, pal_df$Species), , drop = FALSE]
    genera <- sort(unique(pal_df$Genus))
    high_contrast <- c(
      "#3B0F70", "#D55E00", "#0072B2", "#009E73", "#CC79A7",
      "#E69F00", "#56B4E9", "#8C564B", "#BDBD00", "#6A3D9A",
      "#E31A1C", "#1B9E77", "#7570B3", "#A6761D", "#666666"
    )
    if (length(genera) <= length(high_contrast)) {
      base_cols <- high_contrast[seq_along(genera)]
    } else {
      base_cols <- c(
        high_contrast,
        grDevices::hcl.colors(length(genera) - length(high_contrast), "Dynamic")
      )
    }
    genus_base <- stats::setNames(base_cols, genera)
    species_cols <- unlist(lapply(genera, function(genus_i) {
      species_i <- sort(unique(pal_df$Species[pal_df$Genus == genus_i]))
      n_i <- length(species_i)
      base_i <- grDevices::col2rgb(genus_base[[genus_i]])[, 1] / 255
      mix_vals <- if (n_i == 1) 0 else seq(0, 0.42, length.out = n_i)
      cols_i <- vapply(seq_along(species_i), function(j) {
        rgb_i <- base_i * (1 - mix_vals[j]) + c(1, 1, 1) * mix_vals[j]
        grDevices::rgb(rgb_i[1], rgb_i[2], rgb_i[3])
      }, character(1))
      stats::setNames(cols_i, species_i)
    }))
    species_cols
  }

  .compute_centers <- function(sequence_cloud_list_i, cluster_layout_i) {
    cloud_names_i <- names(sequence_cloud_list_i)
    cloud_stats_i <- do.call(rbind, lapply(sequence_cloud_list_i, .sequence_cloud_stats))
    cloud_stats_i <- cloud_stats_i[match(cloud_names_i, cloud_stats_i$ClusterID), , drop = FALSE]
    angles_i <- seq(pi / 2, pi / 2 + 2 * pi, length.out = length(cloud_names_i) + 1)[-length(cloud_names_i) - 1]

    if (identical(cluster_layout_i, "transformed_mds") && length(cloud_names_i) >= 3) {
      all_node_df_i <- do.call(rbind, lapply(cloud_names_i, function(cl) {
        nodes <- sequence_cloud_list_i[[cl]]$nodes
        data.frame(
          ClusterID = cl,
          name = nodes$name,
          sequence = nodes$sequence,
          stringsAsFactors = FALSE
        )
      }))
      all_node_df_i <- all_node_df_i[!is.na(all_node_df_i$sequence) & nchar(all_node_df_i$sequence) > 0, , drop = FALSE]
      seq_vec_i <- stats::setNames(all_node_df_i$sequence, paste(all_node_df_i$ClusterID, all_node_df_i$name, sep = "__"))
      dna_i <- Biostrings::DNAStringSet(seq_vec_i)
      all_dm_i <- as.matrix(DECIPHER::DistanceMatrix(dna_i, includeTerminalGaps = FALSE))
      cl_dist_i <- matrix(0, nrow = length(cloud_names_i), ncol = length(cloud_names_i), dimnames = list(cloud_names_i, cloud_names_i))
      for (i in seq_along(cloud_names_i)) {
        for (j in seq_along(cloud_names_i)) {
          if (i >= j) next
          keys_i <- names(seq_vec_i)[all_node_df_i$ClusterID == cloud_names_i[i]]
          keys_j <- names(seq_vec_i)[all_node_df_i$ClusterID == cloud_names_i[j]]
          d_ij <- mean(all_dm_i[keys_i, keys_j, drop = FALSE], na.rm = TRUE)
          cl_dist_i[i, j] <- d_ij
          cl_dist_i[j, i] <- d_ij
        }
      }
      cl_dist_t_i <- .transform_distance(cl_dist_i)
      diag(cl_dist_t_i) <- 0
      mds_i <- stats::cmdscale(stats::as.dist(cl_dist_t_i), k = 2)
      if (is.null(mds_i) || any(!is.finite(mds_i))) {
        centers_i <- data.frame(ClusterID = cloud_names_i, cx = cluster_center_radius * cos(angles_i), cy = cluster_center_radius * sin(angles_i))
      } else {
        centers_i <- data.frame(ClusterID = rownames(mds_i), cx = mds_i[, 1], cy = mds_i[, 2], stringsAsFactors = FALSE)
        centers_i$cx <- centers_i$cx - mean(centers_i$cx, na.rm = TRUE)
        centers_i$cy <- centers_i$cy - mean(centers_i$cy, na.rm = TRUE)
        max_r_i <- max(sqrt(centers_i$cx^2 + centers_i$cy^2), na.rm = TRUE)
        if (is.finite(max_r_i) && max_r_i > 0) {
          centers_i$cx <- centers_i$cx / max_r_i * cluster_center_radius
          centers_i$cy <- centers_i$cy / max_r_i * cluster_center_radius
        }
        centers_i <- centers_i[match(cloud_names_i, centers_i$ClusterID), , drop = FALSE]
      }
    } else if (identical(cluster_layout_i, "equal_circle")) {
      centers_i <- data.frame(ClusterID = cloud_names_i, cx = cluster_center_radius * cos(angles_i), cy = cluster_center_radius * sin(angles_i))
    } else {
      radius_i <- .rescale(cloud_stats_i$mean_pairwise_distance, to = cluster_radius_range)
      centers_i <- data.frame(ClusterID = cloud_names_i, cx = radius_i * cos(angles_i), cy = radius_i * sin(angles_i))
    }

    list(
      centers = centers_i,
      cloud_names = cloud_names_i,
      cloud_stats = cloud_stats_i
    )
  }

  plot_files <- list()
  plot_objs <- list()
  main_group_levels <- NULL
  if (!is.null(mainGroup)) {
    sample_df_for_levels <- data.frame(phyloseq::sample_data(psIN))
    if (!mainGroup %in% colnames(sample_df_for_levels)) {
      stop("[Go_intoASV_graphs] mainGroup was not found in sample_data: ", mainGroup)
    }
    main_group_levels <- unique(as.character(sample_df_for_levels[[mainGroup]]))
    main_group_levels <- main_group_levels[!is.na(main_group_levels) & nzchar(main_group_levels)]
    if (!is.null(order)) {
      main_group_levels <- c(
        intersect(as.character(order), main_group_levels),
        setdiff(main_group_levels, as.character(order))
      )
    }
    if (length(main_group_levels) > 4) {
      warning("[Go_intoASV_graphs] mainGroup has more than four groups; only the first four will be plotted.")
      main_group_levels <- main_group_levels[seq_len(4)]
    }
  }

  group_cloud_sets <- NULL
  group_stats_df   <- NULL
  if (!is.null(mainGroup)) {
    sample_df  <- data.frame(phyloseq::sample_data(psIN))
    group_levels <- main_group_levels
    group_cloud_sets <- lapply(group_levels, function(grp) {
      keep_samples <- rownames(sample_df)[as.character(sample_df[[mainGroup]]) == grp]
      ps_grp <- phyloseq::prune_samples(keep_samples, psIN)
      ps_grp <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_grp) > 0, ps_grp)
      .build_cloud_set(ps_grp)
    })
    names(group_cloud_sets) <- group_levels
    group_cloud_sets <- group_cloud_sets[!vapply(group_cloud_sets, is.null, logical(1))]
    group_levels <- names(group_cloud_sets)

    stats_per_group <- lapply(names(group_cloud_sets), function(grp) {
      df <- group_cloud_sets[[grp]]$stats_df
      df[[mainGroup]] <- grp
      df
    })
    if (length(stats_per_group) > 0) {
      group_stats_df <- do.call(rbind, stats_per_group)
      utils::write.csv(
        group_stats_df,
        file.path(dir_base, sprintf("sequence_cloud_stats_%s_%s.csv", mainGroup, date_tag)),
        row.names = FALSE
      )
    }
  }

  if (is.null(species_palette)) {
    palette_stats_df <- stats_df
    if (!is.null(group_stats_df) && all(c("Species", "Genus") %in% colnames(group_stats_df))) {
      palette_stats_df <- unique(rbind(
        stats_df[, c("Species", "Genus"), drop = FALSE],
        group_stats_df[, c("Species", "Genus"), drop = FALSE]
      ))
    }
    species_palette <- .build_genus_palette(palette_stats_df)
  }

  .make_legend_labels <- function(sequence_cloud_lists_i) {
    if (is.null(sequence_cloud_lists_i) || length(sequence_cloud_lists_i) == 0) return(character())
    legend_df <- do.call(rbind, lapply(sequence_cloud_lists_i, function(cloud_list_i) {
      do.call(rbind, lapply(names(cloud_list_i), function(cl) {
        stats_i <- .sequence_cloud_stats(cloud_list_i[[cl]])
        data.frame(
          ClusterID = cl,
          SpeciesLabel = stats_i$Species,
          stringsAsFactors = FALSE
        )
      }))
    }))
    legend_df <- unique(legend_df[!is.na(legend_df$SpeciesLabel) & nzchar(legend_df$SpeciesLabel), , drop = FALSE])
    legend_labels <- tapply(legend_df$ClusterID, legend_df$SpeciesLabel, function(x) {
      paste(sort(unique(x)), collapse = ", ")
    })
    legend_labels <- sprintf("%s (%s)", names(legend_labels), unname(legend_labels))
    names(legend_labels) <- names(tapply(legend_df$ClusterID, legend_df$SpeciesLabel, length))
    if (length(legend_labels) > 30) legend_labels <- legend_labels[seq_len(30)]
    legend_labels
  }

  main_group_legend_labels <- NULL
  if (!is.null(mainGroup) && !is.null(group_cloud_sets) && length(group_cloud_sets) > 0) {
    main_group_legend_labels <- .make_legend_labels(lapply(group_cloud_sets, function(x) x$sequence_cloud_list))
  }

  for (layout_i in c("sequence", "cloud")) {
    cluster_layout <- ifelse(layout_i == "sequence", "transformed_mds", "equal_circle")
    plot_title <- ifelse(
      layout_i == "sequence",
      "intoASV sequence-distance layout",
      "intoASV cloud layout"
    )
    plot_subtitle <- ifelse(
      layout_i == "sequence",
      "Cluster positions reflect transformed between-cluster sequence distance\nClouds show intra-cluster ASV topology using sparse MST + kNN edges",
      "Cluster positions are evenly arranged for visual comparison\nClouds show intra-cluster ASV topology using sparse MST + kNN edges"
    )
    weighting_tag <- .extract_weighting(c(project, cluster_map_file))
    similarity_tag <- .extract_similarity(c(project, cluster_map_file))
    name_tag <- .safe_file_token(name)
    main_group_tag <- if (is.null(mainGroup)) {
      NA_character_
    } else {
      paste(c(.safe_file_token(mainGroup), vapply(main_group_levels, .safe_file_token, character(1))), collapse = "_")
    }
    file_tokens <- c("sequence_cloud", layout_i, weighting_tag, similarity_tag, main_group_tag, name_tag, date_tag)
    file_tokens <- file_tokens[!is.na(file_tokens) & nzchar(file_tokens)]
    plot_file <- file.path(dir_base, paste0(paste(file_tokens, collapse = "_"), ".pdf"))
    plot_obj <- NULL

    .draw_one <- function(sequence_cloud_list_i, panel_label = NULL, zoom_side = "right", show_legend = TRUE, center_map = NULL, cluster_order = NULL, return_parts = FALSE, mirror_x = FALSE, legend_labels_override = NULL) {
    sequence_cloud_list <- sequence_cloud_list_i
    inset_clusters <- NULL
    cloud_names <- if (is.null(cluster_order)) names(sequence_cloud_list) else intersect(cluster_order, names(sequence_cloud_list))
    sequence_cloud_list <- sequence_cloud_list[cloud_names]
    cloud_stats <- do.call(rbind, lapply(sequence_cloud_list, .sequence_cloud_stats))
    cloud_stats <- cloud_stats[match(cloud_names, cloud_stats$ClusterID), , drop = FALSE]
    if (is.null(center_map)) {
      centers <- .compute_centers(sequence_cloud_list, cluster_layout)$centers
    } else {
      centers <- center_map[match(cloud_names, center_map$ClusterID), , drop = FALSE]
    }
    if (isTRUE(mirror_x)) {
      centers$cx <- -centers$cx
    }

    overview_nodes <- do.call(rbind, lapply(cloud_names, function(cl) {
      obj <- sequence_cloud_list[[cl]]
      lay <- obj$layout
      nodes <- as.data.frame(igraph::vertex_attr(obj$graph), stringsAsFactors = FALSE)
      lay <- lay[match(nodes$name, lay$name), , drop = FALSE]
      max_span <- max(diff(range(lay$x)), diff(range(lay$y)), 1e-9)
      stats_i <- .sequence_cloud_stats(obj)
      scale_i <- .rescale(stats_i$mean_pairwise_distance, to = cluster_scale_range,
                          from = range(cloud_stats$mean_pairwise_distance, na.rm = TRUE))
      center_i <- centers[centers$ClusterID == cl, ]
      local_x <- (lay$x - mean(lay$x)) / max_span * scale_i
      if (isTRUE(mirror_x)) local_x <- -local_x
      nodes$plot_x <- center_i$cx + local_x
      nodes$plot_y <- center_i$cy + (lay$y - mean(lay$y)) / max_span * scale_i
      nodes$ClusterID <- cl
      nodes$SpeciesLabel <- stats_i$Species
      nodes
    }))

    overview_edges <- do.call(rbind, lapply(cloud_names, function(cl) {
      obj <- sequence_cloud_list[[cl]]
      ed <- obj$edges
      nodes <- overview_nodes[overview_nodes$ClusterID == cl, ]
      ed$x <- nodes$plot_x[match(ed$from, nodes$name)]
      ed$y <- nodes$plot_y[match(ed$from, nodes$name)]
      ed$xend <- nodes$plot_x[match(ed$to, nodes$name)]
      ed$yend <- nodes$plot_y[match(ed$to, nodes$name)]
      ed$SpeciesLabel <- unique(nodes$SpeciesLabel)[1]
      ed
    }))
    if (is.null(legend_labels_override)) {
      legend_labels <- tapply(overview_nodes$ClusterID, overview_nodes$SpeciesLabel, function(x) {
        paste(sort(unique(x)), collapse = ", ")
      })
      legend_labels <- sprintf("%s (%s)", names(legend_labels), unname(legend_labels))
      names(legend_labels) <- names(tapply(overview_nodes$ClusterID, overview_nodes$SpeciesLabel, length))
      if (length(legend_labels) > 30) legend_labels <- legend_labels[seq_len(30)]
    } else {
      legend_labels <- legend_labels_override
    }
    legend_cols <- min(5, max(1, ceiling(length(legend_labels) / 6)))

    main_p <- ggplot2::ggplot() +
      ggplot2::geom_segment(
        data = overview_edges,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = SpeciesLabel, alpha = similarity),
        linewidth = edge_width, lineend = "round"
      ) +
      ggplot2::geom_point(
        data = overview_nodes,
        ggplot2::aes(x = plot_x, y = plot_y, fill = SpeciesLabel, size = total_abundance),
        shape = 21, stroke = node_stroke, color = "white", alpha = 0.96
      )

    if (isTRUE(show_circle)) {
      main_p <- main_p +
        ggplot2::annotate(
          "path",
          x = circle_radius * cos(seq(0, 2 * pi, length.out = 700)),
          y = circle_radius * sin(seq(0, 2 * pi, length.out = 700)),
          linewidth = 0.45, color = "grey20"
        )
    }

    fill_guide <- if (isTRUE(show_legend)) {
      ggplot2::guide_legend(
        ncol = legend_cols,
        byrow = TRUE,
        override.aes = list(size = 4.2, alpha = 1, stroke = 0.15)
      )
    } else {
      "none"
    }

    main_p <- main_p +
      ggplot2::scale_color_manual(values = species_palette, breaks = names(legend_labels), labels = legend_labels, na.value = "grey65") +
      ggplot2::scale_fill_manual(values = species_palette, breaks = names(legend_labels), labels = legend_labels, na.value = "grey65") +
      ggplot2::scale_alpha(range = edge_alpha_range, guide = "none") +
      ggplot2::scale_size_continuous(range = node_size_range, trans = "sqrt", guide = "none") +
      ggplot2::coord_equal(xlim = c(-axis_limit, axis_limit), ylim = c(-axis_limit, axis_limit), clip = "off") +
      ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(
        legend.position = if (isTRUE(show_legend)) "bottom" else "none",
        legend.justification = c(0.5, 0),
        legend.direction = "horizontal",
        legend.background = ggplot2::element_rect(fill = "white", color = NA),
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 6.5, face = "italic"),
        legend.key.size = grid::unit(3.0, "mm"),
        plot.margin = ggplot2::margin(0, 0, 0, 0),
        plot.title = ggplot2::element_text(face = "bold", size = ifelse(is.null(panel_label), 18, 11), hjust = 0.5, margin = ggplot2::margin(b = 0)),
        plot.subtitle = ggplot2::element_text(size = ifelse(is.null(panel_label), 9.5, 6.5), hjust = 0.5, margin = ggplot2::margin(b = 0))
      ) +
      ggplot2::guides(
        fill = fill_guide,
        color = "none"
      ) +
      ggplot2::labs(
        title = ifelse(is.null(panel_label), plot_title, paste0(plot_title, " | ", panel_label)),
        subtitle = plot_subtitle
      )

    if (!is.null(inset_layout) && is.null(inset_clusters)) {
      inset_clusters <- unique(unlist(inset_layout, use.names = FALSE))
    }
    if (is.null(inset_clusters)) {
      inset_clusters <- cloud_stats$ClusterID[order(-cloud_stats$cloud_dispersion, -cloud_stats$modularity)]
      inset_clusters <- head(inset_clusters, min(4, length(inset_clusters)))
    }
    inset_clusters <- intersect(inset_clusters, cloud_names)

    if (length(inset_clusters) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
      zoom_label_map <- if (!is.null(inset_layout)) {
        data.frame(
          ClusterID = unique(unlist(inset_layout, use.names = FALSE)),
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(ClusterID = inset_clusters, stringsAsFactors = FALSE)
      }
      zoom_label_map <- zoom_label_map[zoom_label_map$ClusterID %in% inset_clusters, , drop = FALSE]
      zoom_labels <- merge(
        centers,
        zoom_label_map,
        by = "ClusterID",
        all.y = TRUE,
        sort = FALSE
      )
      zoom_labels$label <- zoom_labels$ClusterID
      zoom_labels$nudge_x <- ifelse(zoom_labels$cx >= 0, 0.22, -0.22)
      zoom_labels$nudge_y <- ifelse(zoom_labels$cy >= 0, 0.12, -0.12)
      main_p <- main_p +
        ggrepel::geom_label_repel(
          data = zoom_labels,
          ggplot2::aes(x = cx, y = cy, label = label),
          nudge_x = zoom_labels$nudge_x,
          nudge_y = zoom_labels$nudge_y,
          size = 2.4,
          fontface = "bold",
          label.size = 0.15,
          label.r = grid::unit(0.04, "lines"),
          label.padding = grid::unit(0.10, "lines"),
          fill = "white",
          color = "grey10",
          segment.color = "grey35",
          segment.size = 0.25,
          box.padding = 0.20,
          point.padding = 0.20,
          min.segment.length = 0,
          seed = 123,
          show.legend = FALSE
        )
    }

    inset_plots <- lapply(inset_clusters, function(cl) {
      obj <- sequence_cloud_list[[cl]]
      stats_i <- .sequence_cloud_stats(obj)
      nodes <- as.data.frame(igraph::vertex_attr(obj$graph), stringsAsFactors = FALSE)
      nodes <- merge(nodes, obj$layout, by = "name", all.x = TRUE, sort = FALSE)
      nodes$x <- nodes$x - mean(nodes$x, na.rm = TRUE)
      nodes$y <- nodes$y - mean(nodes$y, na.rm = TRUE)
      max_r <- max(sqrt(nodes$x^2 + nodes$y^2), na.rm = TRUE)
      if (is.finite(max_r) && max_r > 0) {
        nodes$x <- nodes$x / max_r * 0.92
        nodes$y <- nodes$y / max_r * 0.92
      }
      ed <- obj$edges
      ed$x <- nodes$x[match(ed$from, nodes$name)]
      ed$y <- nodes$y[match(ed$from, nodes$name)]
      ed$xend <- nodes$x[match(ed$to, nodes$name)]
      ed$yend <- nodes$y[match(ed$to, nodes$name)]
      circle_df <- data.frame(
        x = cos(seq(0, 2 * pi, length.out = 360)),
        y = sin(seq(0, 2 * pi, length.out = 360))
      )
      p_inset <- ggplot2::ggplot() +
        ggplot2::geom_path(
          data = circle_df,
          ggplot2::aes(x = x, y = y),
          color = "grey70", linewidth = 0.28
        ) +
        ggplot2::geom_segment(
          data = ed,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend, alpha = similarity, linewidth = similarity),
          color = "grey25", lineend = "round", show.legend = FALSE
        ) +
        ggplot2::geom_point(
          data = nodes,
          ggplot2::aes(x = x, y = y, size = total_abundance, fill = cloud_module),
          shape = 21, color = "white", stroke = node_stroke, alpha = 0.96, show.legend = FALSE
        ) +
        ggplot2::scale_alpha(range = edge_alpha_range) +
        ggplot2::scale_linewidth(range = inset_edge_width_range) +
        ggplot2::scale_size_continuous(range = inset_node_size_range, trans = "sqrt") +
        ggplot2::coord_equal(xlim = c(-1.00, 1.00), ylim = c(-1.00, 1.00), clip = "off") +
        ggplot2::theme_void(base_size = 8) +
        ggplot2::theme(
          plot.background = ggplot2::element_rect(fill = "white", color = NA),
          panel.background = ggplot2::element_rect(fill = "white", color = NA),
          plot.margin = ggplot2::margin(5, 5, 5, 5),
          plot.title = ggplot2::element_text(face = "italic", size = 8.5, hjust = 0.5),
          plot.subtitle = ggplot2::element_text(face = "plain", size = 6.8, hjust = 0.5)
        ) +
        ggplot2::labs(
          title = stats_i$Species,
          subtitle = sprintf("%s (ASV n=%s)\ndisp=%.2f mod=%.2f", cl, stats_i$n_nodes, stats_i$cloud_dispersion, stats_i$modularity)
        )
      p_inset
    })
    names(inset_plots) <- inset_clusters

    zoom_panel <- NULL
    zoom_width <- 0
    if (length(inset_plots) > 0) {
      if (!is.null(inset_layout)) {
        if (!is.list(inset_layout)) {
          stop("[Go_intoASV_graphs] inset_layout must be a list of cluster ID vectors.")
        }
        inset_columns <- lapply(inset_layout, function(x) {
          x <- intersect(x, names(inset_plots))
          if (length(x) == 0) return(NULL)
          patchwork::wrap_plots(inset_plots[x], ncol = 1)
        })
        inset_columns <- inset_columns[!vapply(inset_columns, is.null, logical(1))]
        if (length(inset_columns) == 0) {
          zoom_panel <- patchwork::wrap_plots(inset_plots, ncol = 1)
        } else {
          zoom_panel <- patchwork::wrap_plots(inset_columns, nrow = 1)
        }
      } else {
        zoom_panel <- patchwork::wrap_plots(inset_plots, ncol = 1)
      }
      zoom_width <- if (!is.null(inset_layout)) max(1.50, 1.30 * length(inset_layout)) else 1.50
      if (identical(zoom_side, "left")) {
        plot_obj <- patchwork::wrap_plots(list(zoom_panel, main_p), ncol = 2, widths = c(zoom_width, 4.2))
      } else {
        plot_obj <- patchwork::wrap_plots(list(main_p, zoom_panel), ncol = 2, widths = c(4.2, zoom_width))
      }
    } else {
      plot_obj <- main_p
    }
    if (isTRUE(return_parts)) {
      return(list(
        main = main_p,
        zoom = zoom_panel,
        zoom_width = zoom_width,
        combined = plot_obj
      ))
    }
    plot_obj
    }
    if (is.null(mainGroup)) {
      plot_obj <- .draw_one(sequence_cloud_list, panel_label = NULL, zoom_side = "right", show_legend = TRUE)
    } else {
      shared_center_info <- if (length(group_levels) == 2) .compute_centers(sequence_cloud_list, cluster_layout) else NULL
      if (length(group_levels) == 2) {
        left_parts <- .draw_one(
          group_cloud_sets[[group_levels[1]]]$sequence_cloud_list,
          panel_label = paste0(mainGroup, ": ", group_levels[1]),
          zoom_side = "left",
          show_legend = TRUE,
          center_map = if (is.null(shared_center_info)) NULL else shared_center_info$centers,
          cluster_order = if (is.null(shared_center_info)) NULL else shared_center_info$cloud_names,
          return_parts = TRUE,
          mirror_x = FALSE,
          legend_labels_override = main_group_legend_labels
        )
        right_parts <- .draw_one(
          group_cloud_sets[[group_levels[2]]]$sequence_cloud_list,
          panel_label = paste0(mainGroup, ": ", group_levels[2]),
          zoom_side = "right",
          show_legend = FALSE,
          center_map = if (is.null(shared_center_info)) NULL else shared_center_info$centers,
          cluster_order = if (is.null(shared_center_info)) NULL else shared_center_info$cloud_names,
          return_parts = TRUE,
          mirror_x = TRUE,
          legend_labels_override = main_group_legend_labels
        )

        bilateral_parts <- list()
        bilateral_widths <- numeric()
        if (!is.null(left_parts$zoom)) {
          bilateral_parts <- c(bilateral_parts, list(left_parts$zoom))
          bilateral_widths <- c(bilateral_widths, left_parts$zoom_width)
        }
        bilateral_parts <- c(bilateral_parts, list(left_parts$main, right_parts$main))
        bilateral_widths <- c(bilateral_widths, 4.2, 4.2)
        if (!is.null(right_parts$zoom)) {
          bilateral_parts <- c(bilateral_parts, list(right_parts$zoom))
          bilateral_widths <- c(bilateral_widths, right_parts$zoom_width)
        }

        plot_obj <- patchwork::wrap_plots(
          bilateral_parts,
          ncol = length(bilateral_parts),
          widths = bilateral_widths,
          guides = "collect"
        ) &
          ggplot2::theme(legend.position = "bottom")
      } else {
        group_plots_i <- lapply(seq_along(group_levels), function(i) {
          grp <- group_levels[i]
          .draw_one(
            group_cloud_sets[[grp]]$sequence_cloud_list,
            panel_label = paste0(mainGroup, ": ", grp),
            zoom_side = ifelse(i %% 2 == 1, "left", "right"),
            show_legend = TRUE,
            center_map = if (is.null(shared_center_info)) NULL else shared_center_info$centers,
            cluster_order = if (is.null(shared_center_info)) NULL else shared_center_info$cloud_names,
            legend_labels_override = main_group_legend_labels
          )
        })
        group_plots_i <- group_plots_i[!vapply(group_plots_i, is.null, logical(1))]
        if (length(group_plots_i) == 0) stop("[Go_intoASV_graphs] No mainGroup panel could be built.")
        panel_n <- length(group_plots_i)
        panel_ncol <- ifelse(panel_n == 1, 1, 2)
        panel_nrow <- ceiling(panel_n / panel_ncol)
        plot_obj <- patchwork::wrap_plots(group_plots_i, ncol = panel_ncol, nrow = panel_nrow, guides = "collect") &
          ggplot2::theme(legend.position = "bottom")
      }
    }
    ggplot2::ggsave(plot_file, plot_obj, width = width, height = height)
    plot_files[[layout_i]] <- plot_file
    plot_objs[[layout_i]] <- plot_obj
    message("[Go_intoASV_graphs] Finished: ", plot_file)
  }

  invisible(list(
    target_clusters = target_clusters,
    graph_list = graph_list,
    sequence_cloud_list = sequence_cloud_list,
    sequence_cloud_stats = stats_df,
    group_stats = group_stats_df,
    plot = plot_objs,
    plot_file = plot_files,
    mainGroup_stats = group_stats_df,
    dir_base = dir_base,
    cluster_map_file = cluster_map_file
  ))
}
