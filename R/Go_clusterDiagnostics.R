#' Cluster Diagnostics for Go_intoASVs
#'
#' Computes per-cluster internal diversity (mean, median, max), silhouette width,
#' and generates diagnostic plots + cluster-quality summary table.
#'
#' @author Heekuk Park <hp2523@cumc.columbia.edu>
#' Created on 2025-11-16
#'
#' @param project Project name (used for output folder naming)
#' @param cluster_map CSV file containing ASV and ClusterID columns
#' @param out_dir Optional output directory (default = auto path)
#'
#' @return Saves diagnostics CSV + PDF plots; returns NULL (invisible)
#'
#'
#' @examples
#' \dontrun{
#' Go_clusterDiagnostics(
#'     project = "CervicalMicrobiome",
#'     cluster_map = "cluster_map_similarity_0.99.csv"
#' )
#' }
#'
#' @export

Go_clusterDiagnostics <- function(
    project = "",
    cluster_map,
    out_dir = NULL
){

  if(!is.null(dev.list())) dev.off()
  ###############################################
  # 0. Output directory structure
  ###############################################
  date_tag <- format(Sys.Date(), "%y%m%d")
  if (is.null(out_dir)) {
    out_dir <- sprintf("%s_%s/intoASV/cluster_diagnostics", project, date_tag)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ###############################################
  # 1. Load cluster_map
  ###############################################
  cm <- read.csv(cluster_map, check.names = FALSE)

  if (!all(c("ASV", "ClusterID") %in% colnames(cm))) {
    stop("cluster_map must contain columns: 'ASV' and 'ClusterID'.")
  }

  # ASV sequences must be the ASV names
  seqs_all <- cm$ASV
  names(seqs_all) <- cm$ASV

  ###############################################
  # 2. Prepare new diagnostic columns
  ###############################################
  cm$cluster_size     <- NA_integer_
  cm$internal_mean    <- NA_real_
  cm$internal_median  <- NA_real_
  cm$internal_max     <- NA_real_
  cm$silhouette       <- NA_real_

  ###############################################
  # 3. Precompute global distance matrix (silhouette)
  ###############################################
  dna_all <- DNAStringSet(cm$ASV)
  dm_global <- as.matrix(
    DECIPHER::DistanceMatrix(dna_all, includeTerminalGaps = FALSE)
  )

  cluster_numeric <- as.integer(factor(cm$ClusterID))

  ###############################################
  # 4. Per-cluster calculations
  ###############################################
  all_clusters <- unique(cm$ClusterID)

  for (cl in all_clusters) {
    idx <- which(cm$ClusterID == cl)
    asvs <- cm$ASV[idx]
    n_asv <- length(asvs)

    cm$cluster_size[idx] <- n_asv

    if (n_asv == 1) {
      cm$internal_mean[idx]   <- NA
      cm$internal_median[idx] <- NA
      cm$internal_max[idx]    <- NA
      cm$silhouette[idx]      <- NA
      next
    }

    # internal diversity
    dna_sub <- DNAStringSet(asvs)
    dm_sub <- as.matrix(
      DECIPHER::DistanceMatrix(dna_sub, includeTerminalGaps = FALSE)
    )
    upper <- dm_sub[upper.tri(dm_sub)]

    cm$internal_mean[idx]   <- mean(upper, na.rm = TRUE)
    cm$internal_median[idx] <- median(upper, na.rm = TRUE)
    cm$internal_max[idx]    <- max(upper, na.rm = TRUE)

    # silhouette (may fail → NA)
    sil_values <- tryCatch({
      sil <- cluster::silhouette(cluster_numeric, dist(dm_global))
      sil[match(asvs, rownames(sil)), "sil_width"]
    }, error = function(e) rep(NA, n_asv))

    cm$silhouette[idx] <- sil_values
  }

  ###############################################
  # 5. Save extended diagnostics table
  ###############################################
  diag_file <- sprintf("%s/cluster_diagnostics_%s.csv", out_dir, date_tag)
  write.csv(cm, diag_file, row.names = FALSE)

  message("[Go_clusterDiagnostics] Diagnostics table saved:")
  message(diag_file)

  ###############################################
  # 6. Minimal plots
  ###############################################
  message("[Go_clusterDiagnostics] Generating plots...")

  # 1) internal_mean histogram
  p1 <- ggplot(cm, aes(internal_mean)) +
    geom_histogram(bins=40, fill="steelblue", alpha=0.7) +
    theme_classic(base_size=14) +
    labs(title="Internal Mean Distribution",
         x="Internal Mean", y="Count")
  ggsave(sprintf("%s/internal_mean_distribution.pdf", out_dir), p1, width=6, height=4)

  # 2) internal_max histogram
  p2 <- ggplot(cm, aes(internal_max)) +
    geom_histogram(bins=40, fill="tomato", alpha=0.7) +
    theme_classic(base_size=14) +
    labs(title="Internal Max Distribution",
         x="Internal Max", y="Count")
  ggsave(sprintf("%s/internal_max_distribution.pdf", out_dir), p2, width=6, height=4)

  # 3) cluster_size vs internal_mean
  p3 <- ggplot(cm, aes(cluster_size, internal_mean)) +
    geom_point(alpha=0.7, color="purple") +
    theme_bw(base_size=14) +
    labs(title="Cluster Size vs Internal Mean",
         x="Cluster size (ASV count)", y="Internal mean")
  ggsave(sprintf("%s/cluster_size_vs_internal_mean.pdf", out_dir), p3, width=6, height=4)

  ###############################################
  # 7. NEW: Create cluster quality summary table
  ###############################################
  quality_df <- cm %>%
    group_by(ClusterID) %>%
    summarise(
      cluster_size  = cluster_size[1],
      internal_mean = internal_mean[1],
      internal_max  = internal_max[1],
      Genus         = Genus[1],
      Species       = Species[1]
    ) %>%
    ungroup() %>%
    mutate(
      quality_label = case_when(
        cluster_size >= 10 & internal_mean < 0.01 ~ "High-confidence",
        cluster_size >= 10 & internal_mean >= 0.01 ~ "High-size but diverse",
        cluster_size >= 4  & cluster_size < 10 & internal_mean < 0.02 ~ "Moderate",
        cluster_size >= 4  & cluster_size < 10 & internal_mean >= 0.02 ~ "Moderate but noisy",
        cluster_size %in% c(2,3) ~ "Low-confidence (small cluster)",
        cluster_size == 1 ~ "Singleton",
        TRUE ~ "Unclassified"
      )
    )

  quality_file <- sprintf("%s/cluster_quality_%s.csv", out_dir, date_tag)
  write.csv(quality_df, quality_file, row.names = FALSE)

  message("[Go_clusterDiagnostics] Cluster quality summary saved:")
  message(quality_file)


  ###############################################
  # 8. NEW: Per-cluster alignment visualization
  ###############################################
  message("[Go_clusterDiagnostics] Generating cluster alignment plots...")

  plot_dir <- sprintf("%s/Cluster_plots", out_dir)
  dir.create(plot_dir, showWarnings = FALSE)

  # make subfolders by quality label
  qual_levels <- unique(quality_df$quality_label)
  for (q in qual_levels) {
    if (q == "Singleton") {
      message(sprintf("[Go_clusterDiagnostics] Cluster %s labeled as Singleton — skipping output.", cl))
      next
    }else{
      dir.create(sprintf("%s/%s", plot_dir, q), showWarnings = FALSE)
    }
  }

  max_full_asv <- 40     # full alignment limit

  for (i in seq_len(nrow(quality_df))) {

    cl <- quality_df$ClusterID[i]
    qlab <- quality_df$quality_label[i]

    idx <- which(cm$ClusterID == cl)
    asvs <- cm$ASV[idx]
    n_asv <- length(asvs)

    # ---- EARLY SKIP: singleton ----


    out_subdir <- sprintf("%s/%s", plot_dir, qlab)


    if (n_asv < 2) {
      message(sprintf("[Go_clusterDiagnostics] Cluster %s has %d ASVs — skipping alignment.", cl, n_asv))
      next
    }

    # full alignment
    dna_sub <- DNAStringSet(asvs)
    aln <- DECIPHER::AlignSeqs(dna_sub, iterations = 1, refinements = 1, verbose = FALSE)

    ### Always generate alignment HTML ###
    html_file <- sprintf("%s/%s_alignment.html", out_subdir, cl)
    DECIPHER::BrowseSeqs(aln, htmlFile = html_file, open= FALSE)

    ### Distance matrix ###
    dm <- DECIPHER::DistanceMatrix(aln, includeTerminalGaps = FALSE)
    dm_mat <- as.matrix(dm)
    dm_mat[is.na(dm_mat)] <- 0

    ### CASE 1: ASV = 2 (tree not informative) ###
    if (n_asv == 2) {
      message(sprintf("[Go_clusterDiagnostics] Cluster %s: only 2 ASVs – skipping tree.", cl))

      pdf(sprintf("%s/%s_alignment_tree.pdf", out_subdir, cl), width=7, height=4)
      plot.new()
      text(0.5, 0.6,
           sprintf("Cluster %s has only 2 ASVs.\nTree is not informative.", cl),
           cex = 1)
      text(0.5, 0.4,
           sprintf("Pairwise distance = %.4f", dm_mat[1,2]),
           cex=0.9, col="blue")
      dev.off()

      next
    }

    ### CASE 2: all distances = 0 → identical cluster ###
    if (all(dm_mat == 0)) {
      message(sprintf("[Go_clusterDiagnostics] Cluster %s: all distances zero — identical ASVs.", cl))

      pdf(sprintf("%s/%s_alignment_tree.pdf", out_subdir, cl), width=7, height=4)
      plot.new()
      text(0.5, 0.6,
           sprintf("Cluster %s: All ASVs are identical.\nTree cannot be constructed.", cl),
           cex = 1)
      dev.off()

      next
    }

    ### Add tiny noise to avoid metric issues ###
    epsilon <- 1e-8
    dm_mat <- dm_mat + epsilon
    diag(dm_mat) <- 0

    ### enforce dimnames ###
    if (is.null(rownames(dm_mat))) {
      rownames(dm_mat) <- colnames(dm_mat) <- paste0("ASV", seq_len(nrow(dm_mat)))
    }

    ### Tree construction ###
    d <- as.dist(dm_mat)
    hc <- try(hclust(d), silent = TRUE)

    if (inherits(hc, "try-error")) {
      message(sprintf("[Go_clusterDiagnostics] hclust() failed for %s — skipping tree.", cl))
      next
    }

    ### CASE 3: flat dendrogram ###
    if (length(unique(hc$height)) == 1) {
      message(sprintf("[Go_clusterDiagnostics] Cluster %s: flat dendrogram.", cl))

      pdf(sprintf("%s/%s_alignment_tree.pdf", out_subdir, cl), width=7, height=4)
      plot.new()
      text(0.5, 0.6,
           sprintf("Cluster %s dendrogram is flat.\nNo informative height variation.", cl),
           cex = 1)
      text(0.5, 0.4,
           sprintf("Height = %.4f", hc$height[1]),
           col = "blue", cex = 0.9)
      dev.off()

      next
    }

    ### CASE 4: Normal tree ###
    pdf(sprintf("%s/%s_alignment_tree.pdf", out_subdir, cl), width=10, height=5)
    par(mar=c(10,4,4,2))
    plot(
      hc,
      hang=-1,
      labels=rownames(dm_mat),
      cex=0.6,
      las=2,
      main=sprintf("%s (n=%d ASVs) – distance-based tree", cl, n_asv)
    )
    dev.off()

  }

  message("[Go_clusterDiagnostics] Alignment plots saved in: ", plot_dir)

  ###############################################
  # 9. ASV-level heatmap (pheatmap, final stable)
  ###############################################
  message("[Go_clusterDiagnostics] Generating ASV-level heatmap (pheatmap, final)...")


  # 1) Distance matrix (ASV × ASV)
  dist_mat <- dm_global
  if (is.null(rownames(dist_mat))) {
    rownames(dist_mat) <- cm$ASV
    colnames(dist_mat) <- cm$ASV
  }

  # 2) Annotation (ClusterID only)
  meta_order <- cm[match(rownames(dist_mat), cm$ASV), ]

  anno_df <- data.frame(
    ClusterID = meta_order$ClusterID,
    row.names = rownames(dist_mat)
  )

  # 3) Annotation colors (ClusterID only)
  cluster_ids <- unique(meta_order$ClusterID)
  n_cluster <- length(cluster_ids)

  set.seed(2025)

  cluster_colors <- Polychrome::createPalette(
    n_cluster,
    seedcolors = c("#FF0000", "#00FF00", "#0000FF", "#FFFF00")
  )

  cluster_colors <- setNames(cluster_colors, cluster_ids)

  anno_colors <- list(ClusterID = cluster_colors)

  # 4) FIXED: hclust ordering (pheatmap & annotation 둘 다 동일하게)
  row_hc <- hclust(as.dist(dist_mat), method = "average")
  col_hc <- hclust(as.dist(dist_mat), method = "average")

  row_order <- row_hc$order
  col_order <- col_hc$order

  # heatmap matrix reorder
  dist_ord <- dist_mat[row_order, col_order]

  # annotation reorder – 이게 핵심 수정!
  anno_ord <- anno_df[row_order, , drop = FALSE]

  # 5) Save stable PDF
  heatmap_file_pdf <- sprintf("%s/ASV_distance_heatmap_%s.pdf", out_dir, date_tag)

  pdf(heatmap_file_pdf, width = 14, height = 14)
  pheatmap::pheatmap(
    dist_ord,
    clustering_method = "average",
    cluster_rows = FALSE,   # 이미 정렬됨
    cluster_cols = FALSE,   # 이미 정렬됨
    color = colorRampPalette(c("white", "yellow", "orange", "red"))(256),
    annotation_row = anno_ord,
    annotation_col = anno_ord,
    annotation_colors = anno_colors,
    legend = TRUE,                # scale bar YES
    annotation_legend = FALSE,    # no annotation legend
    show_rownames = FALSE,
    show_colnames = FALSE,
    fontsize = 7,
    main = "ASV-level Pairwise Distance Heatmap"
  )
  dev.off()

  message("[Go_clusterDiagnostics] ASV heatmap saved:")
  message(heatmap_file_pdf)


  # row/column reordered ASV names
  asv_order_row <- rownames(dist_ord)
  asv_order_col <- colnames(dist_ord)

  # create output file paths
  order_row_file <- sprintf("%s/ASV_order_row_%s.txt", out_dir, date_tag)
  order_col_file <- sprintf("%s/ASV_order_col_%s.txt", out_dir, date_tag)

  # save as simple text files (one ASV per line)
  #writeLines(asv_order_row, order_row_file)
  #writeLines(asv_order_col, order_col_file)

  message("[Go_clusterDiagnostics] Saved ASV order for Python:")
  message(order_row_file)
  message(order_col_file)


  ###############################################
  # 10. Export distance matrix & annotation for Python HTML heatmap
  ###############################################
  message("[Go_clusterDiagnostics] Exporting ASV distance and annotation for Python...")

  # dm_global: ASV × ASV distance matrix (이미 위에서 계산됨)
  # cm: data.frame with at least ASV, ClusterID (그리고 있으면 quality_label)

  # 1) Distance matrix (ASV x ASV)
  dist_mat <- dm_global
  if (is.null(rownames(dist_mat))) {
    rownames(dist_mat) <- cm$ASV
    colnames(dist_mat) <- cm$ASV
  }

  dist_file <- sprintf("%s/asv_distance_matrix_%s.csv", out_dir, date_tag)
  write.csv(
    as.data.frame(dist_mat),
    dist_file,
    row.names = TRUE,
    quote = FALSE
  )

  # 2) Annotation (ASV, ClusterID, optional quality_label)
  # anno_df <- cm[, c("ASV", "ClusterID")]
  # if ("quality_label" %in% colnames(cm)) {
  #  anno_df$quality_label <- cm$quality_label
  # } else {
  #   anno_df$quality_label <- NA_character_
  # }

  #anno_file <- sprintf("%s/asv_annotation_%s.csv", out_dir, date_tag)
  # write.csv(anno_df, anno_file, row.names = FALSE, quote = FALSE)


  # Save row/column order used by pheatmap
  row_order <- hclust(as.dist(dist_mat), method = "average")$order
  col_order <- row_order  # symmetric matrix

  ordered_asv <- rownames(dist_mat)[row_order]

  # Save to CSV for Python to read
  write.csv(
    data.frame(ASV = ordered_asv),
    file = sprintf("%s/asv_order_%s.csv", out_dir, date_tag),
    row.names = FALSE
  )

  message("[Go_clusterDiagnostics] Python input files saved:")
  message(dist_file)
  # message(anno_file)

  return(invisible(NULL))
}
