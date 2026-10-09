# Title: evaluate k_param- dim red mitopaint vis (mean per well) v1
# Step: 6.5
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 07-10-2026

# load packages ####
library(data.table)
library(Seurat)
library(tidyverse)
library(cluster)
library(mclust)
library(aricode)
library(ggplot2)
library(ggpubr)
# set variables ####
file_name <- "mPaintFDA_N1_N2_N3_N4_N5_N6_N7_N8"
redu_state <- "redu"
integrate_state <- "integrated"
dims_use <- 1:50
k_values <- c(8, 12, 16, 20, 30, 40)
res <- 1
perplexity <- 20
max_iter <- 500
avg_profile <- TRUE
selected_k <- 8
# load data ####
# load data as df
if (redu_state == "redu") {
  df <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_data_", integrate_state ,"_redu.csv", sep = ""), 
      header = TRUE)
  )
} else {
  df <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_data_", integrate_state ,".csv", sep = ""), 
      header = TRUE)
  )
}
# keep rownames as WELL_BATCH
rownames(df) <- df$V1
df$V1 <- NULL
# load metadata as meta
meta <- as.data.frame(
  fread(
    paste(
      "data/processed/", file_name, "_meta_", integrate_state, ".csv", sep = ""), 
    header = TRUE)
)
# keep rownames as WELL_BATCH
rownames(meta) <- meta$V1
meta$V1 <- NULL
# remove any columns with NA/ non finite values
df <- df[, colSums(!is.finite(as.matrix(df))) == 0, drop = FALSE]
# remove any rows with NA/ non finite values
df <- df[apply(df, 1, function(x) all(is.finite(x))), , drop = FALSE]
# if avg_profile == TRUE ####
if (avg_profile) {
  # combine metadata and feature matrix
  combined <- cbind(meta, df)
  # average all non-DMSO profiles within Batch × Condition
  non_dmso <- combined |>
    dplyr::filter(Compound != "DMSO") |>
    dplyr::group_by(Batch, Condition) |>
    dplyr::summarise(
      Compound = dplyr::first(Compound),
      Concentration = dplyr::first(Concentration),
      Row = NA_integer_,
      Column = NA_integer_,
      Well = NA_character_,
      ID = NA_character_,
      Order = NA_real_,
      dplyr::across(
        all_of(colnames(df)),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      row_id = paste(
        Condition,
        Batch,
        sep = "_"
      )
    )
  # keep all DMSO wells individually
  dmso <- combined |>
    dplyr::filter(Compound == "DMSO") |>
    dplyr::mutate(
      row_id = paste(
        Condition,
        Well,
        Batch,
        sep = "_"
      )
    )
  # combine averaged non-DMSO with raw DMSO
  combined_avg <- dplyr::bind_rows(
    dmso,
    non_dmso
  )
  # assign unique row names
  rownames(combined_avg) <- combined_avg$row_id
  # split back into df and meta
  df_avg <- combined_avg |>
    dplyr::select(
      all_of(colnames(df))
    )
  meta_avg <- combined_avg |>
    dplyr::select(
      Row,
      Column,
      Compound,
      Concentration,
      Well,
      Batch,
      ID,
      Order,
      Condition
    )
  rownames(df_avg) <- rownames(combined_avg)
  rownames(meta_avg) <- rownames(combined_avg)
  # overwrite originals
  df <- df_avg
  meta <- meta_avg
  rm(
    combined_avg,
    non_dmso,
    dmso,
    combined
  )
}
# put data in Seurat object df.seurat ####
# assign data and meta to seurat object
df.seurat <- CreateSeuratObject(
  # transpose zscore matrix
  counts = t(as.matrix(df)),
  meta.data = meta,
  # assay is saved as MP for mitopaint
  assay = "MP"
)
# make mitopaint data default assay
DefaultAssay(df.seurat) <- "MP"
# manually copy counts to data layer to circumvent v5 seurat structure syntax
df.seurat <- SetAssayData(
  df.seurat,
  assay = "MP",
  layer = "data",
  new.data = GetAssayData(
    df.seurat,
    assay = "MP",
    layer = "counts"
  )
)
# scale data
df.seurat <- ScaleData(
  df.seurat,
  features = rownames(df.seurat)
)
# run PCA ####
df.seurat <- RunPCA(
  df.seurat,
  features = rownames(df.seurat),
  seed.use = 42,
  verbose = FALSE
)
# run TSNE ####
df.seurat <- RunTSNE(
  df.seurat,
  reduction = NULL,
  dims = NULL,
  features = rownames(df.seurat),
  perplexity = perplexity,
  max_iter = max_iter,
  reduction.name = "tsne",
  seed.use = 42,
  verbose = FALSE
)
# run UMAP ####
df.seurat <- RunUMAP(
  df.seurat,
  reduction = NULL,
  dims = NULL,
  features = rownames(df.seurat),
  reduction.name = "umap",
  seed.use = 42,
  verbose = FALSE
)
# helper function: cosine similarity matrix ####
calc_cosine_similarity <- function(x) {
  # save profile IDs
  profile_ids <- rownames(x)
  # convert profiles to numeric matrix
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  # calculate Euclidean norm of each profile
  norms <- sqrt(rowSums(x^2))
  # identify valid profiles
  valid <- is.finite(norms) & norms > 0
  # remove profiles with zero / invalid norm if present
  if (!all(valid)) {
    warning(
      sum(!valid),
      " profiles have zero or non-finite norm and will be removed."
    )
    x <- x[valid, , drop = FALSE]
    profile_ids <- profile_ids[valid]
    norms <- norms[valid]
  }
  # normalise each profile to unit length
  x_norm <- x / norms
  # calculate pairwise cosine similarity
  sim <- x_norm %*% t(x_norm)
  # prevent tiny floating-point errors outside [-1, 1]
  # use direct indexing so matrix dimensions are retained
  sim[sim > 1] <- 1
  sim[sim < -1] <- -1
  # explicitly retain profile IDs
  rownames(sim) <- profile_ids
  colnames(sim) <- profile_ids
  return(sim)
}
# calculate cosine similarity from df ####
# similarity is based on the ORIGINAL profiling matrix,
# and not dim red coordinates (only grouping for averaging is informed by dim red)
cosine_sim <- calc_cosine_similarity(df)
# cosine distance
cosine_dist <- as.dist(1 - cosine_sim)
# calculate Pearson similarity from df ####
pearson_sim <- cor(t(as.matrix(df)),
  method = "pearson",
  use = "pairwise.complete.obs"
)
# explicitly retain profile IDs
rownames(pearson_sim) <- rownames(df)
colnames(pearson_sim) <- rownames(df)
# Pearson correlation distance
pearson_dist <- as.dist(
  1 - pearson_sim
)

# helper function: within-cluster profile similarity (cosine and pearson) ####
calc_cluster_similarity <- function(
    similarity_matrix,
    clusters,
    similarity_name
) {
  
  # align cluster labels to similarity matrix
  clusters <- clusters[
    rownames(similarity_matrix)
  ]
  
  cluster_levels <- unique(
    as.character(clusters)
  )
  
  cluster_results <- purrr::map_dfr(
    cluster_levels,
    function(cl) {
      
      ids <- names(clusters)[
        as.character(clusters) == cl
      ]
      
      n <- length(ids)
      
      # cannot calculate pairwise similarity for singleton
      if (n < 2) {
        
        return(
          tibble::tibble(
            cluster = cl,
            n_profiles = n,
            mean_similarity = NA_real_,
            median_similarity = NA_real_,
            sd_similarity = NA_real_
          )
        )
      }
      
      # similarity matrix for profiles in this cluster
      sim_sub <- similarity_matrix[
        ids,
        ids,
        drop = FALSE
      ]
      
      # use each unique pair once
      vals <- sim_sub[
        upper.tri(sim_sub)
      ]
      
      tibble::tibble(
        cluster = cl,
        n_profiles = n,
        mean_similarity = mean(
          vals,
          na.rm = TRUE
        ),
        median_similarity = median(
          vals,
          na.rm = TRUE
        ),
        sd_similarity = sd(
          vals,
          na.rm = TRUE
        )
      )
    }
  )
  
  # give similarity columns metric-specific names
  cluster_results <- cluster_results |>
    dplyr::rename(
      !!paste0("mean_", similarity_name) := mean_similarity,
      !!paste0("median_", similarity_name) := median_similarity,
      !!paste0("sd_", similarity_name) := sd_similarity
    )
  
  return(cluster_results)
}
# helper function: silhouette ####
calc_cluster_silhouette <- function(
    cosine_dist,
    clusters
) {
  clusters <- clusters[attr(cosine_dist, "Labels")]
  clusters <- factor(clusters)
  # silhouette requires at least two clusters
  if (nlevels(clusters) < 2) {
    return(NA_real_)
  }
  # silhouette is not meaningful if every point is its own cluster
  if (nlevels(clusters) >= length(clusters)) {
    return(NA_real_)
  }
  sil <- cluster::silhouette(
    as.integer(clusters),
    cosine_dist
  )
  mean(
    sil[, "sil_width"],
    na.rm = TRUE
  )
}
# helper function: evaluate one cluster assignment ####
evaluate_clusters <- function(
    cosine_sim,
    cosine_dist,
    pearson_sim,
    clusters,
    method,
    k
) {
  # ensure cluster names are present
  if (is.null(names(clusters))) {
    stop("Cluster labels must have profile names.")
  }
  # cosine similarity
  cluster_cos <- calc_cluster_similarity(
    similarity_matrix = cosine_sim,
    clusters = clusters,
    similarity_name = "cosine"
  )
  # Pearson similarity
  cluster_pearson <- calc_cluster_similarity(
    similarity_matrix = pearson_sim,
    clusters = clusters,
    similarity_name = "pearson"
  )
  # combine per-cluster results
  cluster_similarity <- cluster_cos |>
    dplyr::left_join(
      cluster_pearson,
      by = c(
        "cluster",
        "n_profiles"
      )
    ) |>
    dplyr::mutate(
      method = method,
      k = k,
      .before = 1
    )
  # silhouette
  # keep cosine distance as primary silhouette metric
  mean_silhouette <- calc_cluster_silhouette(
    cosine_dist = cosine_dist,
    clusters = clusters
  )
  # cluster sizes
  cluster_sizes <- table(clusters)
  # valid clusters for pairwise similarity
  valid_sim <- cluster_similarity |>
    dplyr::filter(
      n_profiles >= 2
    )
  # unweighted cosine
  mean_cosine_unweighted <- if (nrow(valid_sim) > 0) {
    mean(
      valid_sim$mean_cosine,
      na.rm = TRUE
    )
  } else {
    NA_real_
  }
  # unweighted Pearson
  mean_pearson_unweighted <- if (nrow(valid_sim) > 0) {
    mean(
      valid_sim$mean_pearson,
      na.rm = TRUE
    )
  } else {
    NA_real_
  }
  # pair weights
  pair_weights <- valid_sim$n_profiles *
    (valid_sim$n_profiles - 1) / 2
  # weighted cosine
  mean_cosine_weighted <- if (nrow(valid_sim) > 0) {
    weighted.mean(
      valid_sim$mean_cosine,
      w = pair_weights,
      na.rm = TRUE
    )
  } else {
    NA_real_
  }
  # weighted Pearson
  mean_pearson_weighted <- if (nrow(valid_sim) > 0) {
    weighted.mean(
      valid_sim$mean_pearson,
      w = pair_weights,
      na.rm = TRUE
    )
  } else {
    NA_real_
  }
  # summary
  summary <- tibble::tibble(
    method = method,
    k = k,
    n_clusters = length(cluster_sizes),
    min_cluster_size = min(
      as.numeric(cluster_sizes)
    ),
    median_cluster_size = median(
      as.numeric(cluster_sizes)
    ),
    mean_cluster_size = mean(
      as.numeric(cluster_sizes)
    ),
    max_cluster_size = max(
      as.numeric(cluster_sizes)
    ),
    n_singletons = sum(
      cluster_sizes == 1
    ),
    mean_within_cluster_cosine_unweighted =
      mean_cosine_unweighted,
    mean_within_cluster_cosine_weighted =
      mean_cosine_weighted,
    mean_within_cluster_pearson_unweighted =
      mean_pearson_unweighted,
    mean_within_cluster_pearson_weighted =
      mean_pearson_weighted,
    mean_silhouette =
      mean_silhouette
  )
  return(
    list(
      summary = summary,
      cluster_similarity = cluster_similarity
    )
  )
}
# containers for k tuning ####
cluster_assignments <- list(
  PCA = list(),
  tSNE = list(),
  UMAP = list()
)
k_summary_list <- list()
cluster_similarity_list <- list()
# tune k ####
for (k in k_values) {
  message(
    "Testing k = ",
    k
  )
  # PCA kNN / SNN ####
  pca_nn_name <- paste0("pca_nn_k",k)
  pca_snn_name <- paste0("pca_snn_k",k)
  pca_cluster_name <- paste0("PCA_NN_k",k)
  # find PCA NN
  df.seurat <- FindNeighbors(
    df.seurat,
    reduction = "pca",
    dims = dims_use,
    k.param = k,
    graph.name = c(
      pca_nn_name,
      pca_snn_name
    ),
    verbose = FALSE
  )
  # find PCA clusters
  df.seurat <- FindClusters(
    df.seurat,
    graph.name = pca_snn_name,
    resolution = res,
    cluster.name = pca_cluster_name,
    random.seed = 42,
    verbose = FALSE
  )
  # extract PCA cluster labels
  pca_clusters <- df.seurat@meta.data[[pca_cluster_name]]
  names(pca_clusters) <- rownames(df.seurat@meta.data)
  cluster_assignments$PCA[[as.character(k)]] <- pca_clusters
  # evaluate PCA clusters against original profiles
  pca_eval <- evaluate_clusters(
    cosine_sim = cosine_sim,
    cosine_dist = cosine_dist,
    pearson_sim = pearson_sim,
    clusters = pca_clusters,
    method = "PCA",
    k = k
  )
  k_summary_list[[paste0("PCA_", k)]] <- pca_eval$summary
  cluster_similarity_list[[paste0("PCA_", k)]] <-
    pca_eval$cluster_similarity
  # tSNE kNN / SNN ####
  tsne_nn_name <- paste0("tsne_nn_k",k)
  tsne_snn_name <- paste0("tsne_snn_k",k)
  tsne_cluster_name <- paste0("tSNE_NN_k",k)
  # find tSNE NN
  df.seurat <- FindNeighbors(
    df.seurat,
    reduction = "tsne",
    dims = 1:2,
    k.param = k,
    graph.name = c(
      tsne_nn_name,
      tsne_snn_name
    ),
    verbose = FALSE
  )
  # find tSNE clusters
  df.seurat <- FindClusters(
    df.seurat,
    graph.name = tsne_snn_name,
    resolution = res,
    cluster.name = tsne_cluster_name,
    random.seed = 42,
    verbose = FALSE
  )
  # extract tSNE cluster labels
  tsne_clusters <- df.seurat@meta.data[[tsne_cluster_name]]
  names(tsne_clusters) <- rownames(
    df.seurat@meta.data
  )
  # save assignments
  cluster_assignments$tSNE[[as.character(k)]] <-
    tsne_clusters
  # evaluate tSNE clusters against original profiles
  tsne_eval <- evaluate_clusters(
    cosine_sim = cosine_sim,
    cosine_dist = cosine_dist,
    pearson_sim = pearson_sim,
    clusters = tsne_clusters,
    method = "tSNE",
    k = k
  )
  k_summary_list[[paste0("tSNE_", k)]] <-
    tsne_eval$summary
  cluster_similarity_list[[paste0("tSNE_", k)]] <-
    tsne_eval$cluster_similarity
  # UMAP kNN / SNN ####
  umap_nn_name <- paste0("umap_nn_k",k)
  umap_snn_name <- paste0("umap_snn_k",k)
  umap_cluster_name <- paste0("UMAP_NN_k",k)
  # find UMAP NN
  df.seurat <- FindNeighbors(
    df.seurat,
    reduction = "umap",
    dims = 1:2,
    k.param = k,
    graph.name = c(
      umap_nn_name,
      umap_snn_name
    ),
    verbose = FALSE
  )
  # find UMAP clusters
  df.seurat <- FindClusters(
    df.seurat,
    graph.name = umap_snn_name,
    resolution = res,
    cluster.name = umap_cluster_name,
    random.seed = 42,
    verbose = FALSE
  )
  # extract UMAP cluster labels
  umap_clusters <- df.seurat@meta.data[[umap_cluster_name]]
  names(umap_clusters) <- rownames(df.seurat@meta.data)
  cluster_assignments$UMAP[[as.character(k)]] <- umap_clusters
  # evaluate UMAP clusters against original profiles
  umap_eval <- evaluate_clusters(
    cosine_sim = cosine_sim,
    cosine_dist = cosine_dist,
    pearson_sim = pearson_sim,
    clusters = umap_clusters,
    method = "UMAP",
    k = k
  )
  k_summary_list[[paste0("UMAP_", k)]] <- umap_eval$summary
  cluster_similarity_list[[paste0("UMAP_", k)]] <-
    umap_eval$cluster_similarity
}
# combine k summary ####
k_summary <- dplyr::bind_rows(
  k_summary_list
) |>
  dplyr::arrange(
    method,
    k
  )
# combine per-cluster cosine results ####
cluster_similarity_results <- dplyr::bind_rows(
  cluster_similarity_list
) |>
  dplyr::arrange(
    method,
    k,
    cluster
  )
# function to compare clustering across k ####
calc_cross_k_agreement <- function(
    assignment_list,
    method
) {
  k_names <- names(
    assignment_list
  )
  comparisons <- combn(
    k_names,
    2,
    simplify = FALSE
  )
  purrr::map_dfr(
    comparisons,
    function(pair) {
      k1 <- pair[1]
      k2 <- pair[2]
      cl1 <- assignment_list[[k1]]
      cl2 <- assignment_list[[k2]]
      # make absolutely sure same profiles are compared
      common_ids <- intersect(
        names(cl1),
        names(cl2)
      )
      cl1 <- factor(
        cl1[common_ids]
      )
      cl2 <- factor(
        cl2[common_ids]
      )
      tibble::tibble(
        method = method,
        k1 = as.integer(k1),
        k2 = as.integer(k2),
        ARI = mclust::adjustedRandIndex(
          cl1,
          cl2
        ),
        NMI = aricode::NMI(
          cl1,
          cl2
        )
      )
    }
  )
}
# calculate cross-k ARI / NMI ####
cross_k_results <- dplyr::bind_rows(
  calc_cross_k_agreement(
    cluster_assignments$PCA,
    method = "PCA"
  ),
  calc_cross_k_agreement(
    cluster_assignments$tSNE,
    method = "tSNE"
  ),
  calc_cross_k_agreement(
    cluster_assignments$UMAP,
    method = "UMAP"
  )
) |>
  dplyr::arrange(
    method,
    k1,
    k2
  )
# calculate agreement with adjacent k values only ####
# useful for identifying a stable plateau
adjacent_k_results <- cross_k_results |>
  dplyr::filter(
    purrr::map2_lgl(
      k1,
      k2,
      function(a, b) {
        i <- match(
          a,
          k_values
        )
        j <- match(
          b,
          k_values
        )
        abs(i - j) == 1
      }
    )
  )
# save cluster assignments in long format ####
cluster_assignment_df <- dplyr::bind_rows(
  purrr::imap_dfr(
    cluster_assignments$PCA,
    function(clusters, k) {
      tibble::tibble(
        profile_id = names(clusters),
        method = "PCA",
        k = as.integer(k),
        cluster = as.character(clusters)
      )
    }
  ),
  purrr::imap_dfr(
    cluster_assignments$tSNE,
    function(clusters, k) {
      tibble::tibble(
        profile_id = names(clusters),
        method = "tSNE",
        k = as.integer(k),
        cluster = as.character(clusters)
      )
    }
  ),
  purrr::imap_dfr(
    cluster_assignments$UMAP,
    function(clusters, k) {
      tibble::tibble(
        profile_id = names(clusters),
        method = "UMAP",
        k = as.integer(k),
        cluster = as.character(clusters)
      )
    }
  )
)
# inspect results ####
#View(k_summary)
#View(cluster_similarity_results)
#View(cross_k_results)
#View(adjacent_k_results)
# function to plot k tuning results ####
plots <- list()
plot_k_tuning <- function(
    k_summary,
    method_name,
    selected_k = NULL
) {
  # subset requested dimensionality reduction
  plot_df <- k_summary |>
    dplyr::filter(
      method == method_name
    )
  # optional selected-k line
  selected_k_layer <- if (!is.null(selected_k)) {
    ggplot2::geom_vline(
      xintercept = selected_k,
      linetype = "dashed",
      colour = "grey50"
    )
  } else {
    NULL
  }
  # cosine similarity
  p_cosine <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = k,
      y = mean_within_cluster_cosine_weighted
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      size = 2.5
    ) +
    selected_k_layer +
    ggplot2::scale_x_continuous(
      breaks = k_values
    ) +
    ggplot2::labs(
      title = paste0(
        method_name,
        "\nCosine similarity"
      ),
      x = "k",
      y = "Mean within-cluster\ncosine similarity"
    ) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 10,
        face = "bold"
      ),
      axis.text = ggplot2::element_text(
        size = 8
      ),
      axis.title = ggplot2::element_text(
        size = 9
      )
    )
  # pearson similarity
  p_pearson <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = k,
      y = mean_within_cluster_pearson_weighted
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      size = 2.5
    ) +
    selected_k_layer +
    ggplot2::scale_x_continuous(
      breaks = k_values
    ) +
    ggplot2::labs(
      title = paste0(
        method_name,
        "\nPearson similarity"
      ),
      x = "k",
      y = "Mean within-cluster\npearson similarity"
    ) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 10,
        face = "bold"
      ),
      axis.text = ggplot2::element_text(
        size = 8
      ),
      axis.title = ggplot2::element_text(
        size = 9
      )
    )
  # silhouette
  p_silhouette <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = k,
      y = mean_silhouette
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      size = 2.5
    ) +
    selected_k_layer +
    ggplot2::scale_x_continuous(
      breaks = k_values
    ) +
    ggplot2::labs(
      title = paste0(
        method_name,
        "\nMean silhouette (cosine distance)"
      ),
      x = "k",
      y = "Mean silhouette"
    ) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 10,
        face = "bold"
      ),
      axis.text = ggplot2::element_text(
        size = 8
      ),
      axis.title = ggplot2::element_text(
        size = 9
      )
    )
  # number of clusters
  p_clusters <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = k,
      y = n_clusters
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      size = 2.5
    ) +
    selected_k_layer +
    ggplot2::scale_x_continuous(
      breaks = k_values
    ) +
    ggplot2::labs(
      title = paste0(
        method_name,
        "\nNumber of clusters"
      ),
      x = "k",
      y = "Number of clusters"
    ) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 10,
        face = "bold"
      ),
      axis.text = ggplot2::element_text(
        size = 8
      ),
      axis.title = ggplot2::element_text(
        size = 9
      )
    )
  return(
    list(
      cosine = p_cosine,
      pearson = p_pearson,
      silhouette = p_silhouette,
      clusters = p_clusters
    )
  )
}
# generate k tuning plots for PCA, tSNE and UMAP ####
plots$PCA <- plot_k_tuning(
  k_summary = k_summary,
  method_name = "PCA",
  selected_k = selected_k
)
plots$tSNE <- plot_k_tuning(
  k_summary = k_summary,
  method_name = "tSNE",
  selected_k = selected_k
)
plots$UMAP <- plot_k_tuning(
  k_summary = k_summary,
  method_name = "UMAP",
  selected_k = selected_k
)
# function to plot adjacent-k stability ####
plot_adjacent_k <- function(
    results_df,
    method_name,
    metric
) {
  
  # subset dimensionality reduction
  plot_df <- results_df |>
    dplyr::filter(method == method_name)
  
  # choose y-axis label
  y_lab <- dplyr::case_when(
    metric == "ARI" ~ "Adjusted Rand Index (ARI)",
    metric == "NMI" ~ "Normalized Mutual Information (NMI)",
    TRUE ~ metric
  )
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = k2,
      y = .data[[metric]]
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      size = 2.5
    ) +
    ggplot2::scale_x_continuous(
      breaks = k_values
    ) +
    ggplot2::labs(
      title = paste0(
        method_name,
        "\nAdjacent cross-k ",
        metric
      ),
      x = "k2 in cross comparison",
      y = y_lab
    ) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 10,
        face = "bold"
      ),
      axis.text = ggplot2::element_text(size = 8),
      axis.title = ggplot2::element_text(size = 9)
    )
}
# plots ####
# PCA
plots$PCA$ARI <- plot_adjacent_k(
  adjacent_k_results,
  method_name = "PCA",
  metric = "ARI"
)

plots$PCA$NMI <- plot_adjacent_k(
  adjacent_k_results,
  method_name = "PCA",
  metric = "NMI"
)

# tSNE
plots$tSNE$ARI <- plot_adjacent_k(
  adjacent_k_results,
  method_name = "tSNE",
  metric = "ARI"
)

plots$tSNE$NMI <- plot_adjacent_k(
  adjacent_k_results,
  method_name = "tSNE",
  metric = "NMI"
)

# UMAP
plots$UMAP$ARI <- plot_adjacent_k(
  adjacent_k_results,
  method_name = "UMAP",
  metric = "ARI"
)

plots$UMAP$NMI <- plot_adjacent_k(
  adjacent_k_results,
  method_name = "UMAP",
  metric = "NMI"
)
plots
# save results ####
write.csv(
  k_summary,
  paste0(
    "outputs/data/",
    file_name,
    "_k_tuning_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  cross_k_results,
  paste0(
    "outputs/data/",
    file_name,
    "_k_tuning_cross_k_ARI_NMI.csv"
  ),
  row.names = FALSE
)

# save plots ####
purrr::iwalk(
  plots,
  function(method_plots, method_name) {
    purrr::iwalk(
      method_plots,
      function(plot, plot_name) {
        ggplot2::ggsave(
          filename = paste0(
            "outputs/figures/",
            file_name,
            "_k_tuning_",
            method_name,
            "_",
            plot_name,
            ".pdf"
          ),
          plot = plot,
          width = 3,
          height = 3,
          units = "in"
        )
      }
    )
  }
)
