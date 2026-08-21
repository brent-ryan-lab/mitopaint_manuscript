# Title: batch integration evaluation stats mitopaint data (mean per well) v1
# Step: 4.2 
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 20-08-2026

# load packages ####
library(data.table)
library(colorspace)
library(viridis)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(cluster)
library(igraph)
library(kBET)
# set variables ####
file_name <- "mPaintDR2_N2_N3_N4"
redu_state <- "redu"
integrate_state <- c("integrated", "unintegrated")
ctrl_cond <- c("DMSO_0", "CCCP_30", "ROT_10")
# create function to load data ####
load_data <- function(file_name, integrate_state) {
  # load profiles as df
  df <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_data_", integrate_state, "_", redu_state, ".csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames as WELL_BATCH
  rownames(df) <- df$V1
  df$V1 <- NULL
  # load integrated/unintegrated and redu/nonredu umap embeddings as umap
  umap <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_umap_embeddings.csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames as WELL_BATCH
  rownames(umap) <- umap$V1
  umap$V1 <- NULL
  # load metadata as meta
  meta <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_dimred_meta.csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames as WELL_BATCH
  rownames(meta) <- meta$V1
  meta$V1 <- NULL
  # remove any columns with NA/ non finite values
  umap <- umap[, colSums(!is.finite(as.matrix(umap))) == 0, drop = FALSE]
  # remove any rows with NA/ non finite values
  umap <- umap[apply(umap, 1, function(x) all(is.finite(x))), , drop = FALSE]
  # load umap nn edges
  nn <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_umap_nn_edges.csv", sep = ""), 
      header = TRUE)
  )
  # load umap knn matrix of indices
  knn <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_knn_umap.csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames
  rownames(knn) <- knn$V1
  knn$V1 <- NULL
  # return list of drift corrected data, raw data, metadata, and drift fits
  return(list(
    df = df,
    umap = umap,
    meta = meta,
    nn = nn,
    knn = knn
  ))
}
# run function to load data ####
# data is a large list containing sublists for integrate_state (integrated, unintegrated)
data <- integrate_state |>
  set_names() |>
  map(
    # each sublist contains corresponding umap and meta
    ~ load_data(
      file_name = file_name,
      integrate_state = .x
    )
  )
# create function to calculate asw ####
calc_asw <- function(embeddings, meta, group_var, keep_groups = NULL) {
  meta <- meta[rownames(embeddings), , drop = FALSE]
  if (!is.null(keep_groups)) {
    keep_idx <- meta[[group_var]] %in% keep_groups
    embeddings <- embeddings[keep_idx, , drop = FALSE]
    meta <- meta[keep_idx, , drop = FALSE]
  }
  group <- factor(meta[[group_var]])
  if (nlevels(group) < 2) {
    return(NA_real_)
  }
  sil <- cluster::silhouette(
    as.integer(group),
    dist(as.matrix(embeddings))
  )
  mean(sil[, 3], na.rm = TRUE)
}
# run function to calculate asw ####
asw_results <- purrr::imap_dfr(
  list(
    integrated = list(
      embeddings = data[["integrated"]][["umap"]],
      meta = data[["integrated"]][["meta"]]
    ),
    unintegrated = list(
      embeddings = data[["unintegrated"]][["umap"]],
      meta = data[["unintegrated"]][["meta"]]
    )
  ),
  function(obj, state_name) {
    tibble::tibble(
      state = state_name,
      asw_batch = calc_asw(obj$embeddings, obj$meta, "Batch"),
      asw_condition_ctrl = calc_asw(
        obj$embeddings,
        obj$meta,
        "Condition",
        keep_groups = ctrl_cond
      )
    )
  }
)
# create function to calculate gc ####
calc_gc <- function(nn,
                    meta,
                    label_var = "Batch",
                    keep_groups = NULL) {
  
  # keep rows aligned
  meta <- meta[rownames(meta) %in% unique(c(nn$from, nn$to)), , drop = FALSE]
  
  # optional subset, e.g. only control conditions
  if (!is.null(keep_groups)) {
    keep_idx <- meta[[label_var]] %in% keep_groups
    meta <- meta[keep_idx, , drop = FALSE]
  }
  
  # keep only edges where both endpoints are present after subsetting
  nn <- nn[nn$from %in% rownames(meta) & nn$to %in% rownames(meta), , drop = FALSE]
  
  # build graph
  g <- igraph::graph_from_data_frame(
    d = nn[, c("from", "to"), drop = FALSE],
    directed = FALSE,
    vertices = data.frame(name = rownames(meta))
  )
  
  # simplify to avoid duplicate edges / loops
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  labels <- factor(meta[[label_var]])
  cell_ids <- rownames(meta)
  
  if (nlevels(labels) < 2) {
    return(NA_real_)
  }
  
  # fraction of each label in its largest connected component
  lcc_frac <- vapply(levels(labels), function(lv) {
    cells <- cell_ids[labels == lv]
    
    if (length(cells) < 2) {
      return(1)
    }
    
    subg <- igraph::induced_subgraph(g, vids = cells)
    comps <- igraph::components(subg)
    max(comps$csize) / length(cells)
  }, numeric(1))
  
  mean(lcc_frac, na.rm = TRUE)
}
# run function to calculate gc ####
gc_results <- purrr::imap_dfr(
  data,
  function(obj, state_name) {
    tibble::tibble(
      state = state_name,
      gc_batch = calc_gc(
        nn = obj$nn,
        meta = obj$meta,
        label_var = "Batch"
      ),
      gc_condition_ctrl = calc_gc(
        nn = obj$nn,
        meta = obj$meta,
        label_var = "Condition",
        keep_groups = ctrl_cond
      ),
      gc_condition = calc_gc(
        nn = obj$nn,
        meta = obj$meta,
        label_var = "Condition"
      )
    )
  }
)

# create function to calculate kbet ####
calc_kbet <- function(umap, meta, knn, batch_var = "Batch", k0 = NULL) {
  meta <- meta[rownames(umap), , drop = FALSE]
  
  batch <- meta[[batch_var]]
  
  # ensure knn is a plain matrix of indices
  knn <- as.matrix(knn)
  
  if (is.null(k0)) {
    k0 <- ncol(knn)
  }
  
  res <- kBET::kBET(
    df = as.matrix(umap),
    batch = batch,
    k0 = k0,
    knn = knn,
    do.pca = FALSE,
    plot = FALSE,
    verbose = FALSE
  )
  
  tibble::tibble(
    kbet_rejection_rate = res$summary$kBET.observed[1],
    kbet_expected_rejection_rate = res$summary$kBET.expected[1],
    average_pval = res$average.pval
  )
}
# run function to calculate kbet ####
kbet_results <- purrr::imap_dfr(
  data,
  function(obj, state_name) {
    calc_kbet(
      umap = obj$umap,
      meta = obj$meta,
      knn = obj$knn,
      batch_var = "Batch"
    ) |>
      dplyr::mutate(state = state_name, .before = 1)
  }
)
# create a function to calculate map ####
cosine_similarity_mat <- function(x) {
  x <- as.matrix(x)
  norms <- sqrt(rowSums(x^2))
  norms[norms == 0] <- NA_real_
  x <- x / norms
  sim <- x %*% t(x)
  diag(sim) <- NA_real_  # exclude self-match
  sim
}

average_precision_one <- function(relevant_logical) {
  relevant_logical <- as.logical(relevant_logical)
  n_pos <- sum(relevant_logical, na.rm = TRUE)
  if (n_pos == 0) return(NA_real_)
  
  prec_at_k <- cumsum(relevant_logical) / seq_along(relevant_logical)
  mean(prec_at_k[relevant_logical], na.rm = TRUE)
}

calc_map <- function(df, meta, label_var = "Condition", keep_groups = NULL) {
  meta <- meta[rownames(df), , drop = FALSE]
  
  if (!is.null(keep_groups)) {
    keep_idx <- meta[[label_var]] %in% keep_groups
    df <- df[keep_idx, , drop = FALSE]
    meta <- meta[keep_idx, , drop = FALSE]
  }
  
  labels <- factor(meta[[label_var]])
  
  if (nlevels(labels) < 2) {
    return(tibble(
      mAP = NA_real_,
      n_queries = nrow(df)
    ))
  }
  
  sim <- cosine_similarity_mat(df)
  
  ap_by_query <- vapply(seq_len(nrow(sim)), function(i) {
    ord <- order(sim[i, ], decreasing = TRUE, na.last = NA)
    ord <- ord[ord != i]
    
    if (length(ord) == 0) return(NA_real_)
    
    relevant <- labels[ord] == labels[i]
    average_precision_one(relevant)
  }, numeric(1))
  
  tibble(
    mAP = mean(ap_by_query, na.rm = TRUE),
    n_queries = sum(!is.na(ap_by_query))
  )
}
# run function to calculate map ####
map_results <- imap_dfr(
  list(
    integrated = data[["integrated"]],
    unintegrated = data[["unintegrated"]]
  ),
  function(obj, state_name) {
    calc_map(
      df = obj$df,
      meta = obj$meta,
      label_var = "Condition",
      keep_groups = ctrl_cond
    ) |>
      mutate(state = state_name, .before = 1)
  }
)

# plot all_results ####
all_results <- dplyr::bind_rows(
  asw_results  |> tidyr::pivot_longer(-state, names_to = "metric", values_to = "value"),
  gc_results   |> tidyr::pivot_longer(-state, names_to = "metric", values_to = "value"),
  kbet_results |> tidyr::pivot_longer(-state, names_to = "metric", values_to = "value"),
  map_results  |> tidyr::pivot_longer(-state, names_to = "metric", values_to = "value")
) |>
  tidyr::pivot_wider(
    names_from = state,
    values_from = value
  )

plot_metric_grid <- function(results_df) {
  
  metric_order <- c(
    "gc_batch",
    "kbet_rejection_rate",
    "asw_batch",
    "gc_condition_ctrl",
    "gc_condition",
    "asw_condition_ctrl",
    "mAP"
  )
  
  metric_labels <- c(
    "Graph Connectivity (Batch)",
    "K-Nearest Neighbour Batch Effect Test",
    "Silhouette (Batch)",
    "Graph Connectivity (Control Conditions)",
    "Graph Connectivity (All Conditions)",
    "Silhouette (Control Conditions)",
    "mAP (Control Conditions)"
  )
  
  state_labels <- c(
    "Unintegrated",
    "Integrated"
  )
  
  plot_df <- results_df |>
    dplyr::filter(metric %in% metric_order) |>
    dplyr::mutate(
      metric = factor(
        metric,
        levels = rev(metric_order),
        labels = rev(metric_labels)
      )
    ) |>
    tidyr::pivot_longer(
      cols = c(unintegrated, integrated),
      names_to = "state",
      values_to = "value"
    ) |>
    dplyr::mutate(
      state = factor(
        state,
        levels = c("unintegrated", "integrated"),
        labels = state_labels
      ),
      label = sprintf("%.2f", value)
    )
  
  ggplot2::ggplot(plot_df, ggplot2::aes(x = state, y = metric)) +
    ggplot2::geom_point(
      ggplot2::aes(fill = value),
      shape = 21,
      size = 14,
      colour = "black",
      stroke = 0.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      colour = "black",
      size = 3
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = lighten("steelblue", amount = 0.2),
      na.value = "white"
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::coord_equal() +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(
        size = 10,
        angle = 45,
        hjust = 0,
        face = "bold"
      ),
      legend.position = "right",
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    ) +
    ggplot2::labs(fill = "Metric value") +
    theme(
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
}

plot_metric_grid(all_results)

# save data ####
# save plot ####