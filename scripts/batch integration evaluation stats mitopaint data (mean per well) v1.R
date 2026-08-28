# Title: batch integration evaluation stats mitopaint data (mean per well) v1
# Step: 4.2 
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 28-08-2026

# load packages ####
library(data.table)
library(colorspace)
library(viridis)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(cluster)
library(mclust)
library(igraph)
library(kBET)
library(lisi)
library(aricode)
# set variables ####
file_name <- "mPaintDR2_N2_N3_N4"
redu_state <- "redu"
integrate_states <- c("integrated", "unintegrated")
ctrl_cond <- c("DMSO_0", "CCCP_30", "ROT_10")
# create function to load data ####
load_data <- function(file_name, integrate_state) {
  # avoid bug with numeric coercion with as.numeric helper function
  as_numeric_df <- function(x) {
    x <- as.data.frame(x)
    x[] <- lapply(x, function(col) as.numeric(as.character(col)))
    x
  }
  # set paths
  df_path   <- paste0("data/processed/", file_name, "_data_", integrate_state, "_", redu_state, ".csv")
  umap_path <- paste0("data/processed/", file_name, "_", integrate_state, "_", redu_state, "_umap_embeddings.csv")
  meta_path <- paste0("data/processed/", file_name, "_", integrate_state, "_", redu_state, "_dimred_meta.csv")
  nn_path   <- paste0("data/processed/", file_name, "_", integrate_state, "_", redu_state, "_umap_nn_edges.csv")
  knn_path  <- paste0("data/processed/", file_name, "_", integrate_state, "_", redu_state, "_knn_umap.csv")
  # load profiles as df
  # df is needed for mAP
  df <- as.data.frame(fread(df_path, header = TRUE))
  rownames(df) <- df$V1
  df$V1 <- NULL
  df <- as_numeric_df(df)
  df <- df[, colSums(!is.finite(as.matrix(df))) == 0, drop = FALSE]
  df <- df[apply(df, 1, function(x) all(is.finite(x))), , drop = FALSE]
  # load umap embeddings as umap
  # umap is needed for kbet, silhouette width, and local inverse simpsons index metrics
  umap <- as.data.frame(fread(umap_path, header = TRUE))
  rownames(umap) <- umap$V1
  umap$V1 <- NULL
  umap <- as_numeric_df(umap)
  # remove any columns with NA/ non finite values
  # there ordinarily shouldnt be any NAs
  # check nrows df and umap to double check concordancy
  umap <- umap[, colSums(!is.finite(as.matrix(umap))) == 0, drop = FALSE]
  umap <- umap[apply(umap, 1, function(x) all(is.finite(x))), , drop = FALSE]
  # load metadata as meta
  meta <- as.data.frame(fread(meta_path, header = TRUE))
  rownames(meta) <- meta$V1
  meta$V1 <- NULL
  # load umap nn edges 
  # nn edges needed for graph connectivity and leiden metrics
  nn <- as.data.frame(fread(nn_path, header = TRUE))
  nn$from <- as.character(nn$from)
  nn$to <- as.character(nn$to)
  # load umap knn matrix of indices 
  # knn matrix needed for kbet metric
  knn <- as.data.frame(fread(knn_path, header = TRUE))
  rownames(knn) <- knn$V1
  knn$V1 <- NULL
  knn <- as_numeric_df(knn)
  # keep only rows shared across df/umap/meta if needed
  common_ids <- Reduce(intersect, list(rownames(df), rownames(umap), rownames(meta)))
  df <- df[common_ids, , drop = FALSE]
  umap <- umap[common_ids, , drop = FALSE]
  meta <- meta[common_ids, , drop = FALSE]
  # align knn to the same cells if possible
  common_knn_ids <- intersect(rownames(knn), common_ids)
  knn <- knn[common_knn_ids, , drop = FALSE]
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
data <- integrate_states |>
  set_names() |>
  purrr::imap(~ tryCatch(load_data(file_name = file_name, integrate_state = .x),
                         error = function(e) NULL))
# create function to calculate asw ####
calc_asw <- function(embeddings, meta, group_var, keep_groups = NULL) {
  # match row order of meta and umap embeddings
  meta <- meta[rownames(embeddings), , drop = FALSE]
  # if condition, subset meta and umap embeddings based (eg. ctrl_cond)
  if (!is.null(keep_groups)) {
    keep_idx <- meta[[group_var]] %in% keep_groups
    embeddings <- embeddings[keep_idx, , drop = FALSE]
    meta <- meta[keep_idx, , drop = FALSE]
  }
  # if only one group supplied, asw cannot be computed so returns NA
  group <- factor(meta[[group_var]])
  if (nlevels(group) < 2) {
    return(NA_real_)
  }
  # calculate average silhouette width by on grouping variable
  # asw is the mean of the pairwise distance between embedding points
  sil <- cluster::silhouette(
    as.integer(group),
    dist(as.matrix(embeddings))
  )
  mean(sil[, 3], na.rm = TRUE)
}
# run function to calculate asw ####
asw_results <- purrr::imap_dfr(
  # imap loop calculates asw for integrated and unintegrated data
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
  # computed asw compiled into a data table
  function(obj, state_name) {
    tibble::tibble(
      # column for state (integrated or unintegrated)
      state = state_name,
      # column for batch grouping
      # how does the spread between embeddings compare between integrated and unintegrated
      # lower asw_batch for integrated data means data is better mixed
      asw_batch = calc_asw(obj$embeddings, obj$meta, "Batch"),
      # column for control condition grouping
      # higher asw_cond_ctrl for integrated data means meaningful biological variance is better preserved
      asw_condition_ctrl = calc_asw(obj$embeddings, obj$meta, "Condition", keep_groups = ctrl_cond)
    )
  }
)
# create function to calculate gc ####
calc_gc <- function(nn,
                    meta,
                    label_var = "Batch",
                    keep_groups = NULL) {
  # keep only meta rows which are present in nn edge list
  meta <- meta[rownames(meta) %in% unique(c(nn$from, nn$to)), , drop = FALSE]
  # optional subset (eg control conditions)
  if (!is.null(keep_groups)) {
    keep_idx <- meta[[label_var]] %in% keep_groups
    meta <- meta[keep_idx, , drop = FALSE]
    }
  # keep only edges where both endpoints are present after subsetting
  nn <- nn[nn$from %in% rownames(meta) & nn$to %in% rownames(meta), , drop = FALSE]
  # build an igraph object from nn edge list
  g <- igraph::graph_from_data_frame(
    d = nn[, c("from", "to"), drop = FALSE],
    directed = FALSE,
    vertices = data.frame(name = rownames(meta))
  )
  # simplify to avoid duplicate edges / loops
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  # save labels and cell ids
  labels <- factor(meta[[label_var]])
  cell_ids <- rownames(meta)
  # if there is only one group, graph connectivity is not meaningful so returns NA
  if (nlevels(labels) < 2) {
    return(NA_real_)
  }
  # loops over all labels (eg. each batch, each condition)
  lcc_frac <- vapply(levels(labels), function(lv) {
    cells <- cell_ids[labels == lv]
    # if only one cell in a given label it is "connected"
    if (length(cells) < 2) {
      return(1)
    }
    # computes fraction of label in connected component
    # graph connectivity is the average of these fractions
    # the closer the number is to one, the more connected the points are together 
    # high graph connectivity by batch = points are separating by batch, so poor batch mixing
    # high graph connectivity by condition = points are separating by condition, so preserved real biological variance
    subg <- igraph::induced_subgraph(g, vids = cells)
    comps <- igraph::components(subg)
    max(comps$csize) / length(cells)
  }, numeric(1))
  mean(lcc_frac, na.rm = TRUE)
}
# run function to calculate gc ####
# loop calc_gc over integrated and unintegrated data and combines results in one table
gc_results <- purrr::imap_dfr(
  data,
  function(obj, state_name) {
    tibble::tibble(
      state = state_name,
      # provide variables for gc by batch
      gc_batch = calc_gc(
        nn = obj$nn,
        meta = obj$meta,
        label_var = "Batch"
      ),
      # provide variables for gc by control condition
      gc_condition_ctrl = calc_gc(
        nn = obj$nn,
        meta = obj$meta,
        label_var = "Condition",
        keep_groups = ctrl_cond
      ),
      # provide variables for gc by ALL conditions (no subsetting beforehand)
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
  # align metadata and umap embeddings rownames
  meta <- meta[rownames(umap), , drop = FALSE]
  # pull batch labels
  batch <- meta[[batch_var]]
  # ensure knn is a plain matrix of indices
  knn <- as.matrix(knn)
  # k0 is the number of columns in knn if not explicitly provided 
  # k is usually ~15 (check variables used in 6. dimensionality reduction mitopaint data (mean per well) v1 script)
  if (is.null(k0)) {
    k0 <- ncol(knn)
  }
  # run kBET (k-nearest neighbour batch effect test) on umap embedding and batch labels
  res <- kBET::kBET(
    df = as.matrix(umap),
    batch = batch,
    k0 = k0,
    knn = knn,
    do.pca = FALSE,
    plot = FALSE,
    verbose = FALSE
  )
  # returns a table of the kbet rejection rate, kbet expected rejection rate and average p value
  # kbet rejection rate is the fraction of tested neighborhoods rejects the hypothesis that the batch composition within the neigborhood matches the global batch composition
  # lower kbet rejection rate means less of the neighborhoods are enriched for a batch (better batch mixing)
  # expected rejection rate is assuming random distribution how often would the hypothesis get rejected? well mixed = kbet rate is close to expected
  tibble::tibble(
    kbet_rejection_rate = res$summary$kBET.observed[1],
    kbet_expected_rejection_rate = res$summary$kBET.expected[1],
    # average pval asks across all neigborhoods are batches not mixed well (p < 0.05 = batches are very poorly mixed)
    average_pval = res$average.pval
  )
}
# run function to calculate kbet ####
kbet_results <- purrr::imap_dfr(
  data,
  function(obj, state_name) {
    # run kbet on both integrated and unintegrated data
    calc_kbet(
      umap = obj$umap,
      meta = obj$meta,
      knn = obj$knn,
      batch_var = "Batch"
    ) |>
      # add state_name column if integrated/ unintegrated
      dplyr::mutate(state = state_name, .before = 1)
  }
)
# create a function to calculate map ####
# helper function to compute pairwise cosine similarity matrix between rows
calc_cosine_sim <- function(x) {
  x <- as.matrix(x)
  norms <- sqrt(rowSums(x^2))
  norms[norms == 0] <- NA_real_
  x <- x / norms
  sim <- x %*% t(x)
  diag(sim) <- NA_real_  # exclude self-match
  sim
}
# helper function to calculate average precision for one condition
calc_ap <- function(relevant_logical) {
  relevant_logical <- as.logical(relevant_logical)
  # counts how many replicates of a condition there are
  n_pos <- sum(relevant_logical, na.rm = TRUE)
  # will not compute it there are no replicates
  if (n_pos == 0) return(NA_real_)
  # counts how similar each replicate of a condition is to each other, then averages
  prec_at_k <- cumsum(relevant_logical) / seq_along(relevant_logical)
  mean(prec_at_k[relevant_logical], na.rm = TRUE)
}
# mean average precision is the the average precision for each condition, averaged altogether
# it asks how similar is the profile of each replicate of the same condition 
calc_map <- function(df, meta, label_var = "Condition", keep_groups = NULL) {
  # align rownames of metadata to profile matrix
  meta <- meta[rownames(df), , drop = FALSE]
  # optional subsetting step (eg positive controls)
  if (!is.null(keep_groups)) {
    keep_idx <- meta[[label_var]] %in% keep_groups
    df <- df[keep_idx, , drop = FALSE]
    meta <- meta[keep_idx, , drop = FALSE]
  }
  # if there are less than two unique conditions mAP cant be computed
  labels <- factor(meta[[label_var]])
  if (nlevels(labels) < 2) {
    return(tibble(
      mAP = NA_real_,
      n_queries = nrow(df)
    ))
  }
  # compute similarity matrix from profiling data using helper function
  sim <- calc_cosine_sim(df)
  # loop over each well from each condition, ranking which are most similar to the query well
  ap_by_query <- vapply(seq_len(nrow(sim)), function(i) {
    ord <- order(sim[i, ], decreasing = TRUE, na.last = NA)
    # skip comparing the same well to itself
    ord <- ord[ord != i]
    # skip computing if only one well replicate of condition
    if (length(ord) == 0) return(NA_real_)
    # collate all replicate wells together
    relevant <- labels[ord] == labels[i]
    # compute the precision (cosine similarity ranking) of all replicate wells
    calc_ap(relevant)
  }, numeric(1))
  # collate mean of the average precisions
  tibble(
    mAP = mean(ap_by_query, na.rm = TRUE),
    n_queries = sum(!is.na(ap_by_query))
  )
}
# run function to calculate map ####
# loop over both integrated and unintegrated data
map_results <- imap_dfr(
  list(
    integrated = data[["integrated"]],
    unintegrated = data[["unintegrated"]]
  ),
  # calculate mAP only on positive control conditions 
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
# create function to calculate lisi ####
calc_lisi <- function(embeddings,
                      meta,
                      label_var,
                      keep_groups = NULL,
                      # ideally match perplexity to that used in `6.` dimensionality reduction mitopaint data (mean per well) v1 script
                      perplexity = 30) {
  # align metadata rows and umap embedding rows
  meta <- meta[rownames(embeddings), , drop = FALSE]
  # optional subset (eg positive controls)
  if (!is.null(keep_groups)) {
    keep_idx <- meta[[label_var]] %in% keep_groups
    embeddings <- embeddings[keep_idx, , drop = FALSE]
    meta <- meta[keep_idx, , drop = FALSE]
  }
  # compute lisi using lisi package
  lisi_df <- lisi::compute_lisi(
    X = as.matrix(embeddings),
    meta_data = meta,
    label_colnames = label_var,
    # perplexity is number of each wells/ each cells neighbors
    perplexity = perplexity
  )
  # return mean and sd lisi in a table
  # local inverse simpsons index is a metric of how diverse a given cell/ wells neighbors are
  # higher lisi by batch = more batch mixing, lower lisi by condition = preserve biologically meaningful variance between conditions (ie umap clusters are specific for a condition)
  tibble::tibble(
    lisi_mean = mean(lisi_df[[label_var]], na.rm = TRUE),
    lisi_sd = sd(lisi_df[[label_var]], na.rm = TRUE)
  )
}
# run function to calculate lisi ####
lisi_results <- purrr::imap_dfr(
  # loop lisi calculation for both integrated and unintegrated data
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
  # calculate lisi by batch
  function(obj, state_name) {
    batch_lisi <- calc_lisi(
      obj$embeddings,
      obj$meta,
      label_var = "Batch"
    )
    # calculate lisi by positive control condition
    condition_lisi <- calc_lisi(
      obj$embeddings,
      obj$meta,
      label_var = "Condition",
      keep_groups = ctrl_cond
    )
    # compile all lisi metrics into a single data table
    tibble::tibble(
      state = state_name,
      lisi_batch_mean = batch_lisi$lisi_mean,
      lisi_batch_sd = batch_lisi$lisi_sd,
      lisi_condition_ctrl_mean = condition_lisi$lisi_mean,
      lisi_condition_ctrl_sd = condition_lisi$lisi_sd
    )
  }
)
# create function to calculate leiden ####
# define function to run leiden clustering and compare to true condition annotations
calc_leiden_ari_nmi <- function(nn, meta, label_var = "Compound", keep_groups = NULL) {
  # clean edge list and remove duplicates
  nn <- nn |>
    dplyr::mutate(
      from = as.character(from),
      to   = as.character(to)
    ) |>
    dplyr::filter(
      !is.na(from), from != "",
      !is.na(to),   to != ""
    ) |>
    dplyr::distinct(from, to)
  # keep only metadata of wells present in nn edge list
  meta <- meta[rownames(meta) %in% unique(c(nn$from, nn$to)), , drop = FALSE]
  # optional subsetting (keep only positive control conditions)
  if (!is.null(keep_groups)) {
    meta <- meta[meta[[label_var]] %in% keep_groups, , drop = FALSE]
  }
  # keep only nn edge present in subsetted meta
  cells <- rownames(meta)
  cells <- cells[!is.na(cells) & cells != ""]
  nn <- nn |>
    dplyr::filter(from %in% cells, to %in% cells)
  # if there is only one well in edge list, leiden clusters cant be computed so return NA
  if (nrow(nn) == 0 || length(cells) < 2) {
    return(tibble::tibble(
      leiden_ari = NA_real_,
      leiden_nmi = NA_real_,
      n_cells = length(cells),
      n_leiden_clusters = NA_integer_
    ))
  }
  # convert edge list into a matrix
  edge_mat <- as.matrix(nn[, c("from", "to")])
  storage.mode(edge_mat) <- "character"
  # build igraph object, removing repeated edges
  g <- igraph::graph_from_edgelist(edge_mat, directed = FALSE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  # run leiden clusters on graph
  leiden <- igraph::cluster_leiden(g, resolution = 0.01)
  leiden_lab <- factor(igraph::membership(leiden))
  # match meta labels to graph vertices actually present
  common <- intersect(names(igraph::membership(leiden)), rownames(meta))
  leiden_lab <- factor(igraph::membership(leiden)[common])
  ref_lab <- factor(meta[common, label_var])
  # return leiden adjusted rand index and leiden normalised mutual information
  # ARI = how often pairs of cells/ wells are in the same cluster relative to chance
  # the closer to 1, the closer the cells/ wells are a perfect match
  # NMI = how much information the clustering shares with the reference labels
  # the closer to 1, to closer the metadata labels are an exact overlap to the leiden clusters
  # both are two ways to measure agreement between leiden clusters and actual metadata labels
  # it asks: do the leiden clusters reflect the expected similarities within/ differences between conditions 
  tibble::tibble(
    leiden_ari = mclust::adjustedRandIndex(leiden_lab, ref_lab),
    leiden_nmi = aricode::NMI(leiden_lab, ref_lab),
    n_cells = length(common),
    n_leiden_clusters = nlevels(leiden_lab)
  )
}
# run function to calculate leiden ####
leiden_results <- purrr::imap_dfr(
  # loop over integrated and unintegrated data
  data,
  function(obj, state_name) {
    # calculates leiden clustering on control conditions only
    calc_leiden_ari_nmi(
      nn = obj$nn,
      meta = obj$meta,
      label_var = "Condition",
      keep_groups = ctrl_cond
    ) |>
      # add column for state name (either integrated or unintegrated)
      mutate(state = state_name, .before = 1)
  }
)
# collate all_results ####
# combine all batch integration stat metrics into a single table
all_results <- dplyr::bind_rows(
  asw_results  |> dplyr::select(state, asw_batch, asw_condition_ctrl) |> tidyr::pivot_longer(-state, names_to = "metric", values_to = "value"),
  gc_results   |> dplyr::select(state, gc_batch, gc_condition_ctrl, gc_condition) |> tidyr::pivot_longer(-state, names_to = "metric", values_to = "value"),
  kbet_results |> dplyr::select(state, kbet_rejection_rate) |> tidyr::pivot_longer(-state, names_to = "metric", values_to = "value"),
  lisi_results |> dplyr::select(state, lisi_batch_mean, lisi_condition_ctrl_mean) |> tidyr::pivot_longer(-state, names_to = "metric", values_to = "value"),
  leiden_results |> dplyr::select(state, leiden_ari, leiden_nmi) |> tidyr::pivot_longer(-state, names_to = "metric", values_to = "value"),
  map_results  |> dplyr::select(state, mAP) |> tidyr::pivot_longer(-state, names_to = "metric", values_to = "value")
) |>
  # one row per metric, and one column per state
  tidyr::pivot_wider(
    names_from = state,
    values_from = value
  )
# plot all_results ####
# function to plot all batch correction metrics into single visualisation
plot_metric_grid <- function(results_df) {
  # manual ordering of metrics
  metric_order <- c(
    "gc_batch",
    "kbet_rejection_rate",
    "lisi_batch_mean",
    "asw_batch",
    "lisi_condition_ctrl_mean",
    "leiden_ari",
    "leiden_nmi",
    "gc_condition_ctrl",
    "gc_condition",
    "asw_condition_ctrl",
    "mAP"
  )
  # manually set labels of metrics
  metric_labels <- c(
    "GC (Batch)",
    "KBET",
    "LISI (Batch)",
    "ASW (Batch)",
    "LISI (Controls)",
    "Leiden ARI (Controls)",
    "Leiden NMI (Controls)",
    "GC (Controls)",
    "GC (All Conditions)",
    "ASW (Controls)",
    "mAP (Controls)"
  )
  # plot_df is table of batch correction metrics, reordered
  plot_df <- results_df |>
    dplyr::filter(metric %in% metric_order) |>
    dplyr::mutate(
      metric = factor(
        metric,
        levels = rev(metric_order),
        labels = rev(metric_labels)
      )
    ) |>
    # each metric is on one row (each column in either unintegrated or integrated data)
    tidyr::pivot_longer(
      cols = c(unintegrated, integrated),
      names_to = "state",
      values_to = "value"
    ) |>
    # order columns (unintegrated, integrated)
    dplyr::mutate(
      state = factor(
        state,
        levels = c("unintegrated", "integrated")
      ),
      # metric value to plot is 2dp
      state_pos = as.numeric(state),
      metric_pos = as.numeric(metric),
      label = sprintf("%.2f", value)
    ) |>
    # binary color by row (larger value is colored)
    dplyr::group_by(metric) |>
    dplyr::mutate(
      fill_value = scales::rescale(
        value,
        to = c(0, 1),
        from = range(value, na.rm = TRUE)
      )
    ) |>
    dplyr::ungroup()
  # plot figure
  ggplot2::ggplot(plot_df, ggplot2::aes(x = state_pos, y = metric_pos)) +
    # geom_point (circle with black outline)
    ggplot2::geom_point(
      ggplot2::aes(fill = fill_value),
      shape = 21,
      size = 14,
      colour = "black",
      stroke = 0.5
    ) +
    # text inside circle of metric value (2dp)
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      colour = "black",
      size = 3
    ) +
    # annotation lines on the side say whether metric is batch correction metric or biometric
    # see Arevelo et al., 2024 Nat Comm for more information
    annotate("segment", x = 2.9, xend = 2.9, y = 7.7, yend = 11.5, linewidth = 0.5) +
    annotate("text",    x = 3.18, y = 9.5, label = "Batch correction",
             angle = 90, hjust = 0.5, vjust = 0.5, size = 3.5) +
    annotate("segment", x = 2.9, xend = 2.9, y = 0.5, yend = 7.3, linewidth = 0.5) +
    annotate("text",    x = 3.18, y = 4.0, label = "Bio metrics",
             angle = 90, hjust = 0.5, vjust = 0.5, size = 3.5) +
    # fill color is blue (higher metric value in row is blue)
    ggplot2::scale_fill_gradient(
      low = "white",
      high = colorspace::lighten("steelblue", amount = 0.4),
      na.value = "white",
      limits = c(0, 1),
      guide = "none"
    ) +
    # tidy dims to avoid overlaps
    scale_x_continuous(
      breaks = c(1, 2),
      labels = c("Unintegrated", "Integrated"),
      position = "top",
      expand = expansion(mult = c(0.25, 0.25))
    ) +
    ggplot2::scale_y_continuous(
      breaks = 1:11,
      labels = rev(metric_labels),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::coord_equal(clip = "off") +
    # tidy theme
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
      plot.margin = ggplot2::margin(5, 5, 5, 5),
      legend.position = "none"
    )
}
plot <- plot_metric_grid(all_results)
plot
# save data ####
write.csv(
  all_results,
  paste("outputs/data/", file_name, "_", redu_state, "_batchcorr_stats.csv", sep = "")
)
# save plot ####
ggsave(
  paste0("outputs/figures/", file_name, "_", redu_state, "_batchcorr_stats.pdf", sep = ""),
  plot = plot,
  width = 3,
  height = 6.5,
  units = "in",
  dpi = 300
)
rm(list = ls())