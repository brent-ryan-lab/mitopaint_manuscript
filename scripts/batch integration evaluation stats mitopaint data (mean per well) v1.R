# Title: batch integration evaluation stats mitopaint data (mean per well) v1
# Step: 4.2 
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 20-08-2026

# load packages ####
library(data.table)
#library(colorspace)
library(tidyverse)
library(Seurat)
#library(ggplot2)
#library(ggpubr)
library(cluster)
# set variables ####
file_name <- "mPaintSpace2_N1_N2_N3"
redu_state <- "redu"
integrate_state <- c("integrated", "unintegrated")
ctrl_cond <- c("DMSO_0", "CCCP_20", "ROT_3")
# create function to load data ####
load_data <- function(file_name, integrate_state) {
  # load integrated/unintegrated and redu/nonredu data as df
  df <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_umap_embeddings.csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames as WELL_BATCH
  rownames(df) <- df$V1
  df$V1 <- NULL
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
  df <- df[, colSums(!is.finite(as.matrix(df))) == 0, drop = FALSE]
  # remove any rows with NA/ non finite values
  df <- df[apply(df, 1, function(x) all(is.finite(x))), , drop = FALSE]
  # return list of drift corrected data, raw data, metadata, and drift fits
  return(list(
    df = df,
    meta = meta
  ))
}
# run function to load data ####
# data is a large list containing sublists for integrate_state (integrated, unintegrated)
data <- integrate_state |>
  set_names() |>
  map(
    # each sublist contains corresponding df and meta
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
      embeddings = data[["integrated"]][["df"]],
      meta = data[["integrated"]][["meta"]]
    ),
    unintegrated = list(
      embeddings = data[["unintegrated"]][["df"]],
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
# create function to calculate graph connectivity
calc_gc <- 