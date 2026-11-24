# Title: batch integration mitopaint data (mean per well) v2
# Step: 4
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 05-11-2026

# load packages ####
library(data.table)
library(Seurat)
library(purrr)
# set variables ####
batches_info <- list(
  N1 = list(
    file_name = "SF240627_mPaintDR2_N2"
  ),
  N2 = list(
    file_name = "SF240704_mPaintDR2_N3"
  ),
  N3 = list(
    file_name = "SF240711_mPaintDR2_N4"
  )
)
file_name <- "mPaintDR2_N2_N3_N4"
k_weight <- 50
# create a function to load data ####
load_data <- function(file_name) {
  # load zscore data as df
  df <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_z.csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames as WELL_BATCH
  rownames(df) <- df$V1
  df$V1 <- NULL
  # load metadata as meta
  meta <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_meta_tidy.csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames as WELL_BATCH
  rownames(meta) <- meta$V1
  meta$V1 <- NULL
  # return list of drift corrected data, raw data, metadata, and drift fits
  return(list(
    df = df,
    meta = meta
  ))
}
# run function to load data ####
batches <- purrr::map(
  batches_info,
  function(batch_info) {
    load_data(
      file_name = batch_info$file_name
    )
  }
)
# create function to convert one batch to a Seurat object ####
make_seurat_object <- function(batch_obj) {
  # fetch zscore and metadata
  z_mat <- batch_obj$df
  meta_df <- batch_obj$meta
  # transpose zscore matrix
  # Seurat expects features as rows and wells/samples as columns
  mat <- t(as.matrix(z_mat))
  # create seurat object for each batch
  seurat_obj <- CreateSeuratObject(
    counts = mat,
    assay = "MP",
    meta.data = meta_df,
    min.cells = 0,
    min.features = 0
  )
  # manually copy counts to data layer to circumvent v5 seurat structure syntax
  seurat_obj <- SetAssayData(
    seurat_obj,
    assay = "MP",
    layer = "data",
    new.data = GetAssayData(
      seurat_obj,
      assay = "MP",
      layer = "counts"
    )
  )
  # scale data
  seurat_obj <- ScaleData(
    seurat_obj,
    verbose = FALSE
  )
  return(seurat_obj)
}
# put unintegrated data in a single table ####
unintegrated <- list()
# all unintegrated batch data
unintegrated$df <- purrr::map(
  batches,
  "df"
) |>
  dplyr::bind_rows()
# all unintegrated metadata
unintegrated$meta <- purrr::map(
  batches,
  "meta"
) |>
  dplyr::bind_rows()
# integrate data only when there is more than one batch ####
if (length(batches) > 1) {
  
  message(
    length(batches),
    " batches detected: running batch integration."
  )
  
  # convert each batch to a Seurat object
  df.seurat <- purrr::map(
    batches,
    make_seurat_object
  )
  
  # find integration anchors
  anchors <- FindIntegrationAnchors(
    object.list = df.seurat,
    # only integrate features shared across all batches
    anchor.features = Reduce(
      intersect,
      purrr::map(
        batches,
        ~ colnames(.x$df)
      )
    ),
    scale = FALSE
  )
  
  # integrate batches
  integrated.seurat <- IntegrateData(
    anchorset = anchors,
    k.weight = k_weight
  )
  
  # put integrated data in a single table
  integrated <- list()
  
  integrated$df <- as.data.frame(
    t(
      GetAssayData(
        integrated.seurat,
        assay = "integrated",
        layer = "data"
      )
    )
  )
  
  integrated$meta <- integrated.seurat@meta.data
  
  # safety check
  stopifnot(
    identical(
      rownames(integrated$df),
      rownames(integrated$meta)
    )
  )
  
} else {
  
  message(
    "Only one batch detected: skipping batch integration."
  )
}
# create function to convert one batch to a Seurat object ####
make_seurat_object <- function(batch_obj) {
  # fetch zscore and metadata
  z_mat <- batch_obj$df
  meta_df <- batch_obj$meta
  # transpose zscore matrix
  # Seurat expects features as rows and wells/samples as columns
  mat <- t(as.matrix(z_mat))
  # create seurat object for each batch
  seurat_obj <- CreateSeuratObject(
    counts = mat,
    assay = "MP",
    meta.data = meta_df,
    min.cells = 0,
    min.features = 0
  )
  # manually copy counts to data layer to circumvent v5 seurat structure syntax
  seurat_obj <- SetAssayData(
    seurat_obj,
    assay = "MP",
    layer = "data",
    new.data = GetAssayData(
      seurat_obj,
      assay = "MP",
      layer = "counts"
    )
  )
  # scale data
  seurat_obj <- ScaleData(
    seurat_obj,
    verbose = FALSE
  )
  return(seurat_obj)
}
# integrate data only when there is more than one batch ####
if (length(batches) > 1) {
  message(
    length(batches),
    " batches detected: running batch integration."
  )
  # convert each batch to a Seurat object
  df.seurat <- purrr::map(
    batches,
    make_seurat_object
  )
  # find integration anchors
  anchors <- FindIntegrationAnchors(
    object.list = df.seurat,
    # only integrate features shared across all batches
    anchor.features = Reduce(
      intersect,
      purrr::map(
        batches,
        ~ colnames(.x$df)
      )
    ),
    scale = FALSE
  )
  # integrate batches
  integrated.seurat <- IntegrateData(
    anchorset = anchors,
    k.weight = k_weight
  )
  # put integrated data in a single table
  integrated <- list()
  integrated$df <- as.data.frame(
    t(
      GetAssayData(
        integrated.seurat,
        assay = "integrated",
        layer = "data"
      )
    )
  )
  integrated$meta <- integrated.seurat@meta.data
  # safety check
  stopifnot(
    identical(
      rownames(integrated$df),
      rownames(integrated$meta)
    )
  )
} else {
  message(
    "Only one batch detected: skipping batch integration."
  )
}

# save data ####
write.csv(unintegrated$df,
          paste(
            "data/processed/", file_name, "_data_unintegrated.csv", sep = "")
)
write.csv(unintegrated$meta,
          paste(
            "data/processed/", file_name, "_meta_unintegrated.csv", sep = "")
)
if (length(batches) > 1) {
write.csv(integrated$df,
          paste(
            "data/processed/", file_name, "_data_integrated.csv", sep = "")
)
write.csv(integrated$meta,
          paste(
            "data/processed/", file_name, "_meta_integrated.csv", sep = "")
)
}
rm(list = ls())
