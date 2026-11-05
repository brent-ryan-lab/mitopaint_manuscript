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
    file_name = "SF261028_mPaint2DR2_N1",
    dmso_wells = c(
      "2_9","2_2","3_2","3_9","4_9","5_2","5_9","6_9","6_2","7_2","7_9","8_9",
      "8_2","9_2","9_9","10_9","10_2","11_2","11_9","12_9","12_2","13_2","13_9",
      "14_9","14_2","15_2","15_9","16_9","16_2","17_2","17_9","18_9","18_2",
      "19_2","19_9","20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9"
    )
  )
)
file_name <- "mPaint2DR2_N1"
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
