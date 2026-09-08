# Title: batch integration mitopaint data (mean per well) v1
# Step: 4
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 18-06-2026

# load packages ####
library(data.table)
library(Seurat)
library(purrr)
# set variables ####
batches_info <- list(
  N1 = list(
    file_name = "SF250813_mPaintFDA_1",
    dmso_wells = c("2_2","4_2","6_2","8_2","10_2","12_2",
                   "14_2","16_2","18_2","20_2","22_2",
                   "3_9","5_9","7_9","9_9","11_9","13_9",
                   "15_9","17_9","19_9","21_9","23_9")
  ),
  N2 = list(
    file_name = "SF250813_mPaintFDA_2",
    dmso_wells = c("2_8","4_8","6_8","8_8","10_8","12_8",
                   "14_8","16_8","18_8","20_8","22_8",
                   "3_15","5_15","7_15","9_15","11_15",
                   "13_15","15_15","17_15","19_15","21_15","23_15")
  ),
  N3 = list(
    file_name = "SF250813_mPaintFDA_3",
    dmso_wells = c("2_2","4_2","6_2","8_2","10_2","12_2",
                   "14_2","16_2","18_2","20_2","22_2",
                   "3_9","5_9","7_9","9_9","11_9","13_9",
                   "15_9","17_9","19_9","21_9","23_9")
  ),
  N4 = list(
    file_name = "SF250813_mPaintFDA_4",
    dmso_wells = c("2_2","4_2","6_2","8_2","10_2","12_2",
                   "14_2","16_2","18_2","20_2","22_2",
                   "3_9","5_9","7_9","9_9","11_9","13_9",
                   "15_9","17_9","19_9","21_9","23_9")
  ),
  N5 = list(
    file_name = "SF250813_mPaintFDA_5",
    dmso_wells = c("2_2","4_2","6_2","8_2","10_2","12_2",
                   "14_2","16_2","18_2","20_2","22_2",
                   "3_9","5_9","7_9","9_9","11_9","13_9",
                   "15_9","17_9","19_9","21_9","23_9")
  ),
  N6 = list(
    file_name = "SF250813_mPaintFDA_6",
    dmso_wells = c("2_2","4_2","6_2","8_2","10_2","12_2",
                   "14_2","16_2","18_2","20_2","22_2",
                   "3_9","5_9","7_9","9_9","11_9","13_9",
                   "15_9","17_9","19_9","21_9","23_9")
  ),
  N7 = list(
    file_name = "SF250813_mPaintFDA_7",
    dmso_wells = c("2_2","4_2","6_2","8_2","10_2","12_2",
                   "14_2","16_2","18_2","20_2","22_2",
                   "3_9","5_9","7_9","9_9","11_9","13_9",
                   "15_9","17_9","19_9","21_9","23_9")
  ),
  N8 = list(
    file_name = "SF250813_mPaintFDA_8",
    dmso_wells = c("2_2","4_2","6_2","8_2","10_2","12_2",
                   "14_2","16_2","18_2","20_2","22_2",
                   "3_9","5_9","7_9","9_9","11_9","13_9",
                   "15_9","17_9","19_9","21_9","23_9")
  )
)
file_name <- "mPaintFDA_N1_N2_N3_N4_N5_N6_N7_N8"
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
batches <- map(
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
# run function to make seurat object in df.seurat container ####
df.seurat <- map(
  batches,
  make_seurat_object
)
# integrate data into integrated.seurat ####
# first, dimensionality is reduced using canonical correlation analysis (CCA)
# this turns conserved feature correlation patterns into vectors
# then, finds integration anchors using mutual nearest neighbors (MNN)
# this finds data pairs which overlap in same neighborhoods
# these anchors represent data points with "similar profiles" across the dataset
# and are used as reference points (anchor set) for integration
integrated.seurat <- IntegrateData(
  anchorset = FindIntegrationAnchors(
    object.list = df.seurat,
    # only integrate features shared across all N
    anchor.features = Reduce(
      intersect,
      map(batches, ~ colnames(.x$df))
    ),
    scale = FALSE
  ),
  # k.weight is number of neighbors to consider when weighting anchors
  # k.weight is dependent on approx no. of well replicates
  k.weight = k_weight   
)
# put integrated data in a single table ####
integrated <- list()
# all integrated batch data saved in integrated$df
integrated$df <- as.data.frame(
  t(
    GetAssayData(
      integrated.seurat,
      assay = "integrated",
      layer = "data"
    )
  )
)
# all batch metadata saved in integrated$meta
integrated$meta <- integrated.seurat@meta.data
# put unintegrated data in a single table ####
unintegrated <- list()
# all unintegrated batch data saved in unintegrated$df
unintegrated$df <- purrr::map(
  batches,
  "df"
) |>
  dplyr::bind_rows()

# all unintegrated metadata saved in unintegrated$meta
unintegrated$meta <- purrr::map(
  batches,
  "meta"
) |>
  dplyr::bind_rows()
# save data ####
write.csv(integrated$df,
          paste(
            "data/processed/", file_name, "_data_integrated.csv", sep = "")
)
write.csv(integrated$meta,
          paste(
            "data/processed/", file_name, "_meta_integrated.csv", sep = "")
)
write.csv(unintegrated$df,
          paste(
            "data/processed/", file_name, "_data_unintegrated.csv", sep = "")
)
write.csv(unintegrated$meta,
          paste(
            "data/processed/", file_name, "_meta_unintegrated.csv", sep = "")
)
rm(list = ls())
