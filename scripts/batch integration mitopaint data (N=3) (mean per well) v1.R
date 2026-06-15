# Title: batch integration mitopaint data (N=3) (mean per well) v1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 15-06-2026

# load packages ####
library(data.table)
library(Seurat)
# set variables ####
file_name_N1 <- NULL
file_name_N2 <- NULL
file_name_N3 <- NULL
file_name <- NULL
k_weight <- NULL
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
N1 <- load_data(file_name_N1)
N2 <- load_data(file_name_N2)
N3 <- load_data(file_name_N3)
# make df.seurat container of each N as a seurat object ####
df.seurat <- list()
# create seurat object of N1, N2, N3 zscore and meta in for() loop
for(batch in c("N1", "N2", "N3")) {
  # fetch zscore and metadata
  z_mat   <- get(batch)$df
  meta_df <- get(batch)$meta
  # transpose zscore matrix
  mat <- t(as.matrix(z_mat))
  # create seurat object for each batch
  df.seurat[[batch]] <- CreateSeuratObject(
    counts = mat,
    # assay is saved as MP for mitopaint
    assay = "MP",
    meta.data = meta_df,
    # include all wells (skip usual RNAseq filter)
    min.cells = 0,
    min.features = 0
  )
  # manually copy counts to data layer to circumvent v5 seurat structure syntax
  df.seurat[[batch]] <- SetAssayData(
    df.seurat[[batch]],
    assay = "MP",
    layer = "data",
    new.data = GetAssayData(df.seurat[[batch]], assay = "MP", layer = "counts")
  )
  # scale data - centre normalises, and is necessary for seurat package functions
  df.seurat[[batch]] <- ScaleData(
    df.seurat[[batch]],
    verbose = FALSE
  )
}
# remove temp variables from for() loop
rm(z_mat, meta_df, mat)
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
      list(colnames(N1$df),
           colnames(N2$df),
           colnames(N3$df))
    ),
    scale = FALSE
  ),
  # k.weight is number of neighbors to consider when weighting anchors
  # k.weight is dependent on approx no. of well replicates
  k.weight = k_weight   
)

# put integrated data in a single table ####
# all INTEGRATED N data saved in integrated$df
integrated <- list()
integrated$df <- as.data.frame(
  t(GetAssayData(integrated.seurat, assay = "integrated", slot = "data"))
)
# all N meta saved in integrated$meta
integrated$meta <- integrated.seurat@meta.data
# put unintegrated data in a single table ####
# all UNINTEGRATED N data saved in integrated$df
unintegrated <- list()
unintegrated$df <- rbind(
  N1$df,
  N2$df,
  N3$df
)
unintegrated$meta <- rbind(
  N1$meta,
  N2$meta,
  N3$meta
)
# save data ####
# save integrated data in /data/processed
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
