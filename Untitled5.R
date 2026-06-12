# Title: batch integration (N=3) (mean per well) v1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 12-06-2026

# load packages ####
library(data.table)
library(Seurat)
# set variables ####
file_name_N1 <- "SF240627_mPaintDR2_N2"
file_name_N2 <- "SF240704_mPaintDR2_N3"
file_name_N3 <- "SF240711_mPaintDR2_N4"
file_name <- "mPaint_DR2_N1.2.3"
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
# integrate data (canonical correlation analysis/ CCA) - identifies conserved feature correlation patterns ####
