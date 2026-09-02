# Title: classic readouts pca pearson corr heatmap vis (mean per well) v1
# Step: 9.1.2
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 07-08-26

# load packages ####
library(data.table)
library(stringr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(purrr)
library(colorspace)
library(ComplexHeatmap)
library(circlize)
# set file variables ####
classic_file_name <- "mPaintSpace2_N1_N2_N3_Classic"
file_name <- "mPaintSpace2_N1_N2_N3"
integrate_state <- "integrated"
redu_state <- "redu"
meta_cols <- c("Row", "Column", "Compound", "Concentration", "Well", "Batch", "Condition")
feature_patterns <- c(
  MMP = "tmrm",
  ROS = "cellrox",
  Morph = "mitochondria",
  Spots = "spots|mt-keima|ph4/ph7|spot")
annot_colors <- list(
  Feature = c(
    ROS = "#DC267F",
    MMP = "#FFB000",
    Morph = "#23CC86",
    Spots = "#FE6100",
    Other = "grey"))
col_scale <- c("blue", "white", "red")
dims_plot <- c("PC_1", "PC_2")
# set function to load data ####
load_data <- function(classic_file_name,
                      file_name,
                      integrate_state,
                      redu_state,
                      meta_cols) {
  # load pca_embeddings
  pca_embeddings <- as.data.frame(
    fread(
      paste("data/processed/", file_name, "_", integrate_state, 
            "_" ,redu_state, "_pca_embeddings", ".csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames as WELL_BATCH
  rownames(pca_embeddings) <- pca_embeddings$V1
  pca_embeddings$V1 <- NULL
  # load classic readouts
  df <- as.data.frame(
    fread(
      paste("data/processed/", classic_file_name, ".csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames as WELL_BATCH
  rownames(df) <- df$V1
  df$V1 <- NULL
  df <- df[rownames(df) %in% rownames(pca_embeddings),]
  # separate metadata
  meta <- df[, colnames(df) %in% meta_cols]
  df <- df[, !colnames(df) %in% meta_cols]
  return(list(pca_embeddings = pca_embeddings,
              meta = meta,
              df = df
              ))
}
# run function to load data ####
data <- load_data(classic_file_name,
                  file_name,
                  integrate_state,
                  redu_state,
                  meta_cols)
# calculate pearson corr ####
cor_mat <- cor(data$df,
               data$pca_embeddings, 
               method = "pearson", use = "pairwise.complete.obs")
# create annotation to group feats ####
annotate_feature <- function(x, patterns) {
  for (label in names(patterns)) {
    if (str_detect(x, regex(patterns[[label]], ignore_case = TRUE))) {
      return(label)
    }
  }
  "Other"
}
row_annot <- data.frame(
  Feature = vapply(
    rownames(cor_mat),
    annotate_feature,
    character(1),
    patterns = feature_patterns
  ),
  row.names = rownames(cor_mat)
)
# create function to make and save heatmap ####
cor_heatmap <- function (cor_mat,
                         row_annot,
                         filename) {
  # create row annotation object
  row_annotation <- rowAnnotation(
    Feature = row_annot$Feature,
    col = annot_colors
  )
  # define heatmap colscale
  col_fun <- colorRamp2(c(-1, 0, 1), col_scale)
  # subset data
  cor_mat_subset <- cor_mat[,colnames(cor_mat) %in% dims_plot]
  # plot heatmap
  p <- draw(
    Heatmap(
      cor_mat_subset,
      name = "Pearson r",
      col = col_fun,
      # clustering
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      # annotations
      left_annotation = row_annotation,
      # aesthetics
      show_row_names = FALSE,
      show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = 10)
    ))
  if (!is.null(filename)) {
    pdf(filename, width = 3, height = 4.5)
    ComplexHeatmap::draw(p)
    dev.off()
  }
  return(p)
  }
# run function to make and save heatmap ####
plot <- cor_heatmap(cor_mat,
                    row_annot,
                    paste("outputs/figures/", classic_file_name, "_pca_cor_heatmap", ".pdf", sep = ""))
rm(list = ls())
