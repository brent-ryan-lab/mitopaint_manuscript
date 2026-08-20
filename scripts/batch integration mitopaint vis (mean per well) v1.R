# Title: batch integration mitopaint vis (mean per well) v1
# Step: 4.1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 04-08-2026

# load packages ####
library(data.table)
library(colorspace)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(ggpubr)
# set variables ####
file_name <- "mPaintSpace2_N1_N2_N3"
redu_state <- "redu"
integrate_state <- c("integrated", "unintegrated")
pastel_cols <- lighten(c("#440154FF","#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF"), amount = 0.3)
n_neighbors <- 30
n_epochs <- 500
# create function to load data ####
load_data <- function(file_name, integrate_state) {
  # load integrated/unintegrated and redu/nonredu data as df
  df <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_data_", integrate_state, "_", redu_state, ".csv", sep = ""), 
      header = TRUE)
  )
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
  # return list of drift corrected data, raw data, metadata, and drift fits
  return(list(
    df = df,
    meta = meta
  ))
}
# initialise plots list
plots <- list()
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
# create function to make seurat object ####
seurat_dim_red <- function(feature_matrix,
                           metadata,
                           assay_name = "MP",
                           seed = 42,
                           n_neighbors,
                           n_epochs) {
  # put data frame and meta into a Seurat
  seurat_obj <- CreateSeuratObject(
    # transpose turns rows = features, cols = wells
    counts = t(as.matrix(feature_matrix)),
    meta.data = metadata,
    assay = assay_name
  )
  DefaultAssay(seurat_obj) <- assay_name
  # copy counts layer to data layer
  seurat_obj <- SetAssayData(
    seurat_obj,
    assay = assay_name,
    layer = "data",
    new.data = GetAssayData(
      seurat_obj,
      assay = assay_name,
      layer = "counts"
    )
  )
  # scale data before dim red
  seurat_obj <- ScaleData(seurat_obj)
  # run PCA dim red
  seurat_obj <- RunPCA(
    seurat_obj,
    # PCA dim red on all features 
    features = rownames(seurat_obj),
    seed.use = seed
  )
  # save PCA embeddings
  pca_embeddings <- as.data.frame(
    Embeddings(seurat_obj, "pca")
  )
  # calculate variance explained (for axis labels)
  pc_sd <- seurat_obj[["pca"]]@stdev
  pc_var <- pc_sd^2
  pca_var <- data.frame(
    PC = paste0("PC", seq_along(pc_sd)),
    SD = pc_sd,
    Variance = pc_var,
    Percent_Variance = 100 * pc_var / sum(pc_var),
    Cumulative_Percent_Variance = cumsum(100 * pc_var / sum(pc_var))
  )
  # run UMAP dim red
  seurat_obj <- RunUMAP(
    seurat_obj,
    dims = NULL, reduction = NULL,
    # UMAP dim red on all features 
    features = rownames(seurat_obj),
    n.epochs = n_epochs,
    n.neighbors = n_neighbors,
    seed.use = seed
  )
  # save UMAP embeddings
  umap_embeddings <- as.data.frame(
    Embeddings(seurat_obj, "umap")
  )
  # return
  return(list(
    seurat = seurat_obj,
    pca_embeddings = pca_embeddings,
    pca_var = pca_var,
    umap_embeddings = umap_embeddings
  ))
}
# run function to make seurat object ####
dim_red <- map(
  data,
  function(data_obj) {
    seurat_dim_red(
      feature_matrix = data_obj$df,
      metadata = data_obj$meta,
      assay_name = "MP",
      seed = 42,
      n_neighbors = n_neighbors,
      n_epochs = n_epochs
    )
  }
)
# create function to plot pca ####
plot_pca <- function(data,
                     grouping_var,
                     title_text
                     ) {
  # assumes DMSO, CCCP and ROT are in the dataset
  compound_levels <- c("DMSO", "CCCP", "ROT")
  # plot_df is combined data and metadata
  df <- data[["pca_embeddings"]]
  meta <- data[["seurat"]]@meta.data
  plot_df <- df %>%
    cbind(meta)
  # Set the order of Compound only when plotting by Compound
  if (grouping_var == "Compound") {
    plot_df$Compound <- factor(
      plot_df$Compound,
      levels = compound_levels
    )
  }
  # pull PC variance explained 
  pc1_var <- round(data$pca_var$Percent_Variance[1], 2)
  pc2_var <- round(data$pca_var$Percent_Variance[2], 2)
  x_lab <- paste0("PC_1 (", pc1_var, "%)")
  y_lab <- paste0("PC_2 (", pc2_var, "%)")
  # plot pca
  ggplot(plot_df, aes(x = PC_1, y = PC_2, color = .data[[grouping_var]])) +
    geom_point(shape = 16, size = 1.5) +
    # color using pastel_cols for grouping_var
    scale_color_manual(values = pastel_cols) +
    # tidy theme
    theme_pubr() +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 9,
                                face = "bold"),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      legend.position = "none",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 7),
      panel.grid = element_blank()
    ) +
    labs(title = title_text,
         x = x_lab,
         y = y_lab)
}
# run function to plot pca ####
plots$corr_pca_batch <- plot_pca(dim_red$integrated, 
                                 grouping_var = "Batch",
                                 title_text = "PCA by Batch\n(Seurat CCA Corrected)"
)
plots$corr_pca_batch
plots$uncorr_pca_batch <- plot_pca(dim_red$unintegrated, 
                                   grouping_var = "Batch",
                                   title_text = "PCA by Batch\n(Uncorrected)"
)
plots$uncorr_pca_batch
plots$corr_pca_drug <- plot_pca(dim_red$integrated, 
                               grouping_var = "Compound",
                               title_text = "PCA by Compound\n(Seurat CCA Corrected)"
)
plots$corr_pca_drug
plots$uncorr_pca_drug <- plot_pca(dim_red$unintegrated, 
                                  grouping_var = "Compound",
                                  title_text = "PCA by Compound\n(Uncorrected)"
)
plots$uncorr_pca_drug 
# create function to plot umap ####
plot_umap <- function(data,
                     grouping_var,
                     title_text
) {
  # assumes DMSO, CCCP and ROT are in the dataset
  compound_levels <- c("DMSO", "CCCP", "ROT")
  # plot_df is combined data and metadata
  df <- data[["umap_embeddings"]]
  meta <- data[["seurat"]]@meta.data
  plot_df <- df %>%
    cbind(meta) 
  # Set the order of Compound only when plotting by Compound
  if (grouping_var == "Compound") {
    plot_df$Compound <- factor(
      plot_df$Compound,
      levels = compound_levels
    )
  }
  # plot umap
  ggplot(plot_df, aes(x = umap_1, y = umap_2, color = .data[[grouping_var]])) +
    geom_point(shape = 16, size = 1.5) +
    # color using pastel_cols for grouping_var
    scale_color_manual(values = pastel_cols) +
    # tidy theme
    theme_pubr() +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 9,
                                face = "bold"),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      legend.position = "none",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 7),
      panel.grid = element_blank()
    ) +
    labs(title = title_text)
}
# run function to plot umap ####
plots$corr_umap_batch <- plot_umap(dim_red$integrated, 
                                  grouping_var = "Batch",
                                  title_text = "UMAP by Batch\n(Seurat CCA Corrected)"
)
plots$corr_umap_batch
plots$uncorr_umap_batch <- plot_umap(dim_red$unintegrated, 
                                     grouping_var = "Batch",
                                     title_text = "UMAP by Batch\n(Uncorrected)"
)
plots$uncorr_umap_batch
plots$corr_umap_drug <- plot_umap(dim_red$integrated, 
                                  grouping_var = "Compound",
                                  title_text = "UMAP by Compound\n(Seurat CCA Corrected)"
)
plots$corr_umap_drug
plots$uncorr_umap_drug <- plot_umap(dim_red$unintegrated, 
                                    grouping_var = "Compound",
                                    title_text = "UMAP by Compound\n(Uncorrected)"
)
plots$uncorr_umap_drug
# save plots ####
# loop save all plots in plots list
# create output folders
dir.create(
  "outputs/figures/pca",
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  "outputs/figures/umap",
  recursive = TRUE,
  showWarnings = FALSE
)
# save all plots
iwalk(
  plots,
  function(plot, plot_name) {
    # save plots in corresponding subfolder for pca or umap
    plot_folder <- case_when(
      str_detect(plot_name, regex("pca", ignore_case = TRUE)) ~ "pca",
      str_detect(plot_name, regex("umap", ignore_case = TRUE)) ~ "umap",
      TRUE ~ "other"
    )
    dir.create(
      paste0("outputs/figures/", plot_folder),
      recursive = TRUE,
      showWarnings = FALSE
    )
    ggsave(
      filename = paste0(
        "outputs/figures/",
        plot_folder,
        "/",
        file_name,
        "_",
        plot_name,
        ".pdf"
      ),
      plot = plot,
      width = 2.3,
      height = 2.5
    )
  }
)
rm(list = ls())
