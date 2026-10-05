# Title: batch integration mitopaint vis (mean per well) v2
# Step: 4.1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 01-10-2026

# load packages ####
library(data.table)
library(colorspace)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(ggpubr)
# set variables ####
file_name <- "mPaintFDA_N1_N2_N3_N4_N5_N6_N7_N8"
redu_state <- "redu"
integrate_state <- c("integrated", "unintegrated")
pastel_cols <- lighten(c("#440154FF","#46337EFF","#365C8DFF","#277F8EFF","#1FA187FF","#4AC16DFF","#9FDA3AFF","#FDE725FF"), amount = 0.3)
n_neighbors <- 30
n_epochs <- 500
avg_profile <- TRUE
# create function to load data ####
load_data <- function(file_name, integrate_state) {
  if (avg_profile == TRUE) {
    # load integrated/unintegrated and redu/nonredu pca as pca
    pca <- as.data.frame(
      fread(
        paste(
          "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_pca_embeddings.csv", sep = ""), 
        header = TRUE)
    )
    # keep rownames as WELL_BATCH
    rownames(pca) <- pca$V1
    pca$V1 <- NULL
    # load pca variance
    pca_var <- as.data.frame(
      fread(
        paste(
          "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_pca_var.csv", sep = ""), 
        header = TRUE)
    )
    # keep rownames as WELL_BATCH
    rownames(pca_var) <- pca_var$V1
    pca_var$V1 <- NULL
    # load integrated/unintegrated and redu/nonredu umap as umap
    umap <- as.data.frame(
      fread(
        paste(
          "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_avg_umap_embeddings.csv", sep = ""), 
        header = TRUE)
    )
    rownames(umap) <- umap$V1
    umap$V1 <- NULL
    # load integrated/unintegrated and redu/nonredu meta as meta
    meta <- as.data.frame(
      fread(
        paste(
          "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_avg_dimred_meta.csv", sep = ""), 
        header = TRUE)
    )
  }
  else {
    # load integrated/unintegrated and redu/nonredu pca as pca
    pca <- as.data.frame(
      fread(
        paste(
          "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_pca_embeddings.csv", sep = ""), 
        header = TRUE)
    )
    # keep rownames as WELL_BATCH
    rownames(pca) <- pca$V1
    pca$V1 <- NULL
    # load pca variance
    pca_var <- as.data.frame(
      fread(
        paste(
          "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_pca_var.csv", sep = ""), 
        header = TRUE)
    )
    # keep rownames as WELL_BATCH
    rownames(pca_var) <- pca_var$V1
    pca_var$V1 <- NULL
    # load integrated/unintegrated and redu/nonredu umap as umap
    umap <- as.data.frame(
      fread(
        paste(
          "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_umap_embeddings.csv", sep = ""), 
        header = TRUE)
    )
    rownames(umap) <- umap$V1
    umap$V1 <- NULL
    # load integrated/unintegrated and redu/nonredu meta as meta
    meta <- as.data.frame(
      fread(
        paste(
          "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_dimred_meta.csv", sep = ""), 
        header = TRUE)
    )
    }
  # return list of drift corrected data, raw data, metadata, and drift fits
  return(list(
    pca = pca,
    pca_var = pca_var,
    umap = umap,
    meta = meta
  ))
}
# initialise plots list
plots <- list()
# run function to load data ####
# data is a large list containing sublists for integrate_state (integrated, unintegrated)
data <- integrate_state |>
  set_names() |>
  purrr::map(
    # each sublist contains corresponding df and meta
    ~ load_data(
      file_name = file_name,
      integrate_state = .x
    )
  )

# create function to plot pca ####
plot_pca <- function(data,
                     grouping_var,
                     title_text) {
  # assumes DMSO, CCCP and ROT are in the dataset
  compound_levels <- c("DMSO", "CCCP", "ROT")
  # plot_df is combined data and metadata
  df <- data[["pca"]]
  meta <- data[["meta"]]
  plot_df <- df %>%
    cbind(meta)
  # pull PC variance explained 
  pc1_var <- round(data$pca_var$Percent_Variance[1], 2)
  pc2_var <- round(data$pca_var$Percent_Variance[2], 2)
  x_lab <- paste0("PC_1 (", pc1_var, "%)")
  y_lab <- paste0("PC_2 (", pc2_var, "%)")
  # if grouping_var is compound:
  # set order to compound_levels
  if (grouping_var == "Compound") {
    plot_df$Compound_plot <- ifelse(
      plot_df$Compound %in% compound_levels,
      as.character(plot_df$Compound),
      "Other"
    )
    # set colours of thos not in compound level to grey70
    plot_df$Compound_plot <- factor(
      plot_df$Compound_plot,
      levels = c(compound_levels, "Other")
    )
    cols <- c(
      setNames(lighten(viridis::viridis(3), amount = 0.3), compound_levels),
      Other = "grey70"
    )
    # plot main structure
    p <- ggplot(plot_df, aes(x = PC_1, y = PC_2, color = Compound_plot)) +
      geom_point(shape = 16, size = 1.5) +
      # plot colors (if group_var is compound)
      scale_color_manual(values = cols)
  } else {
    # plot main structure
    p <- ggplot(plot_df, aes(x = PC_1, y = PC_2, color = .data[[grouping_var]])) +
      geom_point(shape = 16, size = 1.5) +
      # plot colors (if group_var is anything else)
      scale_color_manual(values = pastel_cols)
  }
  # tidy theme
  p +
    theme_pubr() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 9,
        face = "bold"
      ),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      legend.position = "none",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 7),
      panel.grid = element_blank()
    ) +
    labs(
      title = title_text,
      x = x_lab,
      y = y_lab
    )
}
# run function to plot pca ####
plots$corr_pca_batch <- plot_pca(data$integrated, 
                                 grouping_var = "Batch",
                                 title_text = "PCA by Batch\n(Seurat CCA Corrected)"
)
plots$corr_pca_batch
plots$uncorr_pca_batch <- plot_pca(data$unintegrated, 
                                   grouping_var = "Batch",
                                   title_text = "PCA by Batch\n(Uncorrected)"
)
plots$uncorr_pca_batch
plots$corr_pca_drug <- plot_pca(data$integrated, 
                               grouping_var = "Compound",
                               title_text = "PCA by Compound\n(Seurat CCA Corrected)"
)
plots$corr_pca_drug
plots$uncorr_pca_drug <- plot_pca(data$unintegrated, 
                                  grouping_var = "Compound",
                                  title_text = "PCA by Compound\n(Uncorrected)"
)
plots$uncorr_pca_drug 
# create function to plot umap ####
plot_umap <- function(data,
                      grouping_var,
                      title_text) {
  compound_levels <- c("DMSO", "CCCP", "ROT")
  df <- data[["umap"]]
  meta <- data[["meta"]]
  plot_df <- df %>%
    cbind(meta)
  if (grouping_var == "Compound") {
    plot_df$Compound_plot <- ifelse(
      plot_df$Compound %in% compound_levels,
      as.character(plot_df$Compound),
      "Other"
    )
    plot_df$Compound_plot <- factor(
      plot_df$Compound_plot,
      levels = c(compound_levels, "Other")
    )
    cols <- c(
      setNames(lighten(viridis::viridis(3), amount = 0.3), compound_levels),
      Other = "grey70"
    )
    p <- ggplot(plot_df, aes(x = umap_1, y = umap_2, colour = Compound_plot)) +
      geom_point(shape = 16, size = 1.5) +
      scale_color_manual(values = cols)
  } else {
    p <- ggplot(plot_df, aes(x = umap_1, y = umap_2, colour = .data[[grouping_var]])) +
      geom_point(shape = 16, size = 1.5) +
      scale_color_manual(values = pastel_cols)
  }
  p +
    theme_pubr() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 9,
        face = "bold"
      ),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      legend.position = "none",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 7),
      panel.grid = element_blank()
    ) +
    labs(
      title = title_text,
      x = "UMAP_1",
      y = "UMAP_2"
    )
}
# run function to plot umap ####
plots$corr_umap_batch <- plot_umap(data$integrated, 
                                  grouping_var = "Batch",
                                  title_text = "UMAP by Batch\n(Seurat CCA Corrected)"
)
plots$corr_umap_batch
plots$uncorr_umap_batch <- plot_umap(data$unintegrated, 
                                     grouping_var = "Batch",
                                     title_text = "UMAP by Batch\n(Uncorrected)"
)
plots$uncorr_umap_batch
plots$corr_umap_drug <- plot_umap(data$integrated, 
                                  grouping_var = "Compound",
                                  title_text = "UMAP by Compound\n(Seurat CCA Corrected)"
)
plots$corr_umap_drug
plots$uncorr_umap_drug <- plot_umap(data$unintegrated, 
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

legends <- purrr::map(
  plots,
  ~ cowplot::get_legend(.x + theme(legend.position = "right"))
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

iwalk(legends, function(leg, plot_name) {
  plot_folder <- case_when(
    str_detect(plot_name, regex("pca", ignore_case = TRUE)) ~ "pca",
    str_detect(plot_name, regex("umap", ignore_case = TRUE)) ~ "umap",
    TRUE ~ "other"
  )
  dir.create(
    paste0("outputs/figures/", plot_folder, "/legends"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  ggsave(
    filename = paste0(
      "outputs/figures/",
      plot_folder,
      "/legends/",
      file_name,
      "_",
      plot_name,
      "_legend.pdf"
    ),
    plot = cowplot::plot_grid(leg),
    width = 2.5,
    height = 3.5
  )
})

rm(list = ls())
