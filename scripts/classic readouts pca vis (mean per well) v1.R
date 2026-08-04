# Title: classic readouts pca vis (mean per well) v1
# Step: 9.1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 04-08-2026

# load packages ####
library(data.table)
library(colorspace)
library(ggplot2)
library(ggpubr)
library(viridis)
library(tidyverse)
library(cowplot)
# set variables ####
file_name_paint <- "mPaintDR2_N2_N3_N4"
file_name_classic <- "mPaint_DR2_Classic_N2_3_4"
redu_state <- "redu"
integrate_state <- "integrated"
# create function to load data ####
load_data <- function(file_name_paint, file_name_classic) {
  # load df
  df <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name_paint, "_data_", integrate_state, "_", redu_state, ".csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames as WELL_BATCH
  rownames(df) <- df$V1
  df$V1 <- NULL
  # load classic
  classic <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name_classic, "_norm", ".csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames as WELL_BATCH
  rownames(classic) <- classic$V1
  classic$V1 <- NULL
  # load pca embeddings
  pca_embeddings <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name_paint, "_", integrate_state, "_", redu_state, "_pca_embeddings.csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames as WELL_BATCH
  rownames(pca_embeddings) <- pca_embeddings$V1
  pca_embeddings$V1 <- NULL
  # load pca_var
  pca_var <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name_paint, "_", integrate_state, "_", redu_state, "_pca_var.csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames
  pca_var$V1 <- NULL
  # load metadata as meta
  meta <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name_paint, "_meta_", integrate_state, ".csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames as WELL_BATCH
  rownames(meta) <- meta$V1
  meta$V1 <- NULL
  # return list of drift corrected data, raw data, metadata, and drift fits
  return(list(
    df = df,
    classic = classic,
    pca_embeddings =  pca_embeddings,
    pca_var = pca_var,
    meta = meta
  ))
}
# initialise plots list
plots <- list()
# run function to load data ####
# data is a large list containing dataframes: df, classic, pca_embeddings, pca_var, meta
data <- load_data(
  file_name_paint,
  file_name_classic
)

# create function to plot pca ####
plot_pca <- function(data,
                     grouping_var,
                     legend_text,
                     title_text) {
  # plot_df is combined data and metadata
  df <- data[["pca_embeddings"]]
  classic <- data[["classic"]]
  plot_df <- df %>%
    cbind(classic)
  # pull PC variance explained 
  pc1_var <- round(data$pca_var$Percent_Variance[1], 2)
  pc2_var <- round(data$pca_var$Percent_Variance[2], 2)
  x_lab <- paste0("PC_1 (", pc1_var, "%)")
  y_lab <- paste0("PC_2 (", pc2_var, "%)")
  # plot pca
  ggplot(plot_df, aes(x = PC_1, y = PC_2, color = .data[[grouping_var]])) +
    # color viridis for grouping_var
    geom_point(shape = 16, size = 2) +
    scale_color_gradientn(
      colours = rev(lighten(
        viridis(100),
        amount = 0.3
      ))
    ) +
    # title and pca_var labels
    labs(title = title_text,
         color = legend_text,
         x = x_lab,
         y = y_lab) +
    # tidy theme
    theme_pubr() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5,
                                size = 9,
                                face = "bold"),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      legend.position = "none",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      panel.grid = element_blank()
    ) +
    theme(legend.position = "right")
}
# create function to add fixed legend space
add_fixed_legend_space <- function(plot,
                                   plot_width = 1,
                                   legend_width = 0.4) {
  legend <- get_legend(
    plot +
      theme(
        legend.position = "right"
      )
  )
  plot_without_legend <- plot +
    theme(
      legend.position = "none"
    )
  plot_grid(
    plot_without_legend,
    legend,
    nrow = 1,
    rel_widths = c(plot_width, legend_width)
  )
}
# run function to plot pca ####
plots$pca_mmp <- plot_pca(data, 
                          "Intensity Cytoplasm TMRM Mean",
                          "Cytoplasm\nMMP Intensity",
                          "PCA by MMP")
plots$pca_mmp
plots$pca_ros <- plot_pca(data, 
                          "Intensity Cytoplasm CellRox Mean",
                          "Cytoplasm\nROS Intensity",
                          "PCA by ROS")
plots$pca_ros
plots$pca_spots <- plot_pca(data, 
                            "Number of Selected Spots/ Selected Cell" ,
                            "Mitophagy Spots\nper Cell",
                            "PCA by Mitophagy Spots")
plots$pca_spots
plots$pca_morph <- plot_pca(data, 
                           "Mitochondria Selected Ratio Width to Length" ,
                           "Mitochondria Ratio\nWidth to Length",
                           "PCA by Mitochondria Morphology")
plots$pca_morph
plots_fixed <- map(
  plots,
  # apply fixed legend space to all plots so that PCA is square (not squished), and legend is consistent width
  add_fixed_legend_space
)
# save plots ####
iwalk(
  plots_fixed,
  function(plot, plot_name) {
    
    ggsave(
      filename = paste0(
        "outputs/figures/",
        file_name_paint,
        "_",
        plot_name,
        ".pdf"
      ),
      plot = plot,
      width = 3.7,
      height = 3
    )
  }
)
rm(list = ls())
