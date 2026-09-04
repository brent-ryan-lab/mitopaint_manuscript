# Title: classic readouts umap vis (mean per well) v1
# Step: 9.3
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 04-09-2026

# load packages ####
library(data.table)
library(colorspace)
library(ggplot2)
library(ggpubr)
library(viridis)
library(tidyverse)
library(cowplot)
# set variables ####
file_name_paint <- "mPaintSpace2_N1_N2_N3"
file_name_classic <- "mPaintSpace2_N1_N2_N3_Classic"
redu_state <- "redu"
integrate_state <- "integrated"
plot_feats <- c("Intensity Cytoplasm TMRM test Mean",
                "Intensity Cytoplasm CellRox Deep Red test Mean",
                "Number of Mitophagy Spots Selected- per Cell",
                "mkeima ph7 mitochondria Ratio Width to Length")
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
  classic <- classic[rownames(classic) %in% rownames(df),]
  classic$V1 <- NULL
  # load umap embeddings
  umap_embeddings <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name_paint, "_", integrate_state, "_", redu_state, "_umap_embeddings.csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames as WELL_BATCH
  rownames(umap_embeddings) <- umap_embeddings$V1
  umap_embeddings$V1 <- NULL
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
    umap_embeddings =  umap_embeddings,
    meta = meta
  ))
}
# initialise plots list
plots <- list()
# run function to load data ####
# data is a large list containing dataframes: df, classic, umap_embeddings, meta
data <- load_data(
  file_name_paint,
  file_name_classic
)
# create function to plot umap ####
plot_umap <- function(data,
                     grouping_var,
                     legend_text,
                     title_text) {
  # plot_df is combined data and metadata
  df <- data[["umap_embeddings"]]
  classic <- data[["classic"]]
  plot_df <- df %>%
    cbind(classic)
  # plot umap
  ggplot(plot_df, aes(x = umap_1, y = umap_2, color = .data[[grouping_var]])) +
    # color viridis for grouping_var
    geom_point(shape = 16, size = 2) +
    scale_color_gradientn(
      colours = rev(lighten(
        viridis(100),
        amount = 0.3
      ))
    ) +
    # title
    labs(title = title_text,
         color = legend_text) +
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
# run function to plot umap ####
plots$umap_mmp <- plot_umap(data, 
                          plot_feats[1],
                          "Cytoplasm\nMMP Intensity",
                          "umap by MMP")
plots$umap_mmp
plots$umap_ros <- plot_umap(data, 
                          plot_feats[2],
                          "Cytoplasm\nROS Intensity",
                          "umap by ROS")
plots$umap_ros
plots$umap_spots <- plot_umap(data, 
                            plot_feats[3],
                            "Mitophagy Spots\nper Cell",
                            "umap by Mitophagy Spots")
plots$umap_spots
plots$umap_morph <- plot_umap(data, 
                            plot_feats[4],
                            "Mitochondria Ratio\nWidth to Length",
                            "umap by Mitochondria Morphology")
plots$umap_morph
plots_fixed <- map(
  plots,
  # apply fixed legend space to all plots so that umap is square (not squished), and legend is consistent width
  add_fixed_legend_space
)
# save plots ####
# create output folders
dir.create(
  "outputs/figures/umap",
  recursive = TRUE,
  showWarnings = FALSE
)
# save all plots
iwalk(
  plots_fixed,
  function(plot, plot_name) {
    
    ggsave(
      filename = paste0(
        "outputs/figures/umap/",
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