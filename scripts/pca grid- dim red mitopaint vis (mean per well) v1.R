# Title: pca grid- dim red mitopaint vis (mean per well) v1
# Step: 6.1.5
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 25-08-2026

# load packages ####
library(data.table)
library(ggplot2)
library(colorspace)
library(ggpubr)
library(tidyverse)
library(viridis)
library(ggrepel)
library(cowplot)
# set variables ####
file_name <- "mPaintSpace2_N1_N2_N3"
redu_state <- "redu"
integrate_state <- "integrated"
plot_width <- 12
plot_height <- 10
grid_width <- 6
# create function to load data ####
load_data <- function(file_name,
                      integrate_state,
                      redu_state) {
  # load PCA embeddings
  pca_embeddings <- as.data.frame(
    fread(paste0("data/processed/",file_name,"_",integrate_state,
                 "_",redu_state,"_pca_embeddings.csv"),
          header = TRUE))
  rownames(pca_embeddings) <- pca_embeddings$V1
  pca_embeddings$V1 <- NULL
  # load metadata
  meta <- as.data.frame(
    fread(paste0("data/processed/",file_name,"_",integrate_state,
                 "_",redu_state,"_dimred_meta.csv"),
          header = TRUE))
  rownames(meta) <- meta$V1
  meta$V1 <- NULL
  # load PCA var
  pca_var <- as.data.frame(
    fread(paste0("data/processed/",file_name,"_",integrate_state,
                 "_",redu_state,"_pca_var.csv"),
          header = TRUE))
  pca_var$V1 <- NULL
  # keep rows aligned between meta and PCA embeddings
  meta <- meta[rownames(pca_embeddings), , drop = FALSE]
  return(list(
    pca_embeddings = pca_embeddings,
    meta = meta,
    pca_var = pca_var))
}
# load data ####
data <- load_data(
  file_name = file_name,
  integrate_state = integrate_state,
  redu_state = redu_state
)
# plot pca_grid: split by compound, colour by concentration ####
# axis labels
x_lab <- paste0("PC_1 (",round(data$pca_var$Percent_Variance[1],2),"%)")
y_lab <- paste0("PC_2 (",round(data$pca_var$Percent_Variance[2],2),"%)")
# initialise plots list
plots <- list()
# find levels for plots (separate plot for each compound)
compound_levels <- unique(data$meta$Compound)
# make individual plots by compound
plot_pca_compound <- function(data, compound_name) {
  # pull data for plotting
  plot_df <- data$pca_embeddings |>
    cbind(data$meta) |>
    mutate(
      Compound = as.character(Compound),
      Concentration = as.numeric(Concentration)
    )
  # foreground = the compound of interest
  fg_df <- plot_df |>
    filter(Compound == compound_name)
  # background = all other compounds
  bg_df <- plot_df |>
    filter(Compound != compound_name)
  # plot
  ggplot() +
    # background points first
    geom_point(
      data = bg_df,
      aes(x = PC_1, y = PC_2),
      colour = "grey70",
      size = 1
    ) +
    # compound of interest on top, coloured by concentration
    geom_point(
      data = fg_df,
      aes(x = PC_1, y = PC_2, colour = Concentration),
      size = 1.5
    ) +
    scale_colour_gradientn(
      colours = rev(
        lighten(
          viridis(100),
          amount = 0.3
        )
      )
    ) +
    labs(
      title = compound_name,
      x = x_lab,
      y = y_lab,
      colour = "[uM]"
    ) +
    guides(
      colour = guide_colorbar(
        barheight = grid::unit(1.5, "cm"),
        barwidth = grid::unit(0.35, "cm")
      )
    ) +
    theme_pubr() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(
        hjust = 0.5,
        size = 9,
        face = "bold"
      ),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.position = "right",
      panel.grid = element_blank()
    )
}
# make all compound plots 
# fix legend spacing so plots all end up same dims on grid
add_fixed_legend_space <- function(p, legend_width = 0.5) {
  leg <- cowplot::get_legend(
    p + theme(legend.position = "right")
  )
  p_no_leg <- p + theme(legend.position = "none")
  cowplot::plot_grid(
    p_no_leg,
    leg,
    ncol = 2,
    rel_widths = c(1, legend_width)
  )
}
# make all individual plots for grid in loop
plots <- purrr::map(
  compound_levels,
  ~ add_fixed_legend_space(plot_pca_compound(data, .x))
)
# arrange into a grid #
pca_grid <- cowplot::plot_grid(
  plotlist = plots,
  ncol = grid_width
)
pca_grid

# save plot ####
ggsave(
  filename = paste0("outputs/figures/", file_name, "_pca_grid.pdf"),
  plot = pca_grid,
  width = plot_width,
  height = plot_height
)
rm(list = ls())