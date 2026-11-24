# Title: plate drift dim red vis (mean per well) v1
# Step: 2.3
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 24-11-2026

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
file_name <- "mPaintDR2_N2_N3_N4"
redu_state <- "redu"
integrate_state <- "integrated"
avg_profile <- FALSE
plot_width <- 12
plot_height <- 3
grid_width <- 3
plot_cond <- c("DMSO_0", "CCCP_30", "ROT_10")
# create function to load data ####
load_data <- function(file_name,
                      integrate_state,
                      redu_state) {
  # load umap embeddings
  umap_embeddings <- as.data.frame(
    fread(paste0("data/processed/",file_name,"_",integrate_state,
                 "_",redu_state,"_umap_embeddings.csv"),
          header = TRUE))
  rownames(umap_embeddings) <- umap_embeddings$V1
  umap_embeddings$V1 <- NULL
  # load PCA embeddings
  pca_embeddings <- as.data.frame(
    fread(paste0("data/processed/",file_name,"_",integrate_state,
                 "_",redu_state,"_pca_embeddings.csv"),
          header = TRUE))
  rownames(pca_embeddings) <- pca_embeddings$V1
  pca_embeddings$V1 <- NULL
  # load PCA var
  pca_var <- as.data.frame(
    fread(paste0("data/processed/",file_name,"_",integrate_state,
                 "_",redu_state,"_pca_var.csv"),
          header = TRUE))
  pca_var$V1 <- NULL
  # load metadata
  meta <- as.data.frame(
    fread(paste0("data/processed/",file_name,"_",integrate_state,
                 "_",redu_state,"_dimred_meta.csv"),
          header = TRUE))
  rownames(meta) <- meta$V1
  meta$V1 <- NULL
  # keep rows aligned between meta and umap embeddings
  meta <- meta[rownames(umap_embeddings), , drop = FALSE]
  return(list(
    meta = meta,
    umap_embeddings = umap_embeddings,
    pca_embeddings = pca_embeddings,
    pca_var = pca_var))
}
# load data ####
data <- load_data(
  file_name = file_name,
  integrate_state = integrate_state,
  redu_state = redu_state
)
# create function to plot pca ####
# plot pca grid: split by condition, colour by order ####

# axis labels
x_lab <- paste0(
  "PC_1 (",
  round(
    data$pca_var$Percent_Variance[1],
    2
  ),
  "%)"
)

y_lab <- paste0(
  "PC_2 (",
  round(
    data$pca_var$Percent_Variance[2],
    2
  ),
  "%)"
)

# initialise plots list
plots <- list()

# combine PCA embeddings and metadata once
plot_df <- data$pca_embeddings |>
  cbind(data$meta) |>
  dplyr::mutate(
    Condition = as.character(Condition),
    Order = as.numeric(Order)
  )

# find conditions for plots
# one separate plot will be made for each condition
condition_levels <- unique(
  plot_df$Condition
)

# consistent Order range across ALL plots
order_limits <- range(
  plot_df$Order,
  na.rm = TRUE
)


# create function to plot one condition ####
plot_pca_condition <- function(
    plot_df,
    condition_name,
    order_limits
) {
  
  # foreground = condition of interest
  fg_df <- plot_df |>
    dplyr::filter(
      Condition == condition_name
    )
  
  # background = every other condition
  bg_df <- plot_df |>
    dplyr::filter(
      Condition != condition_name
    )
  
  # plot
  ggplot2::ggplot() +
    
    # all other conditions plotted first in grey
    ggplot2::geom_point(
      data = bg_df,
      ggplot2::aes(
        x = PC_1,
        y = PC_2
      ),
      colour = "grey70",
      size = 1
    ) +
    
    # condition of interest plotted on top
    # and coloured by imaging order
    ggplot2::geom_point(
      data = fg_df,
      ggplot2::aes(
        x = PC_1,
        y = PC_2,
        colour = Order
      ),
      size = 1.5
    ) +
    
    # same Order colour scale for every subplot
    ggplot2::scale_colour_gradientn(
      colours = rev(
        colorspace::lighten(
          viridis::viridis(100),
          amount = 0.3
        )
      ),
      limits = order_limits
    ) +
    
    # labels
    ggplot2::labs(
      title = condition_name,
      x = x_lab,
      y = y_lab,
      colour = "Order"
    ) +
    
    # compact colour bar
    ggplot2::guides(
      colour = ggplot2::guide_colorbar(
        barheight = grid::unit(
          1.5,
          "cm"
        ),
        barwidth = grid::unit(
          0.35,
          "cm"
        )
      )
    ) +
    
    # theme
    ggpubr::theme_pubr() +
    ggplot2::theme(
      aspect.ratio = 1,
      
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 9,
        face = "bold"
      ),
      
      axis.text = ggplot2::element_text(
        size = 7
      ),
      
      axis.title = ggplot2::element_text(
        size = 7
      ),
      
      legend.title = ggplot2::element_text(
        size = 8
      ),
      
      legend.text = ggplot2::element_text(
        size = 7
      ),
      
      legend.position = "right",
      
      panel.grid = ggplot2::element_blank()
    )
}
# create function to keep subplot dimensions consistent ####
add_fixed_legend_space <- function(
    p,
    legend_width = 0.5
) {
  
  # extract legend
  leg <- cowplot::get_legend(
    p +
      ggplot2::theme(
        legend.position = "right"
      )
  )
  
  # remove legend from main plot
  p_no_leg <- p +
    ggplot2::theme(
      legend.position = "none"
    )
  
  # put plot and legend beside each other
  cowplot::plot_grid(
    p_no_leg,
    leg,
    ncol = 2,
    rel_widths = c(
      1,
      legend_width
    )
  )
}
# make individual plots for each condition ####
plots <- purrr::map(
  condition_levels,
  ~ add_fixed_legend_space(
    plot_pca_condition(
      plot_df = plot_df,
      condition_name = .x,
      order_limits = order_limits
    )
  )
)

# name plots by condition
names(plots) <- condition_levels
# arrange plots into grid ####
if ("all" %in% plot_cond) {
  pca_grid <- cowplot::plot_grid(
    plotlist = plots,
    ncol = grid_width
  )
} else {
  pca_grid <- cowplot::plot_grid(
    plotlist = plots[plot_cond],
    ncol = grid_width
  )
}
pca_grid

# plot umap grid: split by condition, colour by order ####

# combine UMAP embeddings and metadata once
umap_plot_df <- data$umap_embeddings |>
  cbind(data$meta) |>
  dplyr::mutate(
    Condition = as.character(Condition),
    Order = as.numeric(Order)
  )

# find conditions for plots
condition_levels_umap <- unique(
  umap_plot_df$Condition
)

# use same global Order range across all UMAP plots
# this ensures colours mean the same thing in every subplot
umap_order_limits <- range(
  umap_plot_df$Order,
  na.rm = TRUE
)


# create function to plot one condition ####
plot_umap_condition <- function(
    plot_df,
    condition_name,
    order_limits
) {
  
  # foreground = condition of interest
  fg_df <- plot_df |>
    dplyr::filter(
      Condition == condition_name
    )
  
  # background = every other condition
  bg_df <- plot_df |>
    dplyr::filter(
      Condition != condition_name
    )
  
  # plot
  ggplot2::ggplot() +
    
    # all other conditions plotted first in grey
    ggplot2::geom_point(
      data = bg_df,
      ggplot2::aes(
        x = umap_1,
        y = umap_2
      ),
      colour = "grey70",
      size = 1
    ) +
    
    # condition of interest plotted on top
    # and coloured by imaging order
    ggplot2::geom_point(
      data = fg_df,
      ggplot2::aes(
        x = umap_1,
        y = umap_2,
        colour = Order
      ),
      size = 1.5
    ) +
    
    # same Order colour scale for every subplot
    ggplot2::scale_colour_gradientn(
      colours = rev(
        colorspace::lighten(
          viridis::viridis(100),
          amount = 0.3
        )
      ),
      limits = order_limits
    ) +
    
    # labels
    ggplot2::labs(
      title = condition_name,
      x = "UMAP_1",
      y = "UMAP_2",
      colour = "Order"
    ) +
    
    # compact colour bar
    ggplot2::guides(
      colour = ggplot2::guide_colorbar(
        barheight = grid::unit(
          1.5,
          "cm"
        ),
        barwidth = grid::unit(
          0.35,
          "cm"
        )
      )
    ) +
    
    # theme
    ggpubr::theme_pubr() +
    ggplot2::theme(
      aspect.ratio = 1,
      
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 9,
        face = "bold"
      ),
      
      axis.text = ggplot2::element_text(
        size = 7
      ),
      
      axis.title = ggplot2::element_text(
        size = 7
      ),
      
      legend.title = ggplot2::element_text(
        size = 8
      ),
      
      legend.text = ggplot2::element_text(
        size = 7
      ),
      
      legend.position = "right",
      
      panel.grid = ggplot2::element_blank()
    )
}


# make individual plots for each condition ####
umap_plots <- purrr::map(
  condition_levels_umap,
  ~ add_fixed_legend_space(
    plot_umap_condition(
      plot_df = umap_plot_df,
      condition_name = .x,
      order_limits = umap_order_limits
    )
  )
)

# name plots by condition
names(umap_plots) <- condition_levels_umap


# arrange plots into grid ####
if ("all" %in% plot_cond) {
  
  # plot every condition
  umap_grid <- cowplot::plot_grid(
    plotlist = umap_plots,
    ncol = grid_width
  )
  
} else {
  
  # plot selected conditions
  # order follows the order specified in plot_cond
  umap_grid <- cowplot::plot_grid(
    plotlist = umap_plots[plot_cond],
    ncol = grid_width
  )
}

umap_grid
# save plots ####
ggsave(
  filename = paste0("outputs/figures/umap/", file_name, "_umap_drift_grid.pdf"),
  plot = umap_grid,
  width = plot_width,
  height = plot_height
)
ggsave(
  filename = paste0("outputs/figures/pca/", file_name, "_pca_drift_grid.pdf"),
  plot = pca_grid,
  width = plot_width,
  height = plot_height
)
rm(list = ls())
