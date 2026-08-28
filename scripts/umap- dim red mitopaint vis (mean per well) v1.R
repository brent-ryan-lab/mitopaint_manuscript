# Title: umap- dim red mitopaint vis (mean per well) v1
# Step: 6.2
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 05-08-2026

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
pastel_cols <- scales::hue_pal()(33)
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
    umap_embeddings = umap_embeddings,
    meta = meta))
}
# load data ####
data <- load_data(
  file_name = file_name,
  integrate_state = integrate_state,
  redu_state = redu_state
)
# create function to plot umap ####
plot_umap <- function(data,
                     colour_var,
                     title_text,
                     colour_scale,
                     shape_var = NULL,
                     legend_title = NULL) {
  # combine embeddings and meta in plot_df
  plot_df <- data$umap_embeddings |>
    cbind(data$meta)
  # reorder compounds when appropriate
  if (all(unique(plot_df$Compound) %in% c("DMSO", "CCCP", "ROT"))) {
    plot_df$Compound <- factor(
      plot_df$Compound,
      levels = c("DMSO", "CCCP", "ROT")
    )
  }
  # plot umap (aes)
  p <- ggplot(plot_df,aes(x = umap_1,y = umap_2,
                          colour = .data[[colour_var]]))
  # plot umap (shape)
  if (is.null(shape_var)) {
    p <- p +
      geom_point(size = 2)
  } else {
    p <- p +
      geom_point(aes(shape = .data[[shape_var]]),
                 size = 2)
  }
  # plot umap (theme)
  p +
    colour_scale +
    labs(
      title = title_text,
      colour = legend_title,
      shape = shape_var
    ) +
    # tidy theme
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
# p1: umap1 x umap2 by compound ####
# initialise plots list
plots <- list()
plots$umap_compound <- plot_umap(
  data,
  colour_var = "Compound",
  title_text = "UMAP by Compound",
  legend_title = "Compound",
  colour_scale = scale_colour_manual(
    values = pastel_cols
  )
)
plots$umap_compound
# p2: umap1 x umap2 by concentration ####
plots$umap_concentration <- plot_umap(
  data,
  colour_var = "Concentration",
  title_text = "UMAP by Concentration",
  legend_title = "Concentration (uM)",
  colour_scale = scale_colour_gradientn(
    colours = rev(lighten(viridis(100),
                          amount = 0.3))
  )
)
plots$umap_concentration
# p3: umap1 x umap2 by umap nn ####
# assign floating labels positions for NN
plot_df <- data$umap_embeddings |>
  cbind(data$meta)
labels <- plot_df |>
  group_by(UMAP_NN) |>
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  )
plots$umap_nn <- ggplot(
  plot_df,aes(x = umap_1,y = umap_2,
              colour = as.factor(UMAP_NN))
) +
  geom_point(size = 2) +
  geom_label_repel(
    data = labels,
    aes(x = umap_1,y = umap_2,label = UMAP_NN,
        colour = as.factor(UMAP_NN)),
    fill = "white",
    size = 3,
    linewidth = 0.35,
    fontface = "bold",
    box.padding = 0.15,
    point.padding = 0.2,
    max.overlaps = Inf,
    segment.colour = NA,
    show.legend = FALSE
  ) +
  labs(title = "UMAP by UMAP NN",
       colour = "UMAP_NN") +
  # color legend across two cols
  guides(
    colour = guide_legend(
      ncol = 2,
      byrow = TRUE
    )
  ) +
  # tidy theme
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
plots$umap_nn
# p4: umap1 x umap1 by compound and concentration ####
# plots$umap_compound_concentration <- plot_umap(
#   data,
#   colour_var = "Concentration",
#   shape_var = "Compound",
#   title_text = "UMAP by Compound & Concentration",
#   legend_title = "Concentration (uM)",
#   colour_scale = scale_colour_gradientn(
#     colours = rev(lighten(viridis(100),
#                           amount = 0.3))
#   )
# )
# plots$umap_compound_concentration
# create function to add fixed legend space ####
add_fixed_legend_space <- function(plot,
                                   plot_name,
                                   plot_width = 1,
                                   legend_width = 0.4,
                                   legend_cutoff = 11,
                                   legend_file = NULL) {
  # get number of legend entries from the plot data
  n_legend <- NULL
  # detect compound plot specifically
  if (plot_name == "umap_compound") {
    n_legend <- length(unique(data$meta$Compound))
  }
  # extract legend
  legend <- cowplot::get_legend(
    plot +
      theme(legend.position = "right")
  )
  plot_without_legend <- plot +
    theme(legend.position = "none")
  # if the legend is too long, save it separately and leave blank space in the plot
  if (!is.null(n_legend) && n_legend > legend_cutoff) {
    if (!is.null(legend_file)) {
      ggsave(
        filename = legend_file,
        plot = cowplot::plot_grid(legend),
        width = 3,
        height = 4,
        units = "in",
        dpi = 300
      )
    }
    return(
      plot_grid(
        plot_without_legend,
        ggplot() + theme_void(),
        nrow = 1,
        rel_widths = c(plot_width, legend_width)
      )
    )
  }
  # otherwise keep the legend beside the plot
  plot_grid(
    plot_without_legend,
    legend,
    nrow = 1,
    rel_widths = c(plot_width, legend_width)
  )
}
# apply consistent legend space to plots ####
plots_fixed <- map(
  names(plots),
  function(plot_name) {
    add_fixed_legend_space(
      plot = plots[[plot_name]],
      plot_name = plot_name,
      legend_cutoff = 11,
      legend_file = if (plot_name == "umap_compound") {
        paste0("outputs/figures/umap/", file_name, "_umap_compound_legend.pdf")
      } else {
        NULL
      }
    )
  }
)
names(plots_fixed) <- names(plots)
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
        file_name,
        "_",
        plot_name,
        ".pdf"
      ),
      plot = plot,
      width = 3.7,
      height = 3,
      units = "in",
      dpi = 300
    )
  }
)
rm(list = ls())
