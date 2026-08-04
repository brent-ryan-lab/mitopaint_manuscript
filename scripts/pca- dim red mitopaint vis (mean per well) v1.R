# Title: pca- dim red mitopaint vis (mean per well) v1
# Step: 6.1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 18-06-2026

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
pastel_cols <- lighten(c("#440154FF", "#238A8DFF", "#FDE725FF"), amount = 0.3)
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
# create function to plot pca ####
# axis labels
x_lab <- paste0("PC_1 (",round(data$pca_var$Percent_Variance[1],2),"%)")
y_lab <- paste0("PC_2 (",round(data$pca_var$Percent_Variance[2],2),"%)")

plot_pca <- function(data,
                     colour_var,
                     title_text,
                     colour_scale,
                     shape_var = NULL,
                     legend_title = NULL) {
  # combine embeddings and meta in plot_df
  plot_df <- data$pca_embeddings |>
    cbind(data$meta)
  # reorder compounds when appropriate
  if (all(unique(plot_df$Compound) %in% c("DMSO", "CCCP", "ROT"))) {
    plot_df$Compound <- factor(
      plot_df$Compound,
      levels = c("DMSO", "CCCP", "ROT")
    )
  }
  # plot pca (aes)
  p <- ggplot(plot_df,aes(x = PC_1,y = PC_2,
                     colour = .data[[colour_var]]))
  # plot pca (shape)
  if (is.null(shape_var)) {
    p <- p +
      geom_point(size = 2)
  } else {
    p <- p +
      geom_point(aes(shape = .data[[shape_var]]),
        size = 2)
  }
  # plot pca (theme)
  p +
    colour_scale +
    labs(
      title = title_text,
      colour = legend_title,
      shape = shape_var,
      x = x_lab,
      y = y_lab
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
# p1: pc1 x pc2 by compound ####
# initialise plots list
plots <- list()
plots$pca_compound <- plot_pca(
  data,
  colour_var = "Compound",
  title_text = "PCA by Compound",
  legend_title = "Compound",
  colour_scale = scale_colour_manual(
    values = pastel_cols
  )
)
plots$pca_compound
# p2: pc1 x pc2 by concentration ####
plots$pca_concentration <- plot_pca(
  data,
  colour_var = "Concentration",
  title_text = "PCA by Concentration",
  legend_title = "Concentration (μM)",
  colour_scale = scale_colour_gradientn(
    colours = rev(lighten(viridis(100),
        amount = 0.3))
  )
)
plots$pca_concentration
# p3: pc1 x pc2 by pca nn ####
# assign floating labels positions for NN
plot_df <- data$pca_embeddings |>
  cbind(data$meta)
labels <- plot_df |>
  group_by(PCA_NN) |>
  summarise(
    PC_1 = median(PC_1),
    PC_2 = median(PC_2),
    .groups = "drop"
  )
plots$pca_nn <- ggplot(
  plot_df,aes(x = PC_1,y = PC_2,
              colour = as.factor(PCA_NN))
  ) +
  geom_point(size = 2) +
  geom_label_repel(
    data = labels,
    aes(x = PC_1,y = PC_2,label = PCA_NN,
      colour = as.factor(PCA_NN)),
    fill = "white",
    size = 3,
    linewidth = 0.35,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    segment.colour = NA,
    show.legend = FALSE
  ) +
  labs(title = "PCA by PCA NN",
      colour = "PCA_NN",
      x = x_lab,
      y = y_lab) +
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
plots$pca_nn
# p4: pc1 x pc1 by compound and concentration ####
plots$pca_compound_concentration <- plot_pca(
  data,
  colour_var = "Concentration",
  shape_var = "Compound",
  title_text = "PCA by Compound & Concentration",
  legend_title = "Concentration",
  colour_scale = scale_colour_gradientn(
    colours = rev(lighten(viridis(100),
                          amount = 0.3))
  )
)
plots$pca_compound_concentration
# p5: pc1 x pc2 by batch ####
plots$pca_batch <- plot_pca(
  data,
  colour_var = "Batch",
  title_text = "PCA by Batch",
  legend_title = "Batch",
  colour_scale = scale_colour_manual(
    values = pastel_cols
  )
)
plots$pca_batch
# create function to add fixed legend space ####
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
    rel_widths = c(
      plot_width,
      legend_width
    )
  )
}
# apply consistent legend space to plots ####
plots_fixed <- map(
  plots,
  add_fixed_legend_space
)
# save plots ####
iwalk(
  plots_fixed,
  function(plot, plot_name) {
    ggsave(
      filename = paste0(
        "outputs/figures/pca/",
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
