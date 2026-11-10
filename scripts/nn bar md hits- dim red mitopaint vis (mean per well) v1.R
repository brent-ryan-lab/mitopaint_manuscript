# Title: nn bar md hits- dim red mitopaint vis (mean per well) v1
# Step: 6.4.2
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 10-11-2026

# load packages ####
library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyverse)
# set variables ####
file_name <- "mPaintFDA_N1_N2_N3_N4_N5_N6_N7_N8"
redu_state <- "redu"
integrate_state <- "integrated"
# load meta ####
meta <- as.data.frame(
  fread(
    paste(
      "data/processed/", file_name, "_meta_summary.csv", sep = ""),
    header = TRUE
  )
)
rownames(meta) <- meta$V1
meta$V1 <- NULL
# set function to plot NN bar ####
plot_bar <- function(meta,
                     group_by,
                     color_by) {
  meta <- meta
  plot_meta <- meta[, colnames(meta) %in% c(group_by, color_by)]
  
  p <- ggplot(plot_meta, aes(x = factor(.data[[group_by]]), fill = .data[[color_by]])) +
    geom_bar(position = "fill") +
    scale_fill_manual(
      values = c(
        "FALSE" = "grey60",
        "TRUE" = "#00BFC4"
      )
    ) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(
      y = "Percentage",
      x = group_by,
      fill = color_by
    ) +
    theme_pubr() +
    theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8))
  p
}
# run function to plot NN bar ####
pca_nn <- plot_bar(meta = meta,
                   group_by = "PCA_NN",
                   color_by = "chi_sig")
pca_nn
umap_nn <- plot_bar(meta = meta,
                    group_by = "UMAP_NN",
                    color_by = "chi_sig")
umap_nn
# save plots ####
ggsave(
  paste0("outputs/figures/pca/", file_name, "_pca_nn_bar_hits.pdf"),
  plot = pca_nn,
  width = 6, height = 2.5, units = "in", dpi = 300
)
ggsave(
  paste0("outputs/figures/umap/", file_name, "_umap_nn_bar_hits.pdf"),
  plot = umap_nn,
  width = 9, height = 2.5, units = "in", dpi = 300
)
