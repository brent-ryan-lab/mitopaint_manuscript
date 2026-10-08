# Title: nn bar- dim red mitopaint vis (mean per well) v1
# Step: 6.4
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 07-10-2026

# load packages ####
library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyverse)
# set variables ####
file_name <- "mPaintFDA_N1_N2_N3_N4_N5_N6_N7_N8"
redu_state <- "redu"
integrate_state <- "integrated"
avg_profile <- TRUE
pos_ctrls <- c("SodiumArsenite", "ROT", "Oligomycin", "Nocodazole", "MLN4924",
                "MitoQ", "CytochalasinD", "Chloroquine", "CCCP", "AntimycinA")

# load meta ####
if (avg_profile == TRUE) {
  meta <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_avg_dimred_meta.csv", sep = ""),
      header = TRUE
    )
  )
} else {
  meta <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_dimred_meta.csv", sep = ""),
      header = TRUE
    )
  )
}
rownames(meta) <- meta$V1
meta$V1 <- NULL

# add Compound_ID column ####
meta <- meta |>
  dplyr::mutate(
    Compound_ID = dplyr::case_when(
      Compound %in% pos_ctrls ~ Compound,
      Compound == "DMSO" ~ "DMSO",
      TRUE ~ "Library"
    )
  )
# set function to plot NN bar ####
plot_bar <- function(meta,
                     group_by,
                     color_by) {
  meta <- meta
  plot_meta <- meta[, colnames(meta) %in% c(group_by, color_by)]
  plot_meta <- plot_meta |>
    dplyr::mutate(
      Compound_ID = factor(
        Compound_ID,
        levels = c("DMSO", "Library", pos_ctrls)
      )
    )
  
  p <- ggplot(plot_meta, aes(x = factor(.data[[group_by]]), fill = .data[[color_by]])) +
    geom_bar(position = "fill") +
    scale_fill_manual(
      values = c(
        DMSO = "grey40",
        Library = "grey60",
        setNames(
          scales::hue_pal()(length(pos_ctrls)),
          pos_ctrls
        )
      ),
      drop = FALSE
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
                   color_by = "Compound_ID")
pca_nn
umap_nn <- plot_bar(meta = meta,
                   group_by = "UMAP_NN",
                   color_by = "Compound_ID")
umap_nn
pca_nn_by_umap_nn <- ggplot(meta, aes(x = factor(.data[["PCA_NN"]]), fill = factor(.data[["UMAP_NN"]]))) +
  geom_bar(position = "fill") +
  scale_fill_manual(
    values = c(
      DMSO = "grey40",
      Library = "grey60",
      setNames(
        scales::hue_pal()(length(unique(meta$UMAP_NN))),
        unique(meta$UMAP_NN)
      )
    ),
    drop = FALSE
  ) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    y = "Percentage",
    x = "PCA_NN",
    fill = "UMAP_NN"
  ) +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))
pca_nn_by_umap_nn 
# save plots ####
ggsave(
  paste0("outputs/figures/pca/", file_name, "_pca_nn_bar.pdf"),
  plot = pca_nn,
  width = 6, height = 3, units = "in", dpi = 300
)
ggsave(
  paste0("outputs/figures/umap/", file_name, "_umap_nn_bar.pdf"),
  plot = umap_nn,
  width = 9, height = 3, units = "in", dpi = 300
)
