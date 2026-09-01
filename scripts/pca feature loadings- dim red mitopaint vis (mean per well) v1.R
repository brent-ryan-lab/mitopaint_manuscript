# Title: pca feature loadings- dim red mitopaint vis (mean per well) v1
# Step: 6.1.4
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 11-08-2026

# load packages ####
library(data.table)
library(tidyverse)
library(tidytext)
library(ggplot2)
library(ggpubr)
# set file variables ####
file_name <- "mPaintSpace2_N1_N2_N3"
integrate_state <- "integrated"
redu_state <- "redu"
annot_colors <- c(
    CellROX = "#DC267F",
    TMRM = "#FFB000",
    `mt-Keima pH7` = "#23CC86",
    `mt-Keima pH4` = "#FE6100",
    Other = "grey")
feature_patterns <- c(
  CellROX = "cellrox",
  TMRM = "tmrm",
  `mt-Keima pH7` = "mt-keima ph7",
  `mt-Keima pH4` = "mt-keima ph4")
dims_plot <- c("PC_1", "PC_2")
plot_width <- 5
plot_height <- 5
# set function to load data ####
load_data <- function(file_name,
                      integrate_state,
                      redu_state) {
  pca_loadings <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_pca_loadings.csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames
  rownames(pca_loadings) <- pca_loadings$V1
  pca_loadings$V1 <- NULL
  return(pca_loadings)
}
# run function to load data ####
data <- load_data(file_name, integrate_state, redu_state)
# create function to calc % loadings by feature_patterns ####
calc_pca_loadings <- function(data,
                              feature_patterns) {
  # helper to assign feature type from named patterns
  get_feature_type <- function(x, feature_patterns) {
    for (nm in names(feature_patterns)) {
      if (str_detect(x, regex(feature_patterns[[nm]], ignore_case = TRUE))) {
        return(nm)
      }
    }
    "Other"
  }
  # reshape to long format
  data_long <- data %>%
    as.data.frame() %>%
    rownames_to_column("feature") %>%
    pivot_longer(
      cols = everything()[-1],
      names_to = "PC",
      values_to = "loading"
    )
  # classify feature type + direction
  data_long <- data_long %>%
    mutate(
      feature_type = vapply(
        feature,
        get_feature_type,
        character(1),
        feature_patterns = feature_patterns
      ),
      direction = if_else(loading >= 0, "Positive", "Negative"),
      abs_loading = abs(loading)
    )
  # aggregate absolute loadings
  data_long_summary <- data_long %>%
    group_by(PC, direction, feature_type) %>%
    summarise(
      weight = sum(abs_loading),
      .groups = "drop"
    )
  # convert to proportional within each PC-direction
  data_long_summary <- data_long_summary %>%
    group_by(PC, direction) %>%
    mutate(
      prop = weight / sum(weight)
    ) %>%
    ungroup()
  # create PC + directions axis label
  data_long_summary <- data_long_summary %>%
    mutate(
      PC_dir = paste0(PC, " ", if_else(direction == "Positive", "+", "−"))
    )
  # per bar stacking order
  data_long_summary <- data_long_summary %>%
    mutate(
      feature_type_ord = reorder_within(
        feature_type,
        prop,
        PC_dir
      )
    )
  # add percent labels
  # only label >5%
  # label is rounded to nearest whole percent
  data_long_summary <- data_long_summary %>%
    mutate(
      label = ifelse(prop > 0.10,   
                     paste0(round(prop * 100), "%"),
                     NA)
    )
  data_long_summary
}
# run function to calc % loadings by feature_patterns ####
data_loadings <- calc_pca_loadings(data, feature_patterns)
# create function to bar plot % loadings by feature_patterns ####
plot_pca_loadings <- function(data_loadings,
                              annot_colors,
                              dims_plot) {
  # subset plotting_df
  # only plot PCs given in dims_plot
  plotting_df <- data_loadings[data_loadings$PC %in% dims_plot,]
  # plot
  ggplot(
    plotting_df,
    aes(
      x = PC_dir,
      y = prop * 100,
      group = feature_type_ord,   
      fill = feature_type,
      color = feature_type
    )
  ) +
    # geom_col, with reverse stacking order (largest on bottom)
    geom_col(
      width = 0.8,
      colour = "black",
      position = position_stack(reverse = FALSE)
    ) +
    # geom label with percent labels, no label if too small
    geom_label(
      aes(label = ifelse(prop > 0.10,
                         sprintf("%.2f%%", prop * 100),
                         NA)),
      position = position_stack(vjust = 0.5, reverse = FALSE),
      colour = "black",        
      fill = "white",          
      label.size = 0.3,        
      size = 5,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    # color given by annot_colors
    scale_fill_manual(values = annot_colors
    ) +
    scale_color_manual(values = annot_colors
    )  +
    scale_x_reordered() +   
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, 105)
    ) +
    guides(color = "none") +
    labs(
      x = "Principal Component (direction)",
      y = "Percentage total loading (%)",
      fill = "Feature type"
    ) +
    theme_pubr() +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(
      plot.title = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title = element_text(size = 12),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.box.just = "left",
      legend.title = element_text(hjust = 0.5, size = 12),
      legend.text = element_text(size = 12),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 0)
    )
}
# create function to bar plot % loadings by feature_patterns ####
plot <- plot_pca_loadings(data_loadings,
                          annot_colors,
                          dims_plot)
plot
# save plot ####
ggsave(filename = paste0(
  "outputs/figures/",
  file_name,
  "_",
  "pca_loadings_bar",
  ".pdf"
),
plot,
width = plot_width,
height = plot_height,
units = "in",
dpi = 300)
rm(list = ls())
