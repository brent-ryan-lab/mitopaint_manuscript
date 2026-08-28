# Title: pca biplot- dim red mitopaint vis (mean per well) v1
# Step: 6.1.3
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 19-06-2026

# load packages ####
library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyverse)
# set file variables ####
file_name <- "mPaintSpace2_N1_N2_N3"
redu_state <- "redu"
integrate_state <- "integrated"
# load data ####
# load pca variance as var
var <- as.data.frame(
  fread(
    paste(
      "data/processed/", file_name, "_", integrate_state ,"_", redu_state, "_pca_var.csv", sep = ""), 
    header = TRUE)
)
# keep rownames as WELL_BATCH
rownames(var) <- var$V1
var$V1 <- NULL
# load loadings as loadings
loadings <- as.data.frame(
  fread(
    paste(
      "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_pca_loadings.csv", sep = ""), 
    header = TRUE)
)
rownames(loadings) <- loadings$V1
loadings$V1 <- NULL
# load pca embeddings as embeddings
embeddings <- as.data.frame(
  fread(
    paste(
      "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_pca_embeddings.csv", sep = ""), 
    header = TRUE)
)
# keep rownames as WELL_BATCH
rownames(embeddings) <- embeddings$V1
embeddings$V1 <- NULL
# set plot variables ####
# chan_cols color order is TMRM, CellROX, mt-Keima pH7, mt-Keima pH4
chan_cols <- c("#FFB000", "#DC267F", "#23CC86", "#FE6100")
x_lab <- paste0(
  "PC_1 (",
  round(var[1, "Percent_Variance"], 2),
  "%)"
)
y_lab <- paste0(
  "PC_2 (",
  round(var[2, "Percent_Variance"], 2),
  "%)"
)
point_size <- 2
arrow_width <- 1.2
arrowhead_size <- 0.5
label_text_size <- 2.5
label_line_width <- 0.4
size_title <- 13
size_axis <- 10
plot_width <- 7
plot_height <- 7
# create vector of top +/- 5 features for PC1 and PC2 ####
# first collate all feature names for biplot
plot_feats <- unique(c(
  rownames(slice_max(loadings, n = 5, order_by = PC_1)),
  rownames(slice_min(loadings, n = 5, order_by = PC_1)),
  rownames(slice_max(loadings, n = 5, order_by = PC_2)),
  rownames(slice_min(loadings, n = 5, order_by = PC_2))
))
# then wrap feature names
# add \n every 25 characters for legible labels in plot
plot_feats_tidy <- str_wrap(plot_feats, width = 25)
# create df for arrows ####
# collate loadings for each feature in each PC
plot_df <- data.frame(
  feature = plot_feats,
  feature_tidy = plot_feats_tidy,
  PC_1 = loadings[plot_feats, "PC_1"],
  PC_2 = loadings[plot_feats, "PC_2"],
  row.names = NULL
)
# scale arrows for pca space ####
# this calculates the euclidean distance from the origin 
# ie. this tells how strong the feature loads into the PC1-2 plane overall
plot_df$loading_norm <- sqrt(plot_df$PC_1^2 + 
                               plot_df$PC_2^2) 
# the arrow scaling factor is the same value for all arrows
# this is relative to the range of the PCA axis
# therefore the largest loading in the table extends to 60% of the plot radius
# this keeps all arrows size not too big or small to the plot
arrow_scale <- 0.6 * max(
  abs(embeddings$PC_1), 
  abs(embeddings$PC_2)
) / max(plot_df$loading_norm) 
# multiply PC loadings by arrow scale
biplot_df <- plot_df %>% mutate(
  PC_1 = PC_1 * arrow_scale, 
  PC_2 = PC_2 * arrow_scale)
# colour arrows by feature channel ####
biplot_df <- biplot_df %>%
  mutate(
    feature_color = case_when(
      grepl("TMRM test", feature, ignore.case = TRUE) ~ chan_cols[1],
      grepl("CellRox", feature, ignore.case = TRUE) ~ chan_cols[2],
      grepl("mt-Keima pH7", feature, ignore.case = TRUE) ~ chan_cols[3],
      grepl("mt-Keima pH4", feature, ignore.case = TRUE) ~ chan_cols[4],
      TRUE ~ "black"))
# plot biplot ####
plot <- ggplot(embeddings, aes(x = PC_1, y = PC_2)) +
  ## grey out points, point size = 2
  geom_point(color = "grey70", size = point_size) +
  ## add loading arrows
  geom_segment(
    data = biplot_df,
    aes(x = 0, y = 0,
        xend = PC_1, yend = PC_2,
        color = feature_color),
    arrow = arrow(length = unit(arrowhead_size, "cm")),
    linewidth = arrow_width,
    show.legend = FALSE) +
  ## invisible points are labels anchors set to repel from arrow tips
  geom_point(
    data = biplot_df,
    aes(x = PC_1, y = PC_2),
    # alpha = 0, makes anchors invisible
    alpha = 0,          
    size = 10) +
  ## add labels to loading arrows
  geom_label_repel(
    data = biplot_df,
    aes(x = PC_1, y = PC_2,
        # feature label is tidied name with /n every 25 characters
        label = feature_tidy,  
        # feature label is coloured by channel
        color = feature_color   
    ),
    # label text size
    size = label_text_size,
    # label line width
    label.size = label_line_width,
    fill = "white",
    # space between labels
    box.padding = 0.4,        
    # space away from arrow tips
    point.padding = 0.6,      
    # stronger repulsion
    force = 5,               
    # weak pull back to anchor
    force_pull = 0.2,         
    max.overlaps = Inf,
    direction = "both",
    show.legend = FALSE,
    # do not link label and arrow with pointer
    segment.color = NA        
  ) +
  scale_color_identity() +
  ## axis labels
  labs(title = "PCA Biplot (Top 5 ± Features for PC1 & PC2)",
       x = x_lab,
       y = y_lab) +
  # theme
  theme_pubr() +
  # axis + title sizes
  theme(
    axis.text = element_text(size = size_axis),
    axis.title = element_text(size = size_axis),
    plot.title = element_text(hjust = 0.5, size = size_title, face = "bold"),
    panel.grid = element_blank())
plot
# save biplot ####
ggsave(paste("outputs/figures/pca/", file_name, "pca_biplot.pdf"),
       plot,
       width = plot_width,
       height = plot_height,
       units = "in",
       dpi = 300)
rm(list = ls())
