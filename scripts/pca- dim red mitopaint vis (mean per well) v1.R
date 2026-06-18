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
# set file variables ####
file_name <- "mPaintDR2_N2_N3_N4"
redu_state <- "redu"
integrate_state <- "integrated"
# load file ####
# load data as df
df <- as.data.frame(
  fread(
    paste(
      "data/processed/", file_name, "_", integrate_state ,"_", redu_state, "_pca_embeddings.csv", sep = ""), 
    header = TRUE)
)
# keep rownames as WELL_BATCH
rownames(df) <- df$V1
df$V1 <- NULL
# load metadata as meta
meta <- as.data.frame(
  fread(
    paste(
      "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_dimred_meta.csv", sep = ""), 
    header = TRUE)
)
# keep rownames as WELL_BATCH
rownames(meta) <- meta$V1
meta$V1 <- NULL
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
# set plot variables ####
# plotting variables
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
size_title <- 10
size_axis <- 8
size_point <- 2
plot_width <- 4
plot_height <- 3.2
pastel_cols <- lighten(c("#440154FF", "#238A8DFF", "#FDE725FF"), amount = 0.3)
# p1: pc1 x pc2 by compound ####
# make plots container to hold plots
plots <- list()
# make data cbind of meta and df
data <- cbind(meta, df)
# reorder plot data if only treated with DMSO, CCCP, ROT
if (all(unique(meta$Compound) %in% c("DMSO", "CCCP", "ROT"))) {
  data <- data %>%
    mutate(
      Compound = factor(
        Compound,
        levels = c("DMSO", "CCCP", "ROT")
      )
    )
} else {
  data <- data
}
# plot is p1
plots$p1 <- ggplot(data,aes(x = PC_1, y = PC_2)
                   ) +
  geom_point(aes(color = Compound), size = size_point) +
  labs(title = "PCA by Compound",
       x = x_lab,
       y = y_lab) +
  theme_pubr() +
  theme(
    plot.title = element_text(hjust = 0.5, size = size_title, face = "bold"),
    axis.text = element_text(size = size_axis),
    axis.title = element_text(size = size_axis),
    legend.title = element_text(size = size_axis),
    legend.text = element_text(size = size_axis),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_color_manual(
    values = pastel_cols
    ) + 
  theme(legend.position = "right")
plots$p1 
# p2: pc1 x pc2 by concentration ####
plots$p2 <- ggplot(data, aes(x = PC_1, y = PC_2)) +
  geom_point(aes(color = Concentration), size = size_point) +
  labs(title = "PCA by Concentration",
       x = x_lab,
       y = y_lab) +
  theme_pubr() +
  theme(
    plot.title = element_text(hjust = 0.5, size = size_title, face = "bold"),
    axis.text = element_text(size = size_axis),
    axis.title = element_text(size = size_axis),
    legend.title = element_text(size = size_axis),
    legend.text = element_text(size = size_axis),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_color_gradientn(
    colours = rev(lighten(
      viridis(100),
      amount = 0.3
    ))
  ) +
  theme(legend.position = "right")
plots$p2
# p3: pc1 x pc2 by pca nn ####
labels <- data %>%
  group_by(PCA_NN) %>%
  summarise(
    PC_1 = median(PC_1),
    PC_2 = median(PC_2),
    .groups = "drop"
  )
plots$p3 <- ggplot(data, aes(x = PC_1, y = PC_2)) +
  geom_point(aes(color = as.factor(PCA_NN)), size = size_point) +
  labs(title = "PCA by PCA NN",
       x = x_lab,
       y = y_lab,
       color = "PCA_NN") +
  theme_pubr() +
  theme(
    plot.title = element_text(hjust = 0.5, size = size_title, face = "bold"),
    axis.text = element_text(size = size_axis),
    axis.title = element_text(size = size_axis),
    legend.title = element_text(size = size_axis),
    legend.text = element_text(size = size_axis),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  geom_label_repel(
    data = labels,
    aes(
      x = PC_1,
      y = PC_2,
      label = PCA_NN,
      color = as.factor(PCA_NN)   # match label color to group
    ),
    fill = "white",
    size = 3,
    label.size = 0.35,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    segment.color = NA,
    show.legend = FALSE
  ) +
  theme(legend.position = "right")
plots$p3
# p4: pc1 x pc1 by compound and concentration ####
plots$p4 <- ggplot(data,aes(x = PC_1, y = PC_2)) +
  geom_point(aes(color = Concentration, shape = Compound), size = size_point) +
  labs(title = "PCA by Compound & Concentration",
       x = x_lab,
       y = y_lab) +
  theme_pubr() +
  theme(
    plot.title = element_text(hjust = 0.5, size = size_title, face = "bold"),
    axis.text = element_text(size = size_axis),
    axis.title = element_text(size = size_axis),
    legend.title = element_text(size = size_axis),
    legend.text = element_text(size = size_axis),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_color_gradientn(
    colours = rev(lighten(
      viridis(100),
      amount = 0.3
    ))
  ) +
  theme(legend.position = "right")
plots$p4
# p5: pc1 x pc2 by batch ####
plots$p5 <- ggplot(data, aes(x = PC_1, y = PC_2)
) +
  geom_point(aes(color = Batch), size = size_point) +
  labs(title = "PCA by Batch",
       x = x_lab,
       y = y_lab) +
  theme_pubr() +
  theme(
    plot.title = element_text(hjust = 0.5, size = size_title, face = "bold"),
    axis.text = element_text(size = size_axis),
    axis.title = element_text(size = size_axis),
    legend.title = element_text(size = size_axis),
    legend.text = element_text(size = size_axis),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_color_manual(
    values = pastel_cols
  ) + 
  theme(legend.position = "right")
plots$p5
# save plots ####
ggsave(paste("outputs/figures/pca/", file_name, "pca_1.pdf"),
       plots$p1,
       width = plot_width,
       height = plot_height,
       units = "in",
       dpi = 300)
ggsave(paste("outputs/figures/pca/", file_name, "pca_2.pdf"),
       plots$p2,
       width = plot_width,
       height = plot_height,
       units = "in",
       dpi = 300)
ggsave(paste("outputs/figures/pca/", file_name, "pca_3.pdf"),
       plots$p3,
       width = plot_width,
       height = plot_height,
       units = "in",
       dpi = 300)
ggsave(paste("outputs/figures/pca/", file_name, "pca_4.pdf"),
       plots$p4,
       width = plot_width,
       height = plot_height,
       units = "in",
       dpi = 300)
ggsave(paste("outputs/figures/pca/", file_name, "pca_5.pdf"),
       plots$p5,
       width = plot_width,
       height = plot_height,
       units = "in",
       dpi = 300)
rm(list = ls())
