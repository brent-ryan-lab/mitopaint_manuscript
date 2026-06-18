# Title: pca scree- dim red mitopaint vis (mean per well) v1
# Step: 6.1.2
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 18-06-2026

# load packages ####
library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyverse)
# set file variables ####
file_name <- "mPaintDR2_N2_N3_N4"
redu_state <- "redu"
integrate_state <- "integrated"
# load file ####
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
x_lab <- "Number of PCs"
y_lab <- "Cumulative % Variance Explained"
size_axis <- 8
size_point <- 1
plot_width <- 2.5
plot_height <- 2.5
text_nudge <- 6
# pca scree plot ####
var$nPC <- seq(from = 1, to = 50, by = 1)
p <- ggplot(var, aes(y = Cumulative_Percent_Variance, x = nPC)) +   
  geom_point(color = "black", size = size_point) +
  geom_line(colour = "black", linewidth = 0.6) +
  annotate("line", x = 0:10, y = var[10,5], color = "grey", linetype = "dashed") +
  annotate("line", x = 10, y = 0:var[10,5], color = "grey", linetype = "dashed") +
  annotate("text", x = 10+(text_nudge)*1.2, y = var[10,5]-(text_nudge)*0.8, label = paste(
    round(var[10,5], digits = 2),
    "%", sep =""), size = size_axis/3) +
  annotate("line", x = 0:2, y = var[2,5], color = "grey", linetype = "dashed") +
  annotate("line", x = 2, y = 0:var[2,5], color = "grey", linetype = "dashed") +
  annotate("text", x = 2+(text_nudge)*1.2, y = var[2,5]-(text_nudge)*0.8, label = paste(
    round(var[2,5], digits = 2),
    "%", sep =""), size = size_axis/3) +
  annotate("line", x = 0:50, y = var[50,5], color = "grey", linetype = "dashed") +
  annotate("line", x = 50, y = 0:var[50,5], color = "grey", linetype = "dashed") +
  annotate("text", x = 50+(text_nudge)*1.2, y = var[50,5]-(text_nudge)*0.8, label = paste(
    round(var[50,5], digits = 2),
    "%", sep =""), size = size_axis/3) +
  theme_pubr() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0,10,20,30,40,50)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0,25,50,75,100)) +
  coord_cartesian(ylim = c(0, 100+text_nudge*2), xlim = c(0, 50++text_nudge*2)) +
  labs(x = x_lab, y = y_lab) +
  theme(axis.text=element_text(size = size_axis), 
        axis.title=element_text(size = size_axis)) +
  theme(legend.position="none")
# save plot ####
ggsave(paste("outputs/figures/pca/", file_name, "scree.pdf"),
       p,
       width = plot_width,
       height = plot_height,
       units = "in",
       dpi = 300)
rm(list = ls())
