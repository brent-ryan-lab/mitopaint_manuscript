# Title: mahalanobis distance mitopaint data (mean per well) v1
# Step: 7.
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 07-08-26

# load packages ####
library(data.table)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(viridis)
# set file variables ####
file_name <- "mPaintDR2_N2_N3_N4"
integrate_state <- "integrated"
redu_state <- "redu"
pc_use <- 10
# create function to load data ####
load_data <- function(file_name,
                      integrate_state,
                      redu_state) {
  # load pca_embeddings
  pca_embeddings <- as.data.frame(
    fread(
      paste("data/processed/", file_name, "_", integrate_state, 
            "_" ,redu_state, "_pca_embeddings", ".csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames as WELL_BATCH
  rownames(pca_embeddings) <- pca_embeddings$V1
  pca_embeddings$V1 <- NULL
  # load metadata
  meta <- as.data.frame(
    fread(
      paste("data/processed/", file_name, "_", "meta_", 
            integrate_state, ".csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames as WELL_BATCH
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
    pca_var = pca_var
  ))
}
# load data ####
data <- load_data(file_name, integrate_state, redu_state)
# create function to calculate mahalanobis ####
mahalanobis_data <- function(obj,
                             pc_use) {
  pca_embeddings <- obj$pca_embeddings
  meta <- obj$meta
  # subset dmso
  dmso_id <- meta$Compound == "DMSO" 
  pca_use <- pca_embeddings[, 1:pc_use, drop = FALSE]
  dmso_pca <- pca_use[dmso_id, , drop = FALSE]
  # dmso mean vector in PC space
  dmso_mu <- colMeans(dmso_pca)
  # dmso covariance matrix in PC space
  dmso_sigma <- cov(dmso_pca)
  # calculate mahalanobis
  mahalanobis_dist <- as.data.frame(mahalanobis(
    x = pca_use,
    center = dmso_mu,
    cov = dmso_sigma
  ))
  colnames(mahalanobis_dist) <- paste("mahalanobis_distance_", "PC_", pc_use, sep = "")
  return(mahalanobis_dist)
}
# run function to calculate mahalanobis ####
md <- mahalanobis_data(data, pc_use)
# create function to pca plot by mahalanobis ####
# axis labels
x_lab <- paste0("PC_1 (",round(data$pca_var$Percent_Variance[1],2),"%)")
y_lab <- paste0("PC_2 (",round(data$pca_var$Percent_Variance[2],2),"%)")
plots <- list()
plot_pca <- function(data,
                     md,
                     colour_var,
                     title_text,
                     colour_scale,
                     shape_var = NULL,
                     legend_title = NULL) {
  # combine embeddings, meta, md in plot_df
  plot_df <- data$pca_embeddings[,1:2] |>
    cbind(data$meta) |>
    cbind(md)
  # plot pca (aes)
  p <- ggplot(plot_df,aes(x = PC_1,y = PC_2,
                          colour = .data[[colour_var]]))
  # plot pca (shape)
    p <- p +
      geom_point(size = 2)
  # plot pca (theme)
  p <- p +
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
  return(p)
}
# run function to create pca plot by mahalanobis ####
plots$pca_md<- plot_pca(
  data,
  md,
  colour_var = colnames(md)[1],
  title_text = "PCA by Mahalanobis Distance to DMSO",
  legend_title = "Mahalanobis\nDistance",
  colour_scale = scale_colour_gradientn(
    colours = rev(lighten(viridis(100),
                          amount = 0.3))
  )
)
plots$pca_md
# save data ####
write.csv(
  md,
  paste("data/processed/", file_name, "_mahal_", "PC", pc_use,".csv", sep = "")
)
# save plots ####
ggsave(
  filename = paste0(
    "outputs/figures/pca/",
    file_name,
    "_",
    "pca_md",
    ".pdf"
  ),
  plot = plots$pca_md,
  width = 3.7,
  height = 3,
  units = "in",
  dpi = 300
)
rm(list = ls())