# Title: mahalanobis distance mitopaint data (mean per well) v1
# Step: 7.
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 07-08-26

# load packages ####
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(viridis)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(purrr)
library(cowplot)
# set file variables ####
file_name <- "mPaintDR2_N2_N3_N4"
integrate_state <- "integrated"
redu_state <- "redu"
pc_use <- 10
pastel_cols <- lighten(c("#440154FF", "#238A8DFF"), amount = 0.3)
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
  colnames(mahalanobis_dist) <- paste("mahalanobis_d2_", "PC_", pc_use, sep = "")
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
  title_text = "PCA by Mahalanobis D-squared to DMSO",
  legend_title = "Mahalanobis\nD-squared",
  colour_scale = scale_colour_gradientn(
    colours = rev(lighten(viridis(100),
                          amount = 0.3))
  )
)
plots$pca_md
# create function to calculate mahalanobis hits ####
mahalanobis_hits <- function(obj, md, pc_use) {
  # test 1 = mahalanobis d-squared p val to chi squared val
  # test 1 is calculated on a well level
  # test 2 = linear mixed effect model to find sig of mahalanobis d-squared, accounting for batch effect
  # test 2 is calculated on a condition level across batches
  # combine Mahalanobis distance with metadata
  md_df <- cbind(data$meta, md)
  # get the Mahalanobis distance column name
  md_col <- colnames(md)[1]
  # ensure grouping variables are factors
  md_df$Batch <- factor(md_df$Batch)
  md_df$Condition <- factor(md_df$Condition)
  # test 1, p-value for each well based on chi-square distribution
  md_df$p_md <- pchisq(
    md_df[[md_col]],
    df = pc_use,
    lower.tail = FALSE
  )
  # test 1, significance flag
  md_df$md_sig <- md_df$p_md < 0.05
  # test 1, optional chi-square cutoff for plotting / reporting
  md_cutoff <- qchisq(0.95, df = pc_use)
  # test 1, wells above the cutoff
  md_df$md_sig_cutoff <- md_df[[md_col]] > md_cutoff
  # test 1, mahalanobis chi squared summary 
  md_well_results <- md_df |>
    select(
      Batch,
      Condition,
      all_of(md_col),
      p_md,
      md_sig,
      md_sig_cutoff
    )
  # test 2, batch aware linear mixed model on mahalanobis distance
  md_model <- lmer(
    as.formula(
      paste0("`", md_col, "` ~ Condition + (1|Batch)")
    ),
    data = md_df
  )
  # test 2, estimated marginal means
  md_emm <- emmeans(
    md_model,
    ~ Condition
  )
  # test 2, contrasts only DMSO v. other conditions
  md_contrasts <- contrast(
    md_emm,
    method = "trt.vs.ctrl",
    ref = "DMSO_0",
    adjust = "bonferroni"
  )
  # test 2, model-based pairwise comparison table
  md_contrasts_df <- as.data.frame(md_contrasts)
  # test 2, mahalanobis lmme summary
  md_model_results <- md_contrasts_df
  # test 2, populate table with matching metadata
  # make a lookup table from the metadata
  condition_lookup <- md_df |>
    distinct(Condition, Compound, Concentration)
  # add Condition / Compound / Concentration to the contrast table
  md_model_results <- md_contrasts_df |>
    separate(
      contrast,
      into = c("Condition", "Control"),
      sep = " - ",
      remove = FALSE
    ) |>
    left_join(
      condition_lookup,
      by = "Condition"
    ) |>
    select(
      Condition,
      Compound,
      Concentration,
      contrast,
      estimate,
      SE,
      df,
      t.ratio,
      p.value,
      everything()
    )
  # return list
  return(list(
    md = md,
    md_chi_p = md_well_results,
    md_chi_cutoff = md_cutoff,
    md_lmme_p = md_model_results
  ))
}
# run function to calculate mahalanobis hits ####
md_stats <- mahalanobis_hits(data, md, pc_use)
# pca plot by mahalanobis p ####
# colour_var by p value
plots$pca_md_p <- plot_pca(
  data,
  md_stats$md_chi_p[,-c(1:2)],
  colour_var = "p_md",
  title_text = "PCA by Mahalanobis Chi-squared p-val",
  legend_title = "Mahalanobis\nChi-squared p-val",
  colour_scale = scale_colour_gradientn(
    colours = rev(lighten(viridis(100),
                          amount = 0.3))
  )
)
plots$pca_md_p
# colour_var by significance
plots$pca_md_hit <- plot_pca(
  data,
  md_stats$md_chi_p[,-c(1:2)],
  colour_var = "md_sig",
  title_text = "PCA by Mahalanobis Chi-squared Hit",
  legend_title = "Mahalanobis\nChi-squared Hit",
  colour_scale = scale_colour_manual(values = pastel_cols)
)
plots$pca_md_hit
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
    rel_widths = c(plot_width, legend_width)
  )
}
plots_fixed <- map(
  plots,
  # apply fixed legend space to all plots so that PCA is square (not squished), and legend is consistent width
  add_fixed_legend_space
)
plots_fixed
# save data ####
write.csv(
  md,
  paste("data/processed/", file_name, "_mahal_", "PC", pc_use,".csv", sep = "")
)
write.csv(
  md_stats$md_chi_p,
  paste("data/processed/", file_name, "_mahal_chi_p_", "PC", pc_use,".csv", sep = "")
)
write.csv(
  md_stats$md_lmme_p,
  paste("data/processed/", file_name, "_mahal_lmme_p_", "PC", pc_use,".csv", sep = "")
)
# wip: save plots ####
dir.create(
  "outputs/figures/pca",
  recursive = TRUE,
  showWarnings = FALSE
)
# save all plots
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
      height = 3
    )
  }
)
rm(list = ls())