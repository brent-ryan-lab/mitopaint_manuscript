# Title: similarity heatmap mitopaint vis - NN (mean per well) v1
# Step: 8.2
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 07-10-2026

# load packages ####
library(data.table)
library(ComplexHeatmap)
library(colorspace)
library(circlize)
library(dplyr)
library(cluster)
library(tidyverse)
library(vegan)
# set file variables ####
file_name <- "mPaintFDA_N1_N2_N3_N4_N5_N6_N7_N8"
integrate_state <- "integrated"
redu_state <- "redu"
annot_feats_disc <- c("PCA_NN", "UMAP_NN", "Batch", "Compound_ID")
pca_nn_cols <- scales::hue_pal()(11)
umap_nn_cols <- scales::hue_pal()(18)
batch_cols <- lighten(c(viridis(8)), amount = 0.3)
comp_id_cols <- c(Library = "lightgrey", setNames(scales::hue_pal()(10),pos_ctrls))
col_scale <- c("blue", "white", "red")
plot_width <- 20
plot_height <- 15
avg_profile <- TRUE
incl_dmso <- FALSE
pos_ctrls <- c("SodiumArsenite", "ROT", "Oligomycin", "Nocodazole", "MLN4924",
               "MitoQ", "CytochalasinD", "Chloroquine", "CCCP", "AntimycinA")
# set function to load data ####
load_data <- function(file_name,
                      integrate_state,
                      redu_state) {
  # load data
  if (avg_profile == TRUE) {
    data <- as.data.frame(
      fread(
        paste(
          "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_avg_data.csv", sep = ""),
        header = TRUE
      )
    )
  } else {
    data <- as.data.frame(
      fread(
        paste(
          "data/processed/", file_name, "_data_", integrate_state, "_", redu_state, ".csv", sep = ""),
        header = TRUE
      )
    )
  }
  # keep rownames
  rownames(data) <- data$V1
  data$V1 <- NULL
  # load meta
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
          "data/processed/", file_name, "_meta_", integrate_state, ".csv", sep = ""),
        header = TRUE
      )
    )
  }
  # keep rownames
  rownames(meta) <- meta$V1
  meta$V1 <- NULL 
  return(list(
    data = data,
    meta = meta
  ))
}
# run function to load data ####
data <- load_data(file_name, integrate_state, redu_state)
# if incl_dmso == FALSE, exclude DMSO wells ####
if (incl_dmso == FALSE) {
  data$meta <- data$meta[!data$meta$Compound %in% c("DMSO"),]
  data$data <- data$data[rownames(data$data) %in% rownames(data$meta),]
}
# populate Compound_ID if an annotation feature
if ("Compound_ID" %in% annot_feats_disc) {
  data$meta <- data$meta |>
    dplyr::mutate(
      Compound_ID = dplyr::case_when(
        Compound %in% pos_ctrls ~ Compound,
        Compound == "DMSO" ~ "DMSO",
        TRUE ~ "Library"
      )
    )
}
# compute pearson correlation matrix ####
heatmap_matrix <- cor(
  t(data$data),
  method = "pearson",
  use = "pairwise.complete.obs"
)
# build annotation ####
annot_df <- data$meta[
  , colnames(data$meta) %in% c(annot_feats_disc),
  drop = FALSE]
annot_df_stats <- annot_df
annot_df_plot  <- annot_df
# annotation df
annot_df_plot <- as.data.frame(annot_df_plot)
rownames(annot_df_plot) <- paste(rownames(data$meta))
annot_cols <- list()
# PCA_NN cols
pca_nn_levs <- sort(unique(as.character(data$meta$PCA_NN)))
annot_cols$PCA_NN <- setNames(
  pca_nn_cols[seq_along(pca_nn_levs)],
  pca_nn_levs
)
# UMAP_NN cols
umap_nn_levs <- sort(unique(as.character(data$meta$UMAP_NN)))
annot_cols$UMAP_NN <- setNames(
  umap_nn_cols[seq_along(umap_nn_levs)],
  umap_nn_levs
)
# Batch cols
batch_levs <- sort(unique(as.character(data$meta$Batch)))
annot_cols$Batch <- setNames(
  batch_cols[seq_along(batch_levs)],
  batch_levs
)
# comp_id cols
comp_id_levs <- unique(as.character(data$meta$Compound_ID))
# Put Library first if present
if ("Library" %in% comp_id_levs) {
  comp_id_levs <- c("Library", setdiff(comp_id_levs, "Library"))
}
annot_cols$Compound_ID <- comp_id_cols[comp_id_levs]
# full heatmap annotation object
heatmap_annot <- HeatmapAnnotation(
  df = annot_df_plot,
  col = annot_cols,
  annotation_name_gp = gpar(fontsize = 14),
  annotation_legend_param = list(
    Batch = list(
      title_gp = gpar(fontsize = 14, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    ),
    PCA_NN = list(
      title_gp = gpar(fontsize = 14, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    ),
    UMAP_NN = list(
      title_gp = gpar(fontsize = 14, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    ),
    Compound_ID = 
      list(
        title_gp = gpar(fontsize = 14, fontface = "bold"),
        labels_gp = gpar(fontsize = 14)
      )
  )
)
# plot heatmap ####
# clustering dendrogram (hierarchical clustering) using correlation distance #
# correlation distance measures the dissimilarity based on linear relationship between data points
# d = 1 - r, where r is the Pearson correlation coeffiecient
plot <- Heatmap(
  heatmap_matrix,
  name = "Pearson r",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  top_annotation = heatmap_annot,
  show_row_names = FALSE,
  show_column_names = FALSE,
  clustering_distance_rows = as.dist(1 - heatmap_matrix),
  clustering_distance_columns = as.dist(1 - heatmap_matrix),
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  rect_gp = gpar(col = NA),
  show_row_dend = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 14))
)
draw(plot)
# create function to calculate stats from correlation distance ####
calc_heatmap_stats <- function(heatmap_matrix,
                               annot_df,
                               annot_feats_disc
                               n_permanova = 999) {
  
  # make sure annotation rows line up with matrix rows
  annot_df <- annot_df[rownames(heatmap_matrix), , drop = FALSE]
  dist_obj <- as.dist(1 - heatmap_matrix)
  # calculate silhouette score for annot_feats from correlation distance
  # silhouette score quantifies if there are discrete clusters for the given grouping variable
  # value close to 1 = variable explains a high degree of separation, close to 0 = does not explain separation
  silhouette_df <- map_dfr(
    c(annot_feats_disc),
    function(feat) {
      group <- factor(annot_df[[feat]])
      # factor group annot_feats
      if (nlevels(group) < 2) {
        return(tibble(
          feature = feat,
          mean_silhouette = NA_real_,
          n_groups = nlevels(group)
        ))
      }
      # calculate silhouette score
      sil <- silhouette(
        as.integer(group),
        dist_obj
      )
      # populate mean_silhouette in placeholder table
      tibble(
        feature = feat,
        mean_silhouette = mean(sil[, 3], na.rm = TRUE),
        n_groups = nlevels(group)
      )
    }
  )
  # calculate PERMANOVA for annot_feats from correlation distance
  # PERMANOVA tests the variance explained by the given relationship
  # high R2 (1 = 100%) of variance explained by given variable
  # F statistic compares variation between groups and variation within groups
  # large F means the groups are more separated relative to the spread within groups
  # low p value = statistically significant
  permanova_df <- map_dfr(
    c(annot_feats_disc),
    function(feat) {
      ann <- annot_df
      # discrete vars should be factors
      if (feat %in% annot_feats_disc) {
        ann[[feat]] <- factor(ann[[feat]])
      }
      # calculate PERMANOVA
      perm <- vegan::adonis2(
        as.formula(paste("dist_obj ~", feat)),
        data = ann,
        permutations = n_permanova
      )
      # populate PERMANOVA values into a table
      tibble(
        feature = feat,
        variable_type = ifelse(feat %in% annot_feats_disc, "discrete", "continuous"),
        R2 = perm$R2[1],
        `F` = perm$`F`[1],
        p_value = perm$`Pr(>F)`[1]
      )
    }
  )
  # combine silhouette score and PERMANOVA stats into a single table
  combined_df <- left_join(
    silhouette_df,
    permanova_df,
    by = "feature"
  )
  return(
    combined_df
  )
}
# run function to calculate stats from correlation distance ####
stats <- calc_heatmap_stats(
  heatmap_matrix,
  annot_df_stats,
  annot_feats_disc,
  annot_feats_cont,
  n_permanova = 999
)
# save heatmap ####
pdf(paste("outputs/figures/", file_name, "_", "similarity_heatmap.pdf", sep = ""),
    width = plot_width,
    height = plot_height,
    useDingbats = FALSE)
draw(plot)
dev.off()
# save stats ####
write.csv(stats, 
          paste("outputs/data/", file_name, "_", "heatmap_stats.csv", sep = ""))
rm(list = ls())
