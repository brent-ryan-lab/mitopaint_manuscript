# Title: similarity heatmap mitopaint vis (mean per well) v1
# Step: 8
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 12-08-2026

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
file_name <- "mPaintDR2_N2_N3_N4"
integrate_state <- "integrated"
redu_state <- "redu"
pastel_cols <- lighten(c("#440154FF", "#238A8DFF", "#FDE725FF"), amount = 0.3)
annot_feats_disc <- c("Compound", "Batch")
annot_feats_cont <- c("Concentration")
col_scale <- c("blue", "white", "red")
agg_well <- TRUE
# set function to load data ####
load_data <- function(file_name,
                      integrate_state,
                      redu_state) {
  # load data
  data <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_data_", integrate_state, "_", redu_state, ".csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames
  rownames(data) <- data$V1
  data$V1 <- NULL
  # load meta
  meta <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_meta_", integrate_state, ".csv", sep = ""),
      header = TRUE
    )
  )
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
# set function to agg_well ####
# if agg_well = TRUE, then technical replicate wells within batch get averaged together
avg <- function(obj) {
  meta <- obj$meta
  dat  <- obj$data
  # combine meta + data so the grouping variables stay attached to each well
  combined <- cbind(meta, dat)
  # average replicate wells within each Batch x Condition group
  summary_df <- combined |>
    group_by(Batch, Condition) |>
    summarise(
      Compound = first(Compound),
      Concentration = first(Concentration),
      across(
        all_of(colnames(dat)),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) |>
    select(
      Compound,
      Concentration,
      Batch,
      Condition,
      all_of(colnames(dat))
    )
  # create output meta and data
  meta_out <- summary_df |>
    select(Compound, Concentration, Batch, Condition)
  meta_out <- as.data.frame(meta_out)
  data_out <- summary_df |>
    select(all_of(colnames(dat)))
  data_out <- as.data.frame(data_out)
  # row names must match and be Condition_Batch
  rn <- paste(meta_out$Condition, meta_out$Batch, sep = "_")
  rownames(meta_out) <- rn
  rownames(data_out) <- rn
  return(list(
    data = data_out,
    meta = meta_out
  ))
}
# run function to agg_well
if (agg_well) {
  data <- avg(data)
} else {
  tata <- data
}
# compute pearson correlation matrix ####
heatmap_matrix <- cor(
  t(data$data),
  method = "pearson",
  use = "pairwise.complete.obs"
)
# build annotation ####
annot_df <- data$meta[
  , colnames(data$meta) %in% c(annot_feats_disc, annot_feats_cont),
  drop = FALSE]
annot_df <- as.data.frame(annot_df)
rownames(annot_df) <- paste(annot_df$Compound, annot_df$Concentration, annot_df$Batch, sep = "_")
make_disc_cols <- function(meta, var, palette) {
  levs <- sort(unique(as.character(meta[[var]])))
  if ("DMSO" %in% levs) {
    levs <- c("DMSO", setdiff(levs, "DMSO"))
  }
  setNames(
    palette[seq_along(levs)],
    levs
  )
}
  annot_cols <- list()
  # discrete annotation features (vector)
  for (feat in annot_feats_disc) {
    levs <- sort(unique(as.character(data$meta[[feat]])))
    if ("DMSO" %in% levs) {
      levs <- c("DMSO", setdiff(levs, "DMSO"))
    }
    annot_cols[[feat]] <- setNames(
      pastel_cols[seq_along(levs)],
      levs
    )
  }
  # continuous annotation features (function)
  for (feat in annot_feats_cont) {
    annot_cols[[feat]] <- circlize::colorRamp2(
      range(data$meta[[feat]], na.rm = TRUE),
      c("white", "black")
    )
  }
  
  heatmap_annot <- HeatmapAnnotation(
    df = annot_df,
    col = annot_cols,
    annotation_name_gp = gpar(fontsize = 14),
    annotation_legend_param = list(
      Compound = list(
        title_gp = gpar(fontsize = 14, fontface = "bold"),
        labels_gp = gpar(fontsize = 14)
      ),
      Batch = list(
        title_gp = gpar(fontsize = 14, fontface = "bold"),
        labels_gp = gpar(fontsize = 14)
      ),
      Concentration = list(
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
                               annot_feats_disc,
                               annot_feats_cont,
                               n_permanova = 999) {
  
  # make sure annotation rows line up with matrix rows
  annot_df <- annot_df[rownames(heatmap_matrix), , drop = FALSE]
  dist_obj <- as.dist(1 - heatmap_matrix)
  # calculate silhouette score for annot_feats from correlation distance
  # silhouette score quantifies if there are discrete clusters for the given grouping variable
  # value close to 1 = variable explains a high degree of separation, close to 0 = does not explain separation
  silhouette_df <- map_dfr(
    c(annot_feats_disc, annot_feats_cont),
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
    c(annot_feats_disc, annot_feats_cont),
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
stats <- calc_heatmap_stats(heatmap_matrix,
                            annot_df,
                            annot_feats_disc,
                            annot_feats_cont,
                            n_permanova = 999)
# save heatmap ####
pdf(paste("outputs/figures/", file_name, "_", "similarity_heatmap.pdf", sep = ""),
    width =7.4,
    height =6,
    useDingbats = FALSE)
draw(plot)
dev.off()
# save stats ####
write.csv(stats, 
          paste("outputs/data/", file_name, "_", "heatmap_stats.csv", sep = ""))
rm(list = ls())