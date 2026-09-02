# Title: remove redundant mitopaint vis (mean per well) v1
# Step: 5.1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 03-08-26

# load packages ####
library(data.table)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(tidyverse)
library(Seurat)
library(cowplot)
library(cluster)
library(purrr)
library(stringr)
library(tibble)
library(rlang)
library(viridis)
library(cluster)
# set variables ####
file_name <- "mPaintDR2_N2_N3_N4"
integrate_state <- "integrated"
dataset_name <- "mPaintDR2_N2_N3_N4"
cor_thresh <- 0.95
cor_thresh_range_scree <- seq(from = 0.90, to = 1, by = 0.01)
cor_thresh_range_pca <- seq(from = 0.65, to = 1, by = 0.05)
var_tol <- 1e-12
excl_feats <- c("Nucleus", "Nuclei", "mTagBFP2")
pastel_cols <- lighten(c("#440154FF", "#238A8DFF", "#FDE725FF"), amount = 0.3)
# load data ####
# initialise plots list
plots <- list()
# load df
df <- as.data.frame(
  fread(
    paste(
      "data/processed/", file_name, "_data_", integrate_state ,".csv", sep = ""), 
    header = TRUE)
)
# keep rownames as WELL_BATCH
rownames(df) <- df$V1
df$V1 <- NULL
# load meta
meta <- as.data.frame(
  fread(
    paste(
      "data/processed/", file_name, "_meta_", integrate_state ,".csv", sep = ""), 
    header = TRUE)
)
# keep rownames as WELL_BATCH
rownames(meta) <- meta$V1
meta$V1 <- NULL
pastel_cols_5 <- lighten(viridis(n=5), amount = 0.3)
# create function to remove redundant features ####
# this function removes: 
# features with variance below tolerance threshold
# feature names with certain char strings in the name (eg nuclear features)
# features with >0.95 pearson correlation to an existing feature
rem_redundant <- function(df,
                          output_name = NULL,
                          excl_feats = NULL,
                          # default assumes features with variance below 0.000000000001
                          # have negligible variance contributed to data structure
                          var_tol = 1e-12,
                          # default assumes 0.95 pearson correlation threshold
                          # ie. features with >90% shared variance or are near identical, are redundant
                          cor_thresh = 0.95,
                          verbose = TRUE) {
  data <- df
  original_n <- ncol(data)
  # remove unwanted name patterns
  excl_pattern <- paste(excl_feats, collapse = "|")
  keep_pattern <- !grepl(excl_pattern, colnames(data), ignore.case = TRUE)
  data <- data[, keep_pattern, drop = FALSE]
  removed_pattern <- sum(!keep_pattern)
  # remove low variance
  variances <- apply(data, 2, var, na.rm = TRUE)
  # also remove features with NA variance
  # NA = variance unable to be computed due to too many NA values
  keep_variance <- !is.na(variances) & variances > var_tol
  data <- data[, keep_variance, drop = FALSE]
  removed_variance <- sum(!keep_variance)
  # remove highly correlated features
  removed_cor <- 0
  if (ncol(data) > 1) {
    cor_mat <- cor(data, use = "pairwise.complete.obs")
    cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
    high_cor <- which(abs(cor_mat) > cor_thresh, arr.ind = TRUE)
    if (nrow(high_cor) > 0) {
      remove_cols <- unique(colnames(cor_mat)[high_cor[,2]])
      removed_cor <- length(remove_cols)
      data <- data[, !colnames(data) %in% remove_cols, drop = FALSE]
    }
  }
  # dynamically print output
  # just set verbose = FALSE, to skip printing below summary outputs into console
  if (verbose) {
    if (is.null(output_name)) {
      output_name <- "Dataset"
    }
    cat("Dataset:", output_name, "\n",
        "Original features:", original_n, "\n",
        "Removed by pattern:", removed_pattern, "\n",
        "Removed by variance:", removed_variance, "\n",
        "Removed by correlation (>|", cor_thresh, "|):", removed_cor, "\n",
        "Final features retained:", ncol(data), "\n"
    )
  }
  return(data)
}
# run function to remove redundant features (by cor_thresh_range_scree) ####
scree_df <- list()
# loop removes redundant features for range of pearson corr. threshold cutoffs
for (ct in cor_thresh_range_scree) {
  # format cutoff cleanly (each value label is at 2 decimal places)
  ct_label <- sprintf("%.2f", ct)
  output_name <- paste0("redu_", ct_label)
  scree_df[[output_name]] <- rem_redundant(
    df,
    # remove feature names with certain char strings in the name (eg nuclear features)
    excl_feats = excl_feats,
    # remove features with variance below tolerance threshold
    var_tol = 1e-12,
    output_name = output_name,
    # remove features with >ct pearson correlation (eg. ct = 0.95) to an existing feature
    cor_thresh = ct
  )
}
# remove temp variables
rm(output_name, ct, ct_label)
# plot elbow plot of n features x corr thresh ####
# scree_plotting_df makes a simple data frame of the redundancy cutoff (ct) vs. number of features retained
scree_plotting_df <- imap_dfr(
  scree_df,
  function(df, name) {
    tibble(
      # ct cutoff value (as a number to plot on continuous scale)
      ct = as.numeric(
        str_remove(name, "^redu_")
      ),
      # number of features retained (= columns in scree_df)
      features_retained = ncol(df)
    )
  }
)
plots$scree <- ggplot(scree_plotting_df,
       aes(x = ct,
           y = features_retained)) +
  geom_point(color = "black", size = 1) +
  geom_line(colour = "black", linewidth = 0.6) +
  # horizontal annotation line for cor_thresh (eg. 0.95)
  annotate("line",
           x = c(min(scree_plotting_df$ct), cor_thresh),
           y = as.integer(scree_plotting_df[scree_plotting_df$ct == cor_thresh,2]),
           colour = "grey",
           linetype = "dashed") +
  # vertical annotation line for cor_thresh (eg. 0.95)
  annotate("line",
           x = cor_thresh,
           y = c(0, as.integer(scree_plotting_df[scree_plotting_df$ct == cor_thresh,2])),
           colour = "grey",
           linetype = "dashed") +
  # annotation text for number of retained features @ cor_thresh (eg. 0.95)
  annotate("text",
           x = cor_thresh + 0.005,
           y = as.integer(scree_plotting_df[scree_plotting_df$ct == cor_thresh,2]) - 
             (0.04 * max(scree_plotting_df$features_retained)),
           label = as.integer(scree_plotting_df[scree_plotting_df$ct == cor_thresh,2]),
           size = 2.5) +
  # tidying theme
  theme_pubr() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        legend.position = "none") +
  # manual x scale to avoid cutting off line
  scale_x_continuous(
    limits = c(min(cor_thresh_range_scree), max(cor_thresh_range_scree) + 0.01),
    breaks = seq(min(cor_thresh_range_scree), max(cor_thresh_range_scree), by = 0.02),
    labels = function(x) sprintf("%.2f", x),
    expand = c(0, 0)
  ) +
  # manual y scale to avoid cutting off line
  scale_y_continuous(
    limits = c(0, max(scree_plotting_df$features_retained) * 1.1),
    expand = c(0, 0)
  ) +
  labs(x = "Pearson Correlation Cutoff",
       y = "Number of Features Retained")
plots$scree
# create function to plot feature retained pie chart  ####
plot_feature_pie <- function(obj,
                             patterns,
                             plot_title = "",
                             colours) {
  # df_redu counts how many features are of each class 
  # eg. n of retained features per channel/ region
  df_redu <- data.frame(
    feature_name = colnames(obj)
  ) %>%
    mutate(
      feature_type = case_when(
        !!!imap(patterns, 
                       ~ expr(str_detect(feature_name, !!.x) ~ !!.y)),
        TRUE ~ "Other"
      )
    ) %>%
    count(feature_type, name = "counts")
  # df_redu_2 calculates the plotting order in the pie chart
  df_redu_2 <- df_redu %>% 
    mutate(
      csum = rev(cumsum(rev(counts))),
      pos = counts/2 + lead(csum, 1),
      pos = if_else(is.na(pos), counts/2, pos)
    ) %>%
    arrange(csum)
  # match df_redu and df_redu_2
  df_redu <- df_redu %>%
    slice(match(df_redu_2$feature_type, feature_type))
  # reverse plotting order
  df_redu$feature_type <- factor(
    df_redu$feature_type,
    levels = rev(df_redu$feature_type)
  )
  # build pie
  pie_plot <- ggplot(df_redu,
                     aes(x = "", y = counts, fill = feature_type)) +
    geom_bar(stat = "identity", width = 1, 
             colour = "black",
             linewidth = 0.2) +
    coord_polar("y", start = 0) +
    # add label of counts in the centre of pie slice
    geom_label(aes(label = counts),
               position = position_stack(vjust = 0.5),
               colour = "black",
               fill = "white",
               size = 3,
               show.legend = FALSE) +
    # add label of class (eg. channel/ region) outside pie slice
    geom_label(data = df_redu_2,
               aes(y = pos,
                   label = gsub(" ", "\n", feature_type)),
               fill = "white",
               colour = "black",
               size = 3,
               nudge_x = 1.1,
               show.legend = FALSE) +
    # colour fill based on given colours
    scale_fill_manual(values = colours) +
    # tidying theme
    labs(x = NULL, y = NULL, fill = NULL,
         # title based on given title
         title = plot_title) +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        margin = margin(b = -16)
      )
    )
  return(pie_plot)
}
# run function to plot feature retained pie chart (by channel) ####
plots$channel_pie <- plot_feature_pie(
  # uses cor_thresh df (eg. redu_0.95)
  # this is under the assumption that is was already calculate in scree_df
  scree_df[[paste("redu_", cor_thresh, sep = "")]],
  patterns = c(
    "mTagBFP2" = "mTagBFP2",
    "CellROX" = "CellRox",
    "mt-Keima pH7" = "mt-Keima pH7",
    "mt-Keima pH4" = "mt-Keima pH4",
    "TMRM" = "TMRM test"
  ),
  plot_title = "Features Retained by Channel",
  colours = c("TMRM" = "#FFB000", "CellROX" = "#DC267F", "mt-Keima pH7" = "#23CC86", "mt-Keima pH4" = "#FE6100", "Other" = "grey70")
)
plots$channel_pie
# run function to plot feature retained pie chart (by region) ####
plots$region_pie <- plot_feature_pie(
  # uses cor_thresh df (eg. redu_0.95)
  # this is under the assumption that is was already calculate in scree_df
  scree_df[[paste("redu_", cor_thresh, sep = "")]],
  patterns = c(
    "Ring Region" = "Ring Region",
    "Cytoplasm Region" = "Cytoplasm Region",
    "Membrane Region" = "Membrane Region",
    "Cell Region" = "Cell Region"
  ),
  plot_title = "Features Retained by Region",
  colours = pastel_cols_5[c(1,2,3,5)]
)
plots$region_pie
# plot pca x corr thresh (by cor_thresh_range_pca) ####
scree_df <- list()
# loop overwrite scree_df using _pca range instead of _scree range
for (ct in cor_thresh_range_pca) {
  # format cutoff cleanly (each value label is at 2 decimal places)
  ct_label <- sprintf("%.2f", ct)
  output_name <- paste0("redu_", ct_label)
  scree_df[[output_name]] <- rem_redundant(
    df,
    # remove feature names with certain char strings in the name (eg nuclear features)
    excl_feats = excl_feats,
    # remove features with variance below tolerance threshold
    var_tol = 1e-12,
    output_name = output_name,
    # remove features with >ct pearson correlation (eg. ct = 0.95) to an existing feature
    cor_thresh = ct
  )
}
# remove temp variables
rm(output_name, ct, ct_label)
# run_pca function does following:
# 1. put data frame and meta into a Seurat
# 2. run PCA dim red (on scaled data)
# 3. save PCA embeddings
# 4. calculate variance explained (for axis labels)
run_pca <- function(feature_matrix,
                    metadata,
                    assay_name = "MP",
                    dims_use = 1:2,
                    k_param = 15,
                    resolution = 1,
                    seed = 42) {
  # keep only rows where Compound is in the selected conditions 
  keep_rows <- metadata$Compound %in% c("DMSO", "CCCP", "ROT")
  metadata <- metadata[keep_rows, , drop = FALSE]
  feature_matrix <- feature_matrix[rownames(metadata), , drop = FALSE]
  
  # put data frame and meta into a Seurat
  seurat_obj <- CreateSeuratObject(
    counts = t(as.matrix(feature_matrix)),
    meta.data = metadata,
    assay = assay_name
  )
  DefaultAssay(seurat_obj) <- assay_name
  seurat_obj <- SetAssayData(
    seurat_obj,
    assay = assay_name,
    layer = "data",
    new.data = GetAssayData(
      seurat_obj,
      assay = assay_name,
      layer = "counts"
    )
  )
  # run PCA dim red
  seurat_obj <- ScaleData(seurat_obj)
  
  seurat_obj <- RunPCA(
    seurat_obj,
    features = rownames(seurat_obj),
    seed.use = seed
  )
  # save PCA embeddings
  pca_embeddings <- as.data.frame(
    Embeddings(seurat_obj, "pca")
  )
  # calculate variance explained (for axis labels)
  pc_sd <- seurat_obj[["pca"]]@stdev
  pc_var <- pc_sd^2
  pca_var <- data.frame(
    PC = paste0("PC", seq_along(pc_sd)),
    SD = pc_sd,
    Variance = pc_var,
    Percent_Variance = 100 * pc_var / sum(pc_var),
    Cumulative_Percent_Variance = cumsum(100 * pc_var / sum(pc_var))
  )
  return(list(
    seurat = seurat_obj,
    pca_embeddings = pca_embeddings,
    pca_var = pca_var
  ))
}
# vectorise the list of redu df to loop through (based on cor_thresh_range_pca)
redu_names <- grep("^redu_", names(scree_df), value = TRUE)
# lapply() loop run_pca through range, and set names as redu_names
integrated_redu_ct <- setNames(
  lapply(redu_names, function(nm) {
    run_pca(
      feature_matrix = scree_df[[nm]],
      metadata = meta
    )
  }),
  sub("^redu_", "ct_", redu_names)
)
# calculated pca metrics are stored in stats list
stats <- list()
# calculate silhouette score between compound treatments
# this assumes DMSO, CCCP, and ROT are present in the dataset
stats$pca_silhouette <- {
  compound_levels <- c("DMSO", "CCCP", "ROT")
  sil_scores <- sapply(
    names(integrated_redu_ct),
    function(ct) {
      obj <- integrated_redu_ct[[ct]]
      # pca embedding matrix
      emb <- obj$pca_embeddings[, 1:2, drop = FALSE]  # use first 10 PCs
      # compound labels
      comp <- obj$seurat@meta.data$Compound
      comp <- factor(comp, levels = compound_levels)
      # compute silhouette score
      # the closer the values are to 1, the more the groups are relatively separated
      sil <- silhouette(
        as.numeric(comp),
        dist(emb)
      )
      mean(sil[, 3])
    }
  )
  # dataframe for annotation
  data.frame(
    Cutoff = names(integrated_redu_ct),
    Mean_Silhouette = sil_scores
  )
}
# calculate bw ratio between compound treatments
# this assumes DMSO, CCCP, and ROT are present in the dataset
stats$pca_bw_ratio <- {
  compound_levels <- c("DMSO", "CCCP", "ROT")
  bw_vals <- sapply(
    names(integrated_redu_ct),
    function(ct) {
      obj <- integrated_redu_ct[[ct]]
      # uses only PC1 and 2
      df <- obj$pca_embeddings[, c("PC_1", "PC_2"), drop = FALSE]
      meta <- obj$seurat@meta.data
      meta$Compound <- factor(meta$Compound,
                              levels = compound_levels)
      # overall centroid is the variance of all points in the dataset combined
      # ie. not separated by compound
      overall_mean <- colMeans(df)
      # calculate between variance
      # this corresponds to how different each condition is to other conditions
      group_means <- aggregate(df,
                               by = list(meta$Compound),
                               FUN = mean)
      group_sizes <- table(meta$Compound)
      between_var <- sum(
        sapply(seq_len(nrow(group_means)), function(i) {
          ni <- group_sizes[i]
          sum((as.numeric(group_means[i, -1]) - overall_mean)^2) * ni
        })
      ) / sum(group_sizes)
      # calculate within variance
      # this corresponds to how consistent the points are within a condition
      within_var <- mean(
        sapply(split(df, meta$Compound), function(group_df) {
          group_center <- colMeans(group_df)
          mean(apply(group_df, 1, function(row)
            sum((row - group_center)^2)))
        })
      )
      # bw ratio
      # higher bw ratio = more relative separation
      between_var / within_var
    }
  )
  # dataframe for annotation
  data.frame(
    Cutoff = names(integrated_redu_ct),
    BW_Ratio_PC12 = bw_vals
  )
}
# pca ct grid plots pca and adds annotation for stats in cowplot grid of entire range
# this is to visualise the loss is separation by drug as the threshold goes too low
make_pca_cutoff_grid <- function(integrated_redu_ct,
                                 ncol = 4,
                                 nrow = NULL) {
  # assumes DMSO, CCCP and ROT are in the dataset
  compound_levels <- c("DMSO", "CCCP", "ROT")
  cutoff_names <- names(integrated_redu_ct)
  cutoff_vals <- as.numeric(sub("^ct_", "", cutoff_names))
  ordered_names <- cutoff_names[order(cutoff_vals, decreasing = TRUE)]
  plot_list <- lapply(ordered_names, function(ct) {
    obj <- integrated_redu_ct[[ct]]
    df <- obj$pca_embeddings %>%
      cbind(obj$seurat@meta.data) %>%
      mutate(Compound = factor(Compound,
                               levels = compound_levels))
    # pull ct cutoff used
    cutoff_val <- sub("^ct_", "", ct)
    # pull count of n retained features
    n_features <- nrow(obj$seurat)
    # pull PC variance explained 
    pc1_var <- round(obj$pca_var$Percent_Variance[1], 2)
    pc2_var <- round(obj$pca_var$Percent_Variance[2], 2)
    x_lab <- paste0("PC_1 (", pc1_var, "%)")
    y_lab <- paste0("PC_2 (", pc2_var, "%)")
    # pull between:within groups variance ratio
    bw_df <- stats$pca_bw_ratio
    cutoff_val <- sub("^ct_", "", ct)
    n_features <- nrow(obj$seurat)
    bw_ratio <- round(
      bw_df$BW_Ratio[bw_df$Cutoff == ct],
      2
    )
    # combine all above pulled data together into title text
    title_text <- paste0(
      "Redundancy cutoff: ", cutoff_val,
      "\nFeatures retained: ", n_features,
      "\nBetween:Within ratio: ", bw_ratio
    )
    # pca dot plot
    ggplot(df, aes(x = PC_1, y = PC_2, color = Compound)) +
      geom_point(shape = 16, size = 1.5, alpha = 0.5) +
      # add ellipse around each compound group for easier visualisation
      stat_ellipse(
        aes(group = Compound),
        type = "norm",
        linewidth = 0.8
      ) +
      # color using pastel_cols for compound
      scale_color_manual(values = pastel_cols) +
      # tidy theme
      theme_pubr() +
      theme(
        plot.title = element_text(hjust = 0.5,
                                  size = 9,
                                  face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        panel.grid = element_blank()
      ) +
      labs(title = title_text,
           x = x_lab,
           y = y_lab)
  }
  )
  # cowplot grid
  plot_grid(plotlist = plot_list,
            ncol = ncol,
            nrow = nrow)
}
# run function and keep final plot in plots list
plots$pca_cutoff_grid <- make_pca_cutoff_grid(
  integrated_redu_ct,
  ncol = 4,
  nrow = NULL
)
plots$pca_cutoff_grid
# save stats ####
write.csv(
  purrr::reduce(
    stats,
    dplyr::full_join,
    by = "Cutoff"
  ),
  paste0("outputs/data/", file_name, "_", integrate_state, "_redu_stats.csv"),
  row.names = FALSE
)
# save plots ####
ggsave(paste(
  "outputs/figures/", file_name, "_", integrate_state, "_redu_curve.pdf", sep = ""),
       plots$scree,
       width = 3.6,
       height = 3.6,
       units = "in",
       dpi = 300)
ggsave(paste(
  "outputs/figures/", file_name, "_", integrate_state, "_redu_channel_pie.pdf", sep = ""),
       plots$channel_pie,
       width = 3.6,
       height = 3.6,
       units = "in",
       dpi = 300)
ggsave(paste(
  "outputs/figures/", file_name, "_", integrate_state, "_redu_region_pie.pdf", sep = ""),
  plots$region_pie,
  width = 3.6,
  height = 3.6,
  units = "in",
  dpi = 300)
ggsave(paste(
  "outputs/figures/", file_name, "_", integrate_state, "_pca_ct_cutoff_grid.pdf", sep = ""),
       plots$pca_cutoff_grid,
       width = 9.4,
       height = 5.3,
       units = "in",
       dpi = 300)
rm(list = ls())