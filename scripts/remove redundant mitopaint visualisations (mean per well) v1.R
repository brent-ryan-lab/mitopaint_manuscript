# Title: remove redundant mitopaint visualisations (mean per well) v1
# Step: 5.1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 03-08-26

# load packages ####
library(data.table)
library(ggplot2)
library(ggpubr)
# set variables ####
file_name <- "mPaintDR2_N2_N3_N4"
integrate_state <- "integrated"
dataset_name <- "mPaintDR2_N2_N3_N4"
cor_thresh <- 0.95
cor_thresh_range_scree <- seq(from = 0.90, to = 1, by = 0.01)
cor_thresh_range_pca <- seq(from = 0.65, to = 1, by = 0.05)
var_tol <- 1e-12
excl_feats <- c("Nucleus", "Nuclei", "mTagBFP2")
# load data ####
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
# create function to remove redundant features ####
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
  keep_variance <- variances > var_tol
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
# run function to remove redundant features (loop range) ####
scree_df <- list()
for (ct in cor_thresh_range_scree) {
  
  # Format cutoff cleanly
  ct_label <- sprintf("%.2f", ct)
  
  output_name <- paste0("redu_", ct_label)
  
  scree_df[[output_name]] <- rem_redundant(
    df,
    excl_feats = excl_feats,
    var_tol = 1e-12,
    output_name = output_name,
    cor_thresh = ct
  )
}
rm(output_name, ct, ct_label)
# plot elbow plot of n features x corr thresh ####
scree_plotting_df <- purrr::imap_dfr(
  scree_df,
  function(df, name) {
    
    tibble::tibble(
      ct = as.numeric(
        stringr::str_remove(name, "^redu_")
      ),
      features_retained = ncol(df)
    )
  }
)
plots <- list()
plots$scree <- ggplot(scree_plotting_df,
       aes(x = ct,
           y = features_retained)) +
  geom_point(color = "black", size = 1) +
  geom_line(colour = "black", linewidth = 0.6) +
  annotate("line",
           x = c(min(scree_plotting_df$ct), cor_thresh),
           y = as.integer(scree_plotting_df[scree_plotting_df$ct == cor_thresh,2]),
           colour = "grey",
           linetype = "dashed") +
  annotate("line",
           x = cor_thresh,
           y = c(0, as.integer(scree_plotting_df[scree_plotting_df$ct == cor_thresh,2])),
           colour = "grey",
           linetype = "dashed") +
  annotate("text",
           x = cor_thresh + 0.005,
           y = as.integer(scree_plotting_df[scree_plotting_df$ct == cor_thresh,2]) - (0.04 * max(scree_plotting_df$features_retained)),
           label = as.integer(scree_plotting_df[scree_plotting_df$ct == cor_thresh,2]),
           size = 2.5) +
  theme_pubr() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        legend.position = "none") +
  scale_x_continuous(
    limits = c(min(cor_thresh_range_scree), max(cor_thresh_range_scree) + 0.01),
    breaks = seq(min(cor_thresh_range_scree), max(cor_thresh_range_scree), by = 0.02),
    labels = function(x) sprintf("%.2f", x),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, max(scree_plotting_df$features_retained) * 1.1),
    expand = c(0, 0)
  ) +
  labs(x = "Pearson Correlation Cutoff",
       y = "Number of Features Retained")
plots$scree
# create function to plot feature retained pie chart  ####
pie_df
# run function to plot feature retained pie chart (by channel) ####
# run function to plot feature retained pie chart (by region) ####
# create function to plot pca x corr thresh (loop) ####
# run function to plot pca x corr thresh (loop) ####
# save plots ####
plots$scree