# Title: remove redundant mitopaint data (N=3) (mean per well) v1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 15-06-2026

# load packages ####
library(data.table)
# set variables ####
file_name <- "mPaintDR2_N2_N3_N4"
integrate_state <- "integrated"
dataset_name <- "mPaintDR2_N2_N3_N4"
cor_thresh <- 0.95
var_tol <- 1e-12
excl_feats <- c("Nucleus", "Nuclei", "mTagBFP2")
# load data ####
# load data as df
df <- as.data.frame(
  fread(
    paste(
      "data/processed/", file_name, "_data_", integrate_state ,".csv", sep = ""), 
    header = TRUE)
)
# keep rownames as WELL_BATCH
rownames(df) <- df$V1
df$V1 <- NULL
# create function to remove redundant features ####
# this function removes: 
# features with variance below tolerance threshold
# feature names with certain char strings in the name (eg nuclear features)
# features with >0.95 pearson correlation to an existing feature
rem_redundant <- function(df,
                          dataset_name = NULL,
                          excl_feats = NULL,
# default assumes features with variance below 0.000000000001
# have negligable variance contributed to data structure
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
    if (is.null(dataset_name)) {
      dataset_name <- "Dataset"
    }
    cat("Dataset:", dataset_name, "\n",
        "Original features:", original_n, "\n",
        "Removed by pattern:", removed_pattern, "\n",
        "Removed by variance:", removed_variance, "\n",
        "Removed by correlation (>|", cor_thresh, "|):", removed_cor, "\n",
        "Final features retained:", ncol(data), "\n"
    )
  }
  return(data)
}

# run function to remove redundant features ####
redu <- rem_redundant(df,
                      dataset_name,
                      excl_feats,
                      var_tol,
                      cor_thresh)

# save data ####
# save redu data in /data/processed
write.csv(redu,
          paste(
            "data/processed/", file_name, "_data_redu.csv", sep = "")
)
rm(list = ls())