# Title: plate drift correct mitopaint data (mean per well) v1
# Step: 2
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 18-06-2026

# load packages ####
library(data.table)
library(tidyverse)
library(purrr)
# set variables ####
batches_info <- list(
  N1A = list(
    file_name = "SF260604_mPaintSpace2_N1A",
    dmso_wells = c("2_2","3_9","4_2","5_9","6_2",
                    "7_9","8_2","9_9","10_2","11_9",
                    "12_2","13_9","14_2","15_9",
                    "16_2","17_9","18_2","19_9",
                    "20_2","21_9","22_2","23_9")
  ),
  N1B = list(
    file_name = "SF260604_mPaintSpace2_N1B",
    dmso_wells = c("2_2","3_9","4_2","5_9","6_2",
                   "7_9","8_2","9_9","10_2","11_9",
                   "12_2","13_9","14_2","15_9",
                   "16_2","17_9","18_2","19_9",
                   "20_2","21_9","22_2","23_9")
),
    N2A = list(
    file_name = "SF260604_mPaintSpace2_N2A",
    dmso_wells = c("2_2","3_9","4_2","5_9","6_2",
                   "7_9","8_2","9_9","10_2","11_9",
                   "12_2","13_9","14_2","15_9",
                   "16_2","17_9","18_2","19_9",
                   "20_2","21_9","22_2","23_9")
  ),
 N2B = list(
    file_name = "SF260604_mPaintSpace2_N2B",
    dmso_wells = c("2_2","3_9","4_2","5_9","6_2",
                    "7_9","8_2","9_9","10_2","11_9",
                    "12_2","13_9","14_2","15_9",
                    "16_2","17_9","18_2","19_9",
                    "20_2","21_9","22_2","23_9")
  ),
  N3A = list(
    file_name = "SF260701_mPaintSpace2_N3A",
    dmso_wells = c("2_2","3_9","4_2","5_9","6_2",
                    "7_9","8_2","9_9","10_2","11_9",
                    "12_2","13_9","14_2","15_9",
                    "16_2","17_9","18_2","19_9",
                    "20_2","21_9","22_2","23_9")
  ),
  N3B = list(
    file_name = "SF260701_mPaintSpace2_N3B",
    dmso_wells = c("2_2","3_9","4_2","5_9","6_2",
                   "7_9","8_2","9_9","10_2","11_9",
                   "12_2","13_9","14_2","15_9",
                   "16_2","17_9","18_2","19_9",
                   "20_2","21_9","22_2","23_9")
)
)
p_sig <- 0.05
# create function to load data ####
load_data <- function(file_name, dmso_wells) {
  # load tidy data
  df <- as.data.frame(
    fread(
      paste("data/processed/", file_name, "_data_tidy.csv", sep = ""),
      header = TRUE
    )
  )
  rownames(df) <- df$V1
  df$V1 <- NULL
  # load tidy metadata
  meta <- as.data.frame(
    fread(
      paste("data/processed/", file_name, "_meta_tidy.csv", sep = ""),
      header = TRUE
    )
  )
  rownames(meta) <- meta$V1
  meta$V1 <- NULL
  return(list(
    file_name = file_name,
    dmso_wells = dmso_wells,
    df = df,
    meta = meta
  ))
}
# load data ####
batches <- map(
  batches_info,
  function(batch_info) {
    load_data(
      file_name = batch_info$file_name,
      dmso_wells = batch_info$dmso_wells
    )
  }
)
# create function to calculate plate drift ####
# function will calculate possible models of plate drift 
# 0 = intercept only, 1 = linear, 2 = poly, based on DMSO control wells
# then will ANOVA test which is the model of best fit for each feature
calc_drift <- function(df, meta, dmso_wells, p_sig) {
  # subset dmso wells
  meta_dmso <- meta[meta$Well %in% dmso_wells, , drop = FALSE]
  df_dmso <- df[rownames(meta_dmso), , drop = FALSE]
  # create vector of all features
  features <- colnames(df)
  # lapply() calculatea fit table for all features
  fit_table <- lapply(seq_along(features), function(i) {
    # y in each model is the feature value, the loop goes along the sequence of features
    y <- df_dmso[[i]]
    # x in each model calculation is the imaging order
    order <- meta_dmso$Order
    # calculate lm for intercept only (0) = no drift
    fit_0 <- lm(y ~ 1)
    # calculate lm for linear model (1) = straight drift
    fit_1 <- lm(y ~ order)
    # calculate lm for poly model (2) = curved drift
    # raw must = TRUE to use actual polynomial terms and not orthogonal terms
    fit_2 <- lm(y ~ poly(order, 2, raw = TRUE))
    # anova test to compare between each model
    an <- anova(fit_0, fit_1, fit_2)
    # output coefficients for each model
    data.frame(
      ID = i,
      feature = features[i],
      # intercept only coefficients
      intercept_0 = coef(fit_0)[1],
      # linear model coefficients
      intercept_1 = coef(fit_1)[1],
      gradient_1_1 = coef(fit_1)[2],
      # poly model coefficients
      intercept_2 = coef(fit_2)[1],
      gradient_2_1 = coef(fit_2)[2],
      gradient_2_2 = coef(fit_2)[3],
      # anova p values for comparing between models
      # asks: does adding a linear term improve model fit compared to intercept only?
      p_lin = an$`Pr(>F)`[2],
      # asks: does adding a poly term improve model fit compared to linear?
      p_quad = an$`Pr(>F)`[3]
    )
  })
  # each fit table row represents one feature (~3000 rows)
  fit_table <- do.call(rbind, fit_table)
  # choose model with best p
  fit_table_best <- fit_table |>
    # each row is handled independently
    dplyr::rowwise() |>
    dplyr::mutate(
      # use poly intercept if p is best, else if use lm intercept, else use 0
      best_intercept = dplyr::case_when(
        p_quad < p_sig ~ intercept_2,
        p_lin  < p_sig ~ intercept_1,
        TRUE ~ 0
      ),
      # use poly gradient 1 if p is best, else if use lm gradient 1, else use 0
      best_grad1 = dplyr::case_when(
        p_quad < p_sig ~ gradient_2_1,
        p_lin  < p_sig ~ gradient_1_1,
        TRUE ~ 0
      ),
      # use poly gradient 2 if p is best, else use 0
      best_grad2 = dplyr::case_when(
        p_quad < p_sig ~ gradient_2_2,
        TRUE ~ 0
      ),
      # gives p value of the best model for drift
      # asks: how well does the chosen model fit?
      best_p = dplyr::case_when(
        p_quad < p_sig ~ p_quad,
        p_lin  < p_sig ~ p_lin,
        TRUE ~ 0
      )
    ) |>
    dplyr::ungroup()
  # return table of all fits, and best fits for drift corr of each feature
  return(fit_table_best)
}
# run function to calculate plate drift ####
batches <- map(
  batches,
  function(batch_obj) {
    batch_obj$fits <- calc_drift(
      df = batch_obj$df,
      meta = batch_obj$meta,
      dmso_wells = batch_obj$dmso_wells,
      p_sig = p_sig
    )
    batch_obj
  }
)
# create function to apply plate drift correction ####
drift_corr <- function(df, meta, dmso_wells, fits) {
  # compute predicted drift based on best fit model
  # corr_table contains all the correction factors
  corr_table <- sapply(seq_len(ncol(df)), function(i) {
    fits$best_intercept[i] +
      fits$best_grad1[i] * meta$Order +
      fits$best_grad2[i] * (meta$Order^2)
  })
  # match colnames and rownames in corr table to df 
  # used to identify feature names (cols) and well names (rows)
  corr_table <- as.data.frame(corr_table)
  colnames(corr_table) <- colnames(df)
  rownames(corr_table) <- rownames(df)
  # corr_df contains the corrected data
  # it is temporarily populated with the original df data
  # these values eventually gets overwritten with the drift corrected data
  corr_df <- df
  # find the first DMSO well 
  # the first DMSO is a global anchor to which all drift is corrected towards
  # (ie. control well with the least amount of dye accumulation during imaging,
  # so is the most "normal" well)
  meta_dmso <- meta[meta$Well %in% dmso_wells, ]
  first_dmso_row <- meta_dmso %>%
    slice_min(Order, n = 1, with_ties = FALSE)
  first_dmso_id <- rownames(first_dmso_row)
  # apply correction factor from corr_table to ALL (including non DMSO) wells
  # for() loop goes through each feature
  for (i in seq_len(ncol(df))) {
    # skips correction for plate drift if intercept only is best fit
    if (fits$best_p[i] != 0) {
      # vector of predicted drift for all wells along imaging order
      pred_vals <- corr_table[, i]
      # reference value to first DMSO well (ie well with no drift)
      norm_val <- corr_table[first_dmso_id, i]
      # correction formula
      # note norm_val keeps the data anchored to first DMSO well 
      # rather than forcing curve fit to flatten at 0
      corr_df[, i] <- norm_val +
        df[, i] -
        pred_vals
    }
  }
return(list(
  # return 1. df of correction factors, 
  corr_table = corr_table,
  # return 2. df of corrected data
  corr_df = corr_df
))
}
# run function to apply plate drift correction ####
batches <- map(
  batches,
  function(batch_obj) {
    corr <- drift_corr(
      df = batch_obj$df,
      meta = batch_obj$meta,
      dmso_wells = batch_obj$dmso_wells,
      fits = batch_obj$fits
    )
    batch_obj$corr_table <- corr$corr_table
    batch_obj$corr_df <- corr$corr_df
    batch_obj
  }
)

# save corr data ####
walk(batches, function(batch_obj) {
  file_name <- batch_obj$file_name
  write.csv(
    batch_obj$fits,
    paste("outputs/data/", file_name, "_platedrift_fits.csv", sep = "")
  )
  write.csv(
    batch_obj$corr_table,
    paste("outputs/data/", file_name, "_platedrift_corrtable.csv", sep = "")
  )
  write.csv(
    batch_obj$corr_df,
    paste("data/processed/", file_name, "_platedrift_corr.csv", sep = "")
  )
})
rm(list = ls())
