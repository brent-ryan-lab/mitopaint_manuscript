# Title: robust zscore norm mitopaint data (mean per well) v1
# Step: 3
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 18-06-2026

# load packages ####
library(data.table)
library(purrr)
# set variables ####
batches_info <- list(
  N1 = list(
    file_name = "SF240627_mPaintDR2_N2",
    dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")
  ),
  N2 = list(
    file_name = "SF240704_mPaintDR2_N3",
    dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")
  ),
  N3 = list(
    file_name = "SF240711_mPaintDR2_N4",
    dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")
  )
)
# create function to load data ####
load_data <- function(file_name, dmso_wells) {
  # load plate drift corrected data as df
  df <- as.data.frame(
    fread(
      paste("data/processed/", file_name, "_platedrift_corr.csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames as WELL_BATCH
  rownames(df) <- df$V1
  df$V1 <- NULL
  # load metadata as meta
  meta <- as.data.frame(
    fread(
      paste("data/processed/", file_name, "_meta_tidy.csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames as WELL_BATCH
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
# create function to calculate zscore relative to dmso wells ####
calc_zscore <- function(df, meta, dmso_wells) {
  # df_z is where zscore values will be stored
  # to start, it is just populated with df data but will soon be replaced
  df_z <- df
  # subset dmso wells
  meta_dmso <- meta[meta$Well %in% dmso_wells, , drop = FALSE]
  df_dmso <- df[rownames(meta_dmso), , drop = FALSE]
  # for() loop calculates z score across each feature col
  for(col in colnames(df)) {
    # calculate median of DMSO wells
    med <- median(df_dmso[[col]], na.rm = TRUE)
    # calculate mean absolute deviation of DMSO wells
    mad <- mad(df_dmso[[col]], na.rm = TRUE)
    # z = (feature value - med)/mad
    # if statement is to avoid divide by 0 problems, if MAD = 0 in DMSO wells
    if (mad == 0 || is.na(mad)) {
      df_z[[col]] <- NA_real_
    } else {
      df_z[[col]] <- (df[[col]] - med) / mad
    }
  }
  return(df_z)
}
# run function to calculate zscore relative to dmso wells ####
batches <- map(
  batches,
  function(batch_obj) {
    batch_obj$df_z <- calc_zscore(
      df = batch_obj$df,
      meta = batch_obj$meta,
      dmso_wells = batch_obj$dmso_wells
    )
    batch_obj
  }
)
# save zscore data ####
walk(
  batches,
  function(batch_obj) {
    write.csv(
      batch_obj$df_z,
      paste("data/processed/", batch_obj$file_name, "_z.csv", sep = "")
    )
  }
)
rm(list = ls())
