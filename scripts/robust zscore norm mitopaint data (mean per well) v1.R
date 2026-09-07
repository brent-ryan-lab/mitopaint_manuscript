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
    file_name = "SF250813_mPaintFDA_1",
    dmso_wells = c("2_2","4_2","6_2","8_2","10_2","12_2",
                   "14_2","16_2","18_2","20_2","22_2",
                   "3_9","5_9","7_9","9_9","11_9","13_9",
                   "15_9","17_9","19_9","21_9","23_9")
  ),
  N2 = list(
    file_name = "SF250813_mPaintFDA_2",
    dmso_wells = c("2_8","4_8","6_8","8_8","10_8","12_8",
                   "14_8","16_8","18_8","20_8","22_8",
                   "3_15","5_15","7_15","9_15","11_15",
                   "13_15","15_15","17_15","19_15","21_15","23_15")
  ),
  N3 = list(
    file_name = "SF250813_mPaintFDA_3",
    dmso_wells = c("2_2","4_2","6_2","8_2","10_2","12_2",
                   "14_2","16_2","18_2","20_2","22_2",
                   "3_9","5_9","7_9","9_9","11_9","13_9",
                   "15_9","17_9","19_9","21_9","23_9")
  ),
  N4 = list(
    file_name = "SF250813_mPaintFDA_4",
    dmso_wells = c("2_2","4_2","6_2","8_2","10_2","12_2",
                   "14_2","16_2","18_2","20_2","22_2",
                   "3_9","5_9","7_9","9_9","11_9","13_9",
                   "15_9","17_9","19_9","21_9","23_9")
  ),
  N5 = list(
    file_name = "SF250813_mPaintFDA_5",
    dmso_wells = c("2_2","4_2","6_2","8_2","10_2","12_2",
                   "14_2","16_2","18_2","20_2","22_2",
                   "3_9","5_9","7_9","9_9","11_9","13_9",
                   "15_9","17_9","19_9","21_9","23_9")
  ),
  N6 = list(
    file_name = "SF250813_mPaintFDA_6",
    dmso_wells = c("2_2","4_2","6_2","8_2","10_2","12_2",
                   "14_2","16_2","18_2","20_2","22_2",
                   "3_9","5_9","7_9","9_9","11_9","13_9",
                   "15_9","17_9","19_9","21_9","23_9")
  ),
  N7 = list(
    file_name = "SF250813_mPaintFDA_7",
    dmso_wells = c("2_2","4_2","6_2","8_2","10_2","12_2",
                   "14_2","16_2","18_2","20_2","22_2",
                   "3_9","5_9","7_9","9_9","11_9","13_9",
                   "15_9","17_9","19_9","21_9","23_9")
  ),
  N8 = list(
    file_name = "SF250813_mPaintFDA_8",
    dmso_wells = c("2_2","4_2","6_2","8_2","10_2","12_2",
                   "14_2","16_2","18_2","20_2","22_2",
                   "3_9","5_9","7_9","9_9","11_9","13_9",
                   "15_9","17_9","19_9","21_9","23_9")
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
