# Title: robust zscore norm mitopaint data (mean per well) v1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 12-06-2026

# load packages ####
library(data.table)
# set variables ####
file_name <- NULL
dmso_wells <- NULL
# load data ####
# load plate drift corrected data as df
df <- as.data.frame(
  fread(
    paste(
      "data/processed/", file_name, "_platedrift_corr.csv", sep = ""), 
    header = TRUE)
)
# keep rownames as WELL_BATCH
rownames(df) <- df$V1
df$V1 <- NULL
# load metadata as meta
meta <- as.data.frame(
  fread(
    paste(
      "data/processed/", file_name, "_meta_tidy.csv", sep = ""), 
    header = TRUE)
)
# keep rownames as WELL_BATCH
rownames(meta) <- meta$V1
meta$V1 <- NULL
# create function to calculate zscore relative to dmso wells ####
calc_zscore <- function(df, meta, dmso_wells) {
  # df_z is where zscore values will be stored
  # to start, it is just populated with df data but will soon be replaced
  df_z <- df
  # subset dmso wells
  meta_dmso <- meta[meta$Well %in% dmso_wells,]
  df_dmso <- df[rownames(df) %in% rownames(meta_dmso),]
  # for() loop calculates z score across each feature col
  for(col in colnames(df)) {
    # calculate median of DMSO wells
    med <- median(df_dmso[[col]], na.rm = TRUE)
    # calculate mean absolute deviation of DMSO wells
    mad <- mad(df_dmso[[col]], na.rm = TRUE)
    # z = (feature value - med)/mad
    df_z[[col]] <- (df[[col]] - med) / mad
  }
  return(df_z)
}
# run function to calculate zscore relative to dmso wells ####
df_z <- calc_zscore(df, meta, dmso_wells)
# save zscore data ####
write.csv(df_z,
          paste(
            "data/processed/", file_name, "_z.csv", sep = "")
)
rm(list = ls())