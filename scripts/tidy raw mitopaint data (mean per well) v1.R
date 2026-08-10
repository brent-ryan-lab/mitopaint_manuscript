# Title: tidy raw mitopaint data (mean per well) v1
# Step: 1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 18-06-2026

# load packages ####
library(data.table)
library(tidyverse)
library(stringr)
library(purrr)
# set variables ####
batches_info <- list(
  N1A = list(
    file_name = "SF260604_mPaintSpace2_N1A",
    batch_name = "N1A"
  ),
  N1B = list(
    file_name = "SF260604_mPaintSpace2_N1B",
    batch_name = "N1B"
  ),
  N2A = list(
    file_name = "SF260604_mPaintSpace2_N2A",
    batch_name = "N2A"
  ),
  N2B = list(
    file_name = "SF260604_mPaintSpace2_N2B",
    batch_name = "N2B"
  ),
  N3A = list(
    file_name = "SF260701_mPaintSpace2_N3A",
    batch_name = "N3A"
  ),
  N3B = list(
    file_name = "SF260701_mPaintSpace2_N3B",
    batch_name = "N3B"
  )
)
rm_cols = c("Timepoint",
            "Number of Analyzed Fields",
            "Time [s]",
            "Temperature",
            "Target Temperature",
            "CO2",	"Target CO2",
            "Nuclei - Number of Objects",
            "Non-border cells Selected - Number of Objects",
            "Non-border cells Selected - Nucleus Area [¬µm¬≤] - Mean per Well",
            "Non-border cells Selected - Nucleus Roundness - Mean per Well",
            "Non-border cells Selected - Cell Area [¬µm¬≤] - Mean per Well",
            "Non-border cells Selected - Cell Roundness - Mean per Well",
            "Non-border cells Selected - Intensity Cytoplasm mKeima ph7 Mean - Mean per Well",
            "Non-border cells Selected - Intensity Cytoplasm mKeima ph4/TMRM Mean - Mean per Well",
            "Non-border cells Selected - Intensity Cytoplasm TMRM test Mean - Mean per Well",
            "Non-border cells Selected - mKeima ph4/ph7 ratio - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Area [¬µm¬≤] - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Roundness - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Width [¬µm] - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Length [¬µm] - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Ratio Width to Length - Mean per Well",
            "Cell Type",	
            "Cell Count"
)
nuc_count = "Non-border cells Selected - Number of Objects"
meta_cols = c("Row",
              "Column",
              "Compound",	
              "Concentration",
              "Condition")
# create function to load and tidy data ####
load_data <- function(file_name, batch_name, rm_cols, meta_cols, nuc_count) {
  # load df 
  df <- as.data.frame(
    fread(
      paste(
        "data/raw/", file_name, ".csv", sep = ""), 
      skip = "Row", header = TRUE)
    )
  # remove any rows with low nuclei
  df <- df %>% 
    filter(.data[[nuc_count]] > 100)
  # remove unwanted columns
  df <- df %>%
    select(-any_of(rm_cols))
  # separate metadata 
  meta <- df %>%
    select(any_of(meta_cols))
  df <- df %>%
    select(-any_of(meta_cols))
  # populate additional metadata
  meta$Well  <- paste(meta$Column, meta$Row, sep = "_")
  meta$Batch <- batch_name
  meta$ID    <- paste(meta$Well, meta$Batch, sep = "_")
  meta$Order <- c(1:(nrow(meta)))
  meta$Condition <- paste(meta$Compound, meta$Concentration, sep = "_")
  rownames(meta) <- meta$ID
  rownames(df)   <- meta$ID
  # clean column names
  names(df) <- names(df) |>
    str_remove("^Non-border cells Selected - ") |>
    str_remove(" - Mean per Well$") |>
    str_replace("mKeima ph4/TMRM", "mt-Keima pH4") |>
    str_replace("mKeima ph7", "mt-Keima pH7") |>
    str_trim(side = "right")
  # convert all features to numeric
  df[] <- lapply(df, as.numeric)
  # remove any columns with NA
  df <- df[, colSums(is.na(df)) == 0, drop = FALSE]
  # remove any rows with NA
  df <- df[rowSums(is.na(df)) == 0, , drop = FALSE]
  # keep meta and df aligned
  meta <- meta[rownames(df), , drop = FALSE]
  # return data and metadata as a list
  return(list(data = df,
              meta = meta,
              file_name = file_name,
              batch_name = batch_name)
         )
}
# run function to load and tidy data ####
batches <- map(
  batches_info,
  function(batch_info) {
    load_data(
      file_name = batch_info$file_name,
      batch_name = batch_info$batch_name,
      rm_cols = rm_cols,
      meta_cols = meta_cols,
      nuc_count = nuc_count
    )
  }
)
# save tidy data ####
walk(batches, function(batch_obj) {
  write.csv(
    batch_obj$data,
    paste("data/processed/", batch_obj$file_name, "_data_tidy.csv", sep = "")
  )
  write.csv(
    batch_obj$meta,
    paste("data/processed/", batch_obj$file_name, "_meta_tidy.csv", sep = "")
  )
})
rm(list = ls())
