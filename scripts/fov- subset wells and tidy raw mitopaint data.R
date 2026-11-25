# Title: fov- subset wells and tidy raw mitopaint data v1
# Step: 1.1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 25-11-2026

# load packages ####
library(data.table)
library(tidyverse)
library(stringr)
library(purrr)
# set variables ####
batches_info <- list(
  FOV1 = list(
    file_name = "SF240627_mPaintDR2_N2_FOV1",
    batch_name = "FOV1"
  ),
  FOV2 = list(
    file_name = "SF240627_mPaintDR2_N2_FOV2",
    batch_name = "FOV2"
  ),
  FOV3 = list(
    file_name = "SF240627_mPaintDR2_N2_FOV3",
    batch_name = "FOV3"
  ),
  FOV4 = list(
    file_name = "SF240627_mPaintDR2_N2_FOV4",
    batch_name = "FOV4"
  ),
  FOV5 = list(
    file_name = "SF240627_mPaintDR2_N2_FOV5",
    batch_name = "FOV5"
  ),
  FOV6 = list(
    file_name = "SF240627_mPaintDR2_N2_FOV6",
    batch_name = "FOV6"
  ),
  FOV7 = list(
    file_name = "SF240627_mPaintDR2_N2_FOV7",
    batch_name = "FOV7"
  ),
  FOV8 = list(
    file_name = "SF240627_mPaintDR2_N2_FOV8",
    batch_name = "FOV8"
  ),
  FOV9 = list(
    file_name = "SF240627_mPaintDR2_N2_FOV9",
    batch_name = "FOV9"
  )
)
file_name <- "SF240627_mPaintDR2_N2"
rm_cols <- c("Timepoint",
             "Number of Analyzed Fields",
             "Time [s]",
             "Temperature",
             "Target Temperature",
             "CO2",	"Target CO2",
             "Nuclei - Number of Objects",
             "Non-border cells Selected - Number of Objects",
             "Non-border cells Selected - Nucleus Area [µm²] - Mean per Well",
             "Non-border cells Selected - Nucleus Roundness - Mean per Well",
             "Non-border cells Selected - Cell Area [µm²] - Mean per Well",
             "Non-border cells Selected - Cell Roundness - Mean per Well",
             "Non-border cells Selected - Intensity Cytoplasm mKeima ph7 Mean - Mean per Well",
             "Non-border cells Selected - Intensity Cytoplasm mKeima ph4/TMRM Mean - Mean per Well",
             "Non-border cells Selected - mKeima ph4/ph7 ratio - Mean per Well",
             "Non-border cells Selected - mkeima ph7 mitochondria Area [µm²] - Mean per Well",
             "Non-border cells Selected - mkeima ph7 mitochondria Roundness - Mean per Well",
             "Non-border cells Selected - mkeima ph7 mitochondria Width [µm] - Mean per Well",
             "Non-border cells Selected - mkeima ph7 mitochondria Length [µm] - Mean per Well",
             "Non-border cells Selected - mkeima ph7 mitochondria Ratio Width to Length - Mean per Well",
             "Cell Type",	
             "Cell Count"
)
nuc_count <- "Non-border cells Selected - Number of Objects"
meta_cols <- c("Row",
               "Column",
               "Compound",	
               "Concentration",
               "Condition")
rm_cond <- NULL
keep_cond <- NULL
min_nuc <- 100
wells_info <- list(
  WELL1 = c(2,23),
  WELL2 = c(2,9,16,23),
  WELL3 = c(2,6,10,14,18,23),
  WELL4 = c(1,5,8,11,14,17,20,23),
  WELL6 = c(2,4,6,8,10,12,13,15,17,19,21,23),
  WELL11 = NULL
)
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
    filter(.data[[nuc_count]] > min_nuc)
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
    str_remove("Cell Painting ") |>
    str_replace("mKeima ph4(?:/TMRM)?", "mt-Keima pH4") |>
    str_replace("mKeima ph7", "mt-Keima pH7") |>
    str_trim(side = "right")
  # convert all features to numeric
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  # remove any columns with NA/ non finite values
  df <- df[, colSums(!is.finite(as.matrix(df))) == 0, drop = FALSE]
  # remove any rows with NA/ non finite values
  df <- df[apply(df, 1, function(x) all(is.finite(x))), , drop = FALSE]
  # keep meta and df aligned
  meta <- meta[rownames(df), , drop = FALSE]
  # remove any rows with Condition in rm_cond
  if (!is.null(rm_cond)) {
    keep_rows <- !meta$Condition %in% rm_cond
    meta <- meta[keep_rows, , drop = FALSE]
    df <- df[rownames(meta), , drop = FALSE]
  }
  # only keep an rows with Condition in keep_cond
  if (!is.null(keep_cond)) {
    keep_rows <- meta$Condition %in% keep_cond
    meta <- meta[keep_rows, , drop = FALSE]
    df <- df[rownames(meta), , drop = FALSE]
  }
  # return data and metadata as a list
  return(list(data = df,
              meta = meta,
              file_name = file_name,
              batch_name = batch_name)
  )
}
# run function to load and tidy data ####
batches <- purrr::map(
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
# create function to subset wells in data ####
subset_wells <- function(batch_obj, wells_info) {
  # pull data and metadata
  df <- batch_obj$data
  meta <- batch_obj$meta
  # loop through each WELL group defined in wells_info
  well_subsets <- purrr::imap(
    wells_info,
    function(well_cols, well_name) {
      # if NULL, retain all columns
      if (is.null(well_cols)) {
        meta_well <- meta
      } else {
        # subset metadata by selected plate columns
        meta_well <- meta[
          meta$Column %in% well_cols,
          ,
          drop = FALSE
        ]
      }
      # subset data using matching rownames
      # explicitly preserve the order of meta_well
      data_well <- df[
        rownames(meta_well),
        ,
        drop = FALSE
      ]
      # return one WELL sublist
      list(
        meta = meta_well,
        data = data_well
      )
    }
  )
  return(well_subsets)
}
# run function to subset wells for every FOV ####
batches_wells <- purrr::map(
  batches,
  function(batch_obj) {
    subset_wells(
      batch_obj = batch_obj,
      wells_info = wells_info
    )
  }
)
# save tidy data for each FOV and WELL subset ####
# loop through each FOV
purrr::iwalk(
  batches_wells,
  function(fov_obj, fov_name) {
    # loop through each WELL subset
    purrr::iwalk(
      fov_obj,
      function(well_obj, well_name) {
        # unique filename for this FOV and WELL
        output_name <- paste0(
          file_name,
          "_",
          fov_name,
          "_",
          well_name
        )
        # save data
        write.csv(
          well_obj$data,
          paste0(
            "data/processed/",
            output_name,
            "_data_tidy.csv"
          )
        )
        # save metadata
        write.csv(
          well_obj$meta,
          paste0(
            "data/processed/",
            output_name,
            "_meta_tidy.csv"
          )
        )
      }
    )
  }
)