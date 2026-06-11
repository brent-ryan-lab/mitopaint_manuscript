# Title: tidy raw mitopaint data (mean per well) v1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 11-06-2026

# set variables ####
file_name <- "SF240627_mPaintDR2_N2"
batch_name <- "N2"
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
               "Concentration"
               )
# load packages ####
library(data.table)
library(tidyverse)
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
    filter(nuc_count > 100)
  
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
    str_remove(" - Mean per Well$")
  
  # convert all features to numeric
  df[] <- lapply(df, as.numeric)
  
  # remove any columns with NA
  df <- df[, colSums(is.na(df)) == 0]
  
  # remove any rows with NA
  df <- df[rowSums(is.na(df)) == 0, ]
  
  # keep meta and df aligned
  meta <- meta[rownames(df), ]
  
  # return data and metadata as a list
  return(list(data = df, meta = meta))
}
# run function to load and tidy data ####
df_tidy <- load_data(file_name, batch_name, rm_cols, meta_cols, nuc_count)
# open data to inspect
View(df_tidy[["data"]])
# open meta to inspect
View(df_tidy[["meta"]])
# save tidy data ####
write.csv(df_tidy$data,
          paste(
            "data/processed/", file_name, "_data_tidy.csv", sep = "")
          )

write.csv(df_tidy$meta,
          paste(
            "data/processed/", file_name, "_meta_tidy.csv", sep = "")
)

