# Title: library meta- merge library and hits metadata (mean per well) v1
# Step: 7.3
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 09-11-26

# load packages ####
library(data.table)
library(tidyverse)
# set file variables ####
meta_name <- "mPaintFDA_N1_N2_N3_N4_N5_N6_N7_N8_mahal_summary_PC10"
library_name <- "SF250909_mPaintFDA_Library"
file_name <- "mPaintFDA_N1_N2_N3_N4_N5_N6_N7_N8_meta_summary"
# load data ####
meta <- as.data.frame(
  fread(
    paste("data/processed/", meta_name, ".csv", sep = ""),
    header = TRUE
  )
)
rownames(meta) <- meta$V1
meta$V1 <- NULL
library <- as.data.frame(
  fread(
    paste("data/raw/", library_name, ".csv", sep = ""),
    header = TRUE
  )
)
rownames(library) <- library$V1
library$V1 <- NULL
# join meta and library ####
df <- left_join(meta, library, by = "Condition") |>
  arrange(desc(md2))
# save df ####
write.csv(
  df,
  paste("data/processed/", file_name,".csv", sep = "")
)
rm(list = ls())
