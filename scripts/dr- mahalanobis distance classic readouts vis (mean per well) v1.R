# Title: # dr- mahalanobis distance classic readouts vis (mean per well) v1
# Step: 9.3
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 07-08-26

# load packages ####
library(data.table)
library(ggplot2)
library(ggpubr)
library(colorspace)
# set file variables ####
file_name <- "mPaintDR2_N2_N3_N4"
md_file_name <- "mPaintDR2_N2_N3_N4_mahal_PC10"
classic_file_name <- "mPaint_DR2_Classic_N2_N3_N4"
integrate_state <- "integrated"
redu_state <- "redu"
pc_use <- 10
plot_feats <- c("Intensity Cytoplasm CellRox Mean",
                "Intensity Cytoplasm TMRM Mean",
                "Mitochondria Selected Ratio Width to Length",
                "Number of Selected Spots/ Selected Cell")
y_lab <- c("Cytoplasm ROS Intensity (Z-score)",
           "Cytoplasm MMP Intensity (Z-score)",
           "Mitochondria Width:Length (Z-score)",
           "Mitophagy Spots (Z-score)")
plot_cond <- c("CCCP", "ROT")
pastel_cols <- lighten(c("#238A8DFF", "#FDE725FF"), amount = 0.3)
# create function to load data ####
y_lab_lookup <- setNames(
  y_lab,
  plot_feats
)
load_data <- function(file_name,
                      md_file_name,
                      classic_file_name,
                      integrate_state) {
  # load md
  md <- as.data.frame(
    fread(
      paste("data/processed/", md_file_name, ".csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames as WELL_BATCH
  rownames(md) <- md$V1
  md$V1 <- NULL
  # load classic readout
  classic <- as.data.frame(
    fread(
      paste("data/processed/", classic_file_name, ".csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames as WELL_BATCH
  rownames(classic) <- classic$V1
  classic$V1 <- NULL
  # load metadata
  meta <- as.data.frame(
    fread(
      paste("data/processed/", file_name, "_", "meta_", 
            integrate_state, ".csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames as WELL_BATCH
  rownames(meta) <- meta$V1
  meta$V1 <- NULL
  # remove meta cols from classic
  classic <- classic[,!colnames(classic) %in% colnames(meta)]
  # return list
  return(list(
    md = md,
    classic = classic,
    meta = meta
  ))
}
# load data ####
data <- load_data(file_name,
                  md_file_name,
                  classic_file_name,
                  integrate_state)
# create function to calculate classic z score to dmso by batch ####
zscore <- function(data) {
  df <- data$classic
  meta <- data$meta
  df_dmso <- df[rownames(df) %in% rownames(meta[meta$Compound == "DMSO",]),]

  
  
  
}