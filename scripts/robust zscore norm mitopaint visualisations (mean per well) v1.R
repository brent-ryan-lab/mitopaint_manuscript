# Title: robust zscore norm mitopaint visualisations (mean per well) v1
# Step: 3.1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 31-07-26

# load packages ####
library(data.table)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
# set variables ####
batches_info <- list(
  N1 = list(
    file_name = "SF240627_mPaintDR2_N2",
    dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")
  )
)
# create a function to load data ####
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
  # load zscore data as z
  z <- as.data.frame(
    fread(
      paste("data/processed/", file_name, "_z.csv", sep = ""),
      header = TRUE
    )
  )
    # keep rownames as WELL_BATCH
    rownames(z) <- z$V1
   z$V1 <- NULL
  return(list(
    file_name = file_name,
    dmso_wells = dmso_wells,
    df = df,
    meta = meta,
    z = z
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
# subset dmso wells for all Ns in batches #
batches <- map(
  batches,
  function(batch_obj) {
    # df is subsetted to just wells present in dmso_wells
    batch_obj$df_dmso <- batch_obj$df[
      rownames(batch_obj$df) %in% rownames(batch_obj$meta[batch_obj$meta$Compound == "DMSO",]),
      ,
      drop = FALSE
    ]
    # z is subsetted to just wells present in dmso_wells
    batch_obj$z_dmso <- batch_obj$z[
      rownames(batch_obj$z)  %in% rownames(batch_obj$meta[batch_obj$meta$Compound == "DMSO",]),
      ,
      drop = FALSE
    ]
    # return batch_obj
    batch_obj
  }
)
# create function to plot density of feature values before z score correction ####
plot_density <- function(df, plot_title) {
  # df is the input data matrix of rows = wells, cols = features
  df |>
    # converts df into long format for density curve plot
    pivot_longer(
      cols = everything(),
      names_to = "Feature",
      values_to = "Value"
    ) |>
    # removes NAs and infinite values
    filter(is.finite(Value)) |>
    ggplot(
      aes(
        x = Value,
        # each line is one feature
        group = Feature
      )
    ) +
    geom_density(
      # density lines are transparent black
      colour = scales::alpha("black", 0.1),
      fill = NA,
      linewidth = 0.5
    ) +
    labs( 
      # title is given at beginning of function with plot_title
      title = plot_title,
      x = "Feature value",
      y = "Density"
    ) +
    theme_pubr()+
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    )
}
# run function to plot density of feature values before and after z score correction ####
plots <- map(
  batches,
  function(batch_obj) {
    # map() loops plot_density for df and z for any number of Ns in batches
    list(
      # first plot uses df 
      raw = plot_density(
        df = batch_obj$df_dmso,
        plot_title = paste(
          batch_obj$batch_name,
          "DMSO Distributions Before \nZ-score Scaling"
        )
      ), # close bracket for first plot
      # second plot uses z
      zscore = plot_density(
        df = batch_obj$z_dmso,
        plot_title = paste(
          batch_obj$batch_name,
          "DMSO Distributions After \nZ-score Scaling"
        )
      ) # close bracket for second plot
    ) # close bracket for list
  } # close bracket for function
) # close bracket for map() loop

# print plots into viewer (NOTE: may take a long time as 3000+ lines per plot)  
walk(
  plots,
  ~ walk(.x, print)
)
# save plots ####
iwalk(
  plots,
  function(batch_plots, batch_name) {
    # loops through each plot (before and after z-score) for each N
    # batch_file_name is the file_name for that N
    batch_file_name <- batches[[batch_name]]$file_name
    # therefore filename is something like: 
    # "SF240627_mPaintDR2_N2_zscore_density"
    # and is saved in outputs/figures/
    iwalk(
      batch_plots,
      function(plot, plot_name) {
        
        ggsave(
          filename = paste0(
            "outputs/figures/",
            batch_file_name,
            "_",
            plot_name,
            "_density.pdf"
          ),
          plot = plot,
          width = 3,
          height = 3
        )
      } # end of loop for each plot_var
    )
  } # end of loop for each N
)
rm(list = ls())
