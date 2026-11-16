# Title: robust zscore norm mitopaint vis (mean per well) v1
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
  N2 = list(
    file_name = "SF261109_mPaint2DR2_N2",
    dmso_wells = c(
      "2_15","2_8","3_8","3_15","4_15","4_8","5_8","5_15","6_8","7_8","7_15",
      "8_15","8_8","9_8","9_15","10_15","10_8","11_8","11_15","12_15","12_8",
      "13_8","13_15","14_15","14_8","15_8","15_15","16_15","16_8","17_8",
      "17_15","18_15","18_8","19_8","19_15","20_15","20_8","21_8","22_15",
      "22_8","23_8","23_15"
    )
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
batches <- purrr::map(
  batches_info,
  function(batch_info) {
    load_data(
      file_name = batch_info$file_name,
      dmso_wells = batch_info$dmso_wells
    )
  }
)
# subset dmso wells for all Ns in batches #
batches <- purrr::map(
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
# make format x labels helper function (so that it is scientific notation just for raw values)
format_x_labels <- function(x) {
  if (max(abs(x), na.rm = TRUE) < 100) {
    scales::label_number()(x)
  } else {
    scales::label_scientific()(x)
  }
}

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
    # format x labels to scientific notation if large range
    scale_x_continuous(
      labels = format_x_labels
    ) +
    labs( 
      # title is given at beginning of function with plot_title
      title = plot_title,
      x = "Feature value",
      y = "Density"
    ) +
    theme_pubr()+
    theme(
      aspect.ratio = 1,
      plot.margin = margin(t = 5, r = 5, b = 15, l = 5),
      panel.grid = element_blank(),
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}
# run function to plot density of feature values before and after z score correction ####
plots <- purrr::map(
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
# walk(plots,~ walk(.x, print))

# align plots
plots <- purrr::map(plots, ~ cowplot::align_plots(
    plotlist = .x,
    align = "v",
    axis = "l"
  )
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
          height = 3.2
        )
      } # end of loop for each plot_var
    )
  } # end of loop for each N
)
rm(list = ls())
