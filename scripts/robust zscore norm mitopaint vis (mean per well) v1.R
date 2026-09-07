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
