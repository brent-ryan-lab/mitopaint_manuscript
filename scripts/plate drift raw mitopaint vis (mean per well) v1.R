# Title: plate drift raw mitopaint vis (mean per well) v1
# Step: 2.2
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 06-08-26

# load packages ####
library(data.table)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(ggplot2)
library(colorspace)
library(viridis)
library(cowplot)
# set file variables ####
batches_info <- list(
  N1 = list(
    file_name = "SF240215_mPaintDrift_DMSO_N1",
    batch_name = "N1"
  )
)
meta_cols <- c("Row",
                "Column",
                "Compound",	
                "Concentration")
rm_cols <- c("Timepoint",
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
            "Non-border cells Selected - mKeima ph4/ph7 ratio - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Area [¬µm¬≤] - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Roundness - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Width [¬µm] - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Length [¬µm] - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Ratio Width to Length - Mean per Well",
            "Cell Type",	
            "Cell Count"
)
plot_var <- c("Intensity Cytoplasm Region TMRM test Mean",
              "Intensity Cytoplasm Region CellRox Deep Red test Mean")
legend_titles <- c(
  "TMRM Intensity",
  "CellROX Intensity"
)
# create a function to load data ####
load_data <- function(file_name, batch_name, rm_cols, meta_cols) {
  # load df 
  df <- as.data.frame(
    fread(
      paste(
        "data/raw/", file_name, ".csv", sep = ""), 
      skip = "Row", header = TRUE)
  )
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
# run function to load data
batches <- map(
  # map loops over all elements in batches (if there is more than one N)
  batches_info,
  function(batch_info) {
    load_data(
      file_name = batch_info$file_name,
      batch_name = batch_info$batch_name,
      rm_cols = rm_cols,
      meta_cols = meta_cols
    )
  }
)
# open data to inspect
View(batches$N1$data)
# open meta to inspect
View(batches$N1$meta)
# create a function to plot plate overview ####
plot_plate <- function(feature, legend_title, batch_obj) {
  # load data and metadata
  meta <- batch_obj$meta
  df <- batch_obj$data
  # plotting data subsets columns to just the row, column and plot_var
  plot_df <- meta |>
    mutate(
      Row_num = as.integer(Row),
      Col_num = as.integer(Column),
      value = df[[feature]]
    )
  # plotting data made into a grid plate map of 384w plate (including edge wells)
  plate_df <- expand.grid(
    Row_num = 1:16,
    Col_num = 1:24
  ) |>
    as_tibble() |>
    # left_join connects row and column to plot_var
    # tibble now includes rows for edge wells (which are NA if only inner 308w were seeded)
    left_join(
      plot_df |>
        select(Row_num, Col_num, value),
      by = c("Row_num", "Col_num")
    ) |>
   mutate(
      x = Col_num,
      # rows are ordered top to bottom
      y = 17 - Row_num
    )
  # ggplot plate plot overview
  ggplot(plate_df, aes(x = x, y = y)) +
    # each well is a circle
    geom_point(
      aes(fill = value),
      shape = 21,
      size = 6,
      colour = "black",
      stroke = 0.4
    ) +
    # fill is viridis
    scale_fill_gradientn(
      name = legend_title,
      colours = rev(
        lighten(
          viridisLite::viridis(100),
          amount = 0.3
        )
      ),
      # NA (edge wells) are white
      na.value = "white"
    ) +
    scale_x_continuous(
      breaks = 1:24,
      labels = 1:24,
      position = "top",
      expand = c(0.02, 0.02)
    ) +
    scale_y_continuous(
      # rows are ordered top to bottom
      breaks = 1:16,
      labels = 16:1,
      expand = c(0.02, 0.02)
    ) +
    coord_equal() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10)
    )
}
# run function to plot plate overview ####
plots <- map(
  # loop through each N in batches
  batches,
  function(batch_obj) {
    # loop through each feature in plot_var (full name) and legend_titles (shorthand)
    # length of plot_var and legend_titles must match, and be in matching order
    map2(
      plot_var,
      legend_titles,
      function(feature, legend_title) {
        # run plot_plate function
        plot_plate(
          feature = feature,
          legend_title = legend_title,
          batch_obj = batch_obj
        )
      }
    ) |>
      # set plot names to match legend_titles
      set_names(legend_titles)
  }
)

# create function to add fixed legend space ####
add_fixed_legend_space <- function(plot,
                                   plot_width = 1.2,
                                   legend_width = 0.2) {
  legend <- get_legend(
    plot +
      theme(
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0)
      )
  )
  plot_without_legend <- plot +
    theme(
      legend.position = "none",
      plot.margin = margin(5, 0, 5, 5)
    )
  plot_grid(
    plot_without_legend,
    legend,
    nrow = 1,
    rel_widths = c(
      plot_width,
      legend_width
    ),
    align = "h",
    axis = "tb"
  )
}
# apply consistent legend space to plots ####
plots_fixed <- map(
  plots,
  function(batch_plots) {
    map(
      batch_plots,
      add_fixed_legend_space
    )
  }
)
# print plots (walk through each N and each plot_var)
walk(
  plots_fixed,
  ~ walk(.x, print)
)
# save plots ####
iwalk(
  # loops through each N
  plots_fixed,
  function(batch_plots, batch_name) {
    # loops through each plot_var
    iwalk(
      batch_plots,
      function(plot, plot_name) {
        # safe_name is the legend_titles with spaces replaced with _
        safe_name <- str_replace_all(
          plot_name,
          "[^[:alnum:]]+",
          "_"
        )
        # batch_file_name is the file_name for that N
        batch_file_name <- batches[[batch_name]]$file_name
        # therefore filename is something like: 
        # "SF240215_mPaintDrift_DMSO_N1_CellROX_Intensity_raw_drift"
        # and is saved in outputs/figures/
        ggsave(
          filename = paste0(
            "outputs/figures/",
            batch_file_name,
            "_",
            safe_name,
            "_raw_drift.pdf"),
          plot = plot,
          # plot size manually set here
          width = 9,
          height = 6)
      }
    ) # end of loop for each plot_var
  }
) # end of loop for each N
rm(list = ls())
