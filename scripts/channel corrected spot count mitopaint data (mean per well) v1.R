# Title: channel corrected spot count mitopaint data (mean per well) v1
# Step: 9.2
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 06-08-26

# load packages ####
library(data.table)
library(colorspace)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(cowplot)
library(purrr)
library(grid)
library(stringr)
# set file variables ####
batches_info <- list(
  N1 = list(
    file_name = "SF240627_mPaintDR2_N2_SpotComparison",
    dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9"),
    batch_name = "N1"
  )
)
nuc_count = "Non-border cells Selected - Number of Objects"
meta_cols = c("Row",
              "Column",
              "Compound",	
              "Concentration")
file_name <- "mPaint_DR2_SpotComparison"
# set plot variables ####
# make sure order of plot_cond and plot_lab are the same
# make sure order of plot_feats and y_lab are the same
plot_cond <- c("DMSO_0", "CCCP_30")
plot_lab <- c("DMSO", "CCCP")
plot_feats <- c("Number of Selected Spots/ Selected Cell",
                "Number of Raw pH4 Spots/ Selected Cell"
                )
y_lab <- c("Selected Spots/ Cell",
           "Raw pH4 Spots/ Cell"
           )
x_lab <- "Compound"
pastel_cols <- lighten(c("#440154FF", "#238A8DFF"), amount = 0.3)
point_size <- 3
size_axis <- 10
# pair up plot_feats and y_lab ####
# pairs link the shorthand (y_lab) for each feature to the full feature name
# make sure input order matches in both
y_lab_lookup <- setNames(
  y_lab,
  plot_feats
)
# create function to load data ####
load_data <- function(file_name, batch_name, dmso_wells) {
  # load classic data as df
  df <- as.data.frame(
    fread(
      paste("data/raw/", file_name, ".csv", sep = ""),
      skip = "Row", header = TRUE
    )
  )
  # keep rownames as WELL_BATCH
  rownames(df) <- df$V1
  df$V1 <- NULL
  # separate metadata 
  meta <- df %>%
    select(any_of(meta_cols))
  df <- df %>%
    select(-any_of(meta_cols))
  # populate additional metadata
  meta$Well  <- paste(meta$Column, meta$Row, sep = "_")
  meta$Batch <- batch_name
  meta$ID    <- paste(meta$Well, meta$Batch, sep = "_")
  meta$Condition <- paste(meta$Compound, meta$Concentration, sep = "_")
  rownames(meta) <- meta$ID
  rownames(df)   <- meta$ID
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
              batch_name = batch_name,
              dmso_wells = dmso_wells)
  )
}
# run function to load data ####
# loop through loading all Ns given in batch_info
batches <- map(
  batches_info,
  function(batch_info) {
    load_data(
      file_name = batch_info$file_name,
      batch_name = batch_info$batch_name,
      dmso_wells = batch_info$dmso_wells
    )
  }
)
# create function to calculate z prime ####
calculate_zprime <- function(obj,
                             pos_ctrl,
                             neg_ctrl) {
  df <- obj$data
  meta <- obj$meta
  # ensure data and metadata have matching row order
  meta <- meta[rownames(df), , drop = FALSE]
  # combine feature data and metadata
  zprime_df <- cbind(df, meta) |>
    filter(
      Condition %in% c(
        neg_ctrl,
        pos_ctrl
      )
    )
  # convert all feature columns into long format
  zprime_long <- zprime_df |>
    pivot_longer(
      cols = all_of(colnames(df)),
      names_to = "Feature",
      values_to = "Value"
    )
  # calculate control means, SDs and well counts
  control_summary <- zprime_long |>
    group_by(
      Feature,
      Condition
    ) |>
    summarise(
      mean_value = mean(Value, na.rm = TRUE),
      sd_value = sd(Value, na.rm = TRUE),
      n_wells = sum(!is.na(Value)),
      .groups = "drop"
    ) |>
    pivot_wider(
      names_from = Condition,
      values_from = c(
        mean_value,
        sd_value,
        n_wells
      ),
      names_glue = "{.value}_{Condition}"
    )
  # dynamically construct control column names
  pos_mean_col <- paste0("mean_value_", pos_ctrl)
  neg_mean_col <- paste0("mean_value_", neg_ctrl)
  pos_sd_col <- paste0("sd_value_", pos_ctrl)
  neg_sd_col <- paste0("sd_value_", neg_ctrl)
  # calculate Z-prime
  control_summary |>
    mutate(
      Batch = obj$batch_name,
      file_name = obj$file_name,
      mean_difference = abs(
        .data[[pos_mean_col]] -
          .data[[neg_mean_col]]
      ),
      z_prime = 1 -
        (
          3 * (
            .data[[pos_sd_col]] +
              .data[[neg_sd_col]]
          )
        ) /
        mean_difference
    ) |>
    relocate(
      Batch,
      file_name,
      Feature,
      z_prime,
      mean_difference
    )
}
# run function to calculate z prime ####
zprime_results <- map_dfr(
  batches,
  function(batch_obj) {
    calculate_zprime(
      obj = batch_obj,
      pos_ctrl = plot_cond[2],
      neg_ctrl = plot_cond[1]
    )
  }
)
# create a function to make barplot of plot_feat ####
plot_batch_feature <- function(batch_obj,
                               feature,
                               zprime_results,
                               plot_cond,
                               plot_lab,
                               y_lab_lookup,
                               x_lab = "Compound",
                               pastel_cols,
                               point_size = 3,
                               size_axis = 12) {
  # keep metadata and data aligned
  meta <- batch_obj$meta
  df <- batch_obj$data
  meta <- meta[rownames(df), , drop = FALSE]
  # combine meta with only the df col of the selected feature
  plot_df <- cbind(
    meta,
    value = df[[feature]]
  )
  # keep only the plot_cond conditions
  plot_df <- plot_df[
    plot_df$Condition %in% plot_cond,
    ,
    drop = FALSE
  ]
  # force conditions into the order given in plot_cond
  plot_df$Condition <- factor(
    plot_df$Condition,
    levels = plot_cond
  )
  # retrieve the Z-prime value for this batch and feature
  zprime_value <- zprime_results |>
    filter(
      Batch == batch_obj$batch_name,
      Feature == feature
    ) |>
    pull(z_prime)
  # use shorthand feature label where available
  feature_label <- unname(y_lab_lookup[feature])
  # title is feature label and Z-prime value
  title_text <- paste0(
    feature_label,
    "\nZ' = ",
    ifelse(
      is.na(zprime_value),
      "NA",
      sprintf("%.2f", zprime_value)
    )
  )
  # calculate upper y-axis limit with room above error bars
  summary_df <- plot_df |>
    group_by(Condition) |>
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      sem_value = sd(value, na.rm = TRUE) /
        sqrt(sum(!is.na(value))),
      .groups = "drop"
    )
  y_max <- max(
    c(
      plot_df$value,
      summary_df$mean_value + summary_df$sem_value
    ),
    na.rm = TRUE
  )
  # make plot
  ggplot(
    plot_df,
    aes(
      x = Condition,
      y = value,
      fill = Condition
    )
  ) +
    # bar height is the mean across replicate wells
    stat_summary(
      fun = mean,
      geom = "bar",
      colour = "black",
      width = 0.8
    ) +
    # error bar is SEM across replicate wells
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.5
    ) +
    # points are individual replicate wells
    geom_point(
      position = position_jitter(
        width = 0.2,
        height = 0,
        seed = 42
      ),
      shape = 16,
      colour = "black",
      size = point_size,
      alpha = 0.1
    ) +
    # fill is pastel_cols
    scale_fill_manual(
      values = setNames(
        pastel_cols[seq_along(plot_cond)],
        plot_cond
      )
    ) +
    scale_x_discrete(
      labels = setNames(
        plot_lab,
        plot_cond
      ),
      drop = FALSE
    ) +
    # y limits is 1.15*y_max to avoid cutting off error bars
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      limits = c(0, y_max * 1.15)
    ) +
    # tidy theme
    labs(
      title = title_text,
      y = feature_label,
      x = x_lab
    ) +
    theme_pubr() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = size_axis,
        face = "bold"
      ),
      axis.text = element_text(size = size_axis),
      axis.title = element_text(size = size_axis),
      panel.grid = element_blank(),
      legend.position = "none"
    )
}
# run function ####
# make one plot per feature per batch #
plots <- map(
  batches,
  function(batch_obj) {
    map(
      plot_feats,
      function(feature) {
        plot_batch_feature(
          batch_obj = batch_obj,
          feature = feature,
          zprime_results = zprime_results,
          plot_cond = plot_cond,
          plot_lab = plot_lab,
          y_lab_lookup = y_lab_lookup,
          x_lab = x_lab,
          pastel_cols = pastel_cols,
          point_size = point_size,
          size_axis = size_axis
        )
      }
    ) |>
      set_names(plot_feats)
  }
)
# give every plot the same y-axis spacing ####
# flatten nested plot list:
# plots$N1$feature1, plots$N1$feature2, plots$N2$feature1, etc.
flat_plots <- unlist(
  plots,
  recursive = FALSE
)
# align the left axes across every plot
flat_plots_aligned <- align_plots(
  plotlist = flat_plots,
  align = "v",
  axis = "l"
)
# restore plot names
names(flat_plots_aligned) <- names(flat_plots)
flat_plots_aligned
walk(flat_plots_aligned, ~ {grid.newpage(); grid.draw(.x) })
# save plots ####
iwalk(
  flat_plots_aligned,
  function(plot, plot_name) {
    safe_name <- str_replace_all(plot_name, "[^[:alnum:]]+", "_")
    ggsave(
      filename = paste0(
        "outputs/figures/",
        file_name,
        "_",
        safe_name,
        "_raw_drift.pdf"
      ),
      plot = plot,
      width = 2.2,
      height = 3
    )
  }
)
rm(list = ls())