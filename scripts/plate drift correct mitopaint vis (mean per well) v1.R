# Title: plate drift correct mitopaint vis (mean per well) v1
# Step: 2.1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 18-06-2026

# load packages ####
library(data.table)
library(tidyverse)
library(colorspace)
library(ggplot2)
library(cowplot)
library(stringr)
library(ggpubr)
library(purrr)
# set variables ####
batches_info <- list(
  N1A = list(
    file_name = "SF260604_mPaintSpace2_N1A",
    dmso_wells = c("2_2","3_9","4_2","5_9","6_2",
                   "7_9","8_2","9_9","10_2","11_9",
                   "12_2","13_9","14_2","15_9",
                   "16_2","17_9","18_2","19_9",
                   "20_2","21_9","22_2","23_9")
  ),
  N1B = list(
    file_name = "SF260604_mPaintSpace2_N1B",
    dmso_wells = c("2_2","3_9","4_2","5_9","6_2",
                   "7_9","8_2","9_9","10_2","11_9",
                   "12_2","13_9","14_2","15_9",
                   "16_2","17_9","18_2","19_9",
                   "20_2","21_9","22_2","23_9")
  ),
  N2A = list(
    file_name = "SF260604_mPaintSpace2_N2A",
    dmso_wells = c("2_2","3_9","4_2","5_9","6_2",
                   "7_9","8_2","9_9","10_2","11_9",
                   "12_2","13_9","14_2","15_9",
                   "16_2","17_9","18_2","19_9",
                   "20_2","21_9","22_2","23_9")
  ),
  N2B = list(
    file_name = "SF260604_mPaintSpace2_N2B",
    dmso_wells = c("2_2","3_9","4_2","5_9","6_2",
                   "7_9","8_2","9_9","10_2","11_9",
                   "12_2","13_9","14_2","15_9",
                   "16_2","17_9","18_2","19_9",
                   "20_2","21_9","22_2","23_9")
  ),
  N3A = list(
    file_name = "SF260701_mPaintSpace2_N3A",
    dmso_wells = c("2_2","3_9","4_2","5_9","6_2",
                   "7_9","8_2","9_9","10_2","11_9",
                   "12_2","13_9","14_2","15_9",
                   "16_2","17_9","18_2","19_9",
                   "20_2","21_9","22_2","23_9")
  ),
  N3B = list(
    file_name = "SF260701_mPaintSpace2_N3B",
    dmso_wells = c("2_2","3_9","4_2","5_9","6_2",
                   "7_9","8_2","9_9","10_2","11_9",
                   "12_2","13_9","14_2","15_9",
                   "16_2","17_9","18_2","19_9",
                   "20_2","21_9","22_2","23_9")
  )
)
perc_width <- 6.2
perc_height <- 6.4
file_name <- "mPaintSpace2_N1_N2_N3"
pos_control <- "CCCP_20"
pastel_cols <- lighten(c("#440154FF", "#238A8DFF", "#FDE725FF"), amount = 0.3)
# create a function to load data ####
load_data <- function(batch_name, file_name, dmso_wells) {
  # load plate drift corrected data as df
  df <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_platedrift_corr.csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames as WELL_BATCH
  rownames(df) <- df$V1
  df$V1 <- NULL
  # load raw data as df_raw
  df_raw <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_data_tidy.csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames as WELL_BATCH
  rownames(df_raw) <- df_raw$V1
  df_raw$V1 <- NULL
  # load metadata as meta
  meta <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_meta_tidy.csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames as WELL_BATCH
  rownames(meta) <- meta$V1
  meta$V1 <- NULL
  # load best fits table
  fits <- as.data.frame(
    fread(
      paste(
        "outputs/data/", file_name, "_platedrift_fits.csv", sep = ""), 
      header = TRUE)
  )
  # keep rownames as WELL_BATCH
  rownames(fits) <- fits$V1
  fits$V1 <- NULL
  # make sure ID is numeric for column indexing later
  fits$ID <- as.integer(fits$ID)
  # return list of drift corrected data, raw data, metadata, and drift fits
  return(list(
    batch_name = batch_name,
    file_name = file_name,
    df = df,
    df_raw = df_raw,
    meta = meta,
    fits = fits,
    dmso_wells = dmso_wells
  ))
}
# run function to load data ####
# make a super list of all N, to be used in lapply() loops
batches <- imap(
  batches_info,
  function(batch_info, batch_name) {
    load_data(
      batch_name = batch_name,
      file_name = batch_info$file_name,
      # use explicit DMSO wells variable in case DMSO wells are different 
      dmso_wells = batch_info$dmso_wells
    )
  }
)
# create function to select feature IDs for plate drift QC plots ####
# this function will actually be run wrapped in the later plotting function
select_ids <- function(fit_table, n_random = 4, seed = 42) {
  # set seed so "randomly" selected features are reproducible
  set.seed(seed)
  # subset feature IDs with poly drift
  lm2 <- filter(fit_table, best_grad2 != 0)
  # slice select 4 examples with poly drift
  lm2_slice <- bind_rows(
    # 1. decreasing slope
    slice_min(lm2, best_grad2, n = 1, with_ties = FALSE),
    # 2. increasing slope
    slice_max(lm2, best_grad2, n = 1, with_ties = FALSE),
    # 3. strong fit (low p)
    slice_min(lm2, best_p, n = 1, with_ties = FALSE),
    # 4. weak fit (p just below p_sig/ 0.05)
    slice_max(lm2, best_p, n = 1, with_ties = FALSE)
  )
  # select feature IDs with linear drift
  lm1 <- filter(fit_table, best_grad2 == 0, best_p != 0)
  # slice select 4 examples with linear drift
  lm1_slice <- bind_rows(
    # 1. decreasing slope
    slice_min(lm1, best_grad1, n = 1, with_ties = FALSE),
    # 2. increasing slope
    slice_max(lm1, best_grad1, n = 1, with_ties = FALSE),
    # 3. strong fit (low p)
    slice_min(lm1, best_p, n = 1, with_ties = FALSE),
    # 4. weak fit (p just below p_sig/ 0.05)
    slice_max(lm1, best_p, n = 1, with_ties = FALSE)
  )
  # select feature IDs with no drift
  lm0 <- filter(fit_table, best_p == 0)
  # slice select 4 examples with no drift
  lm0_slice <- slice_sample(lm0, n = min(n_random, nrow(lm0)))
  # return vector of unique feature IDs to plot
  unique(c(lm2_slice$ID, lm1_slice$ID, lm0_slice$ID))
}
# create function to generate plate drift QC plots of dmso + pos_control####
qc_plots <- function(batch_obj, plot_condition = "DMSO_0") {
  # df is data after plate drift correction
  df <- batch_obj$df
  # df_raw is raw data before plate drift correction
  df_raw <- batch_obj$df_raw
  meta <- batch_obj$meta
  fits <- batch_obj$fits
  dmso_wells <- batch_obj$dmso_wells
  # select for dmso
  # , , preserves the row order in the selection of data from meta
  meta_dmso <- filter(meta, Well %in% dmso_wells)
  df_dmso <- df[rownames(meta_dmso), , drop = FALSE]
  df_raw_dmso <- df_raw[rownames(meta_dmso), , drop = FALSE]
  # select for plot condition 
  # _p will be identical to _dmso if plot condition is just DMSO
  # _p grey points will get covered by black _dmso points if identical
  meta_p <- filter(meta, Condition == plot_condition)
  df_p <- df[rownames(meta_p), , drop = FALSE]
  df_raw_p <- df_raw[rownames(meta_p), , drop = FALSE]
  # select for IDs to plot
  qc_ids <- select_ids(fits)
  # make plots list
  pdf_list <- lapply(qc_ids, function(ID) {
    feature_name <- colnames(df)[ID]
    # before correction
    df_before_dmso <- data.frame(Order = meta_dmso$Order, 
                                 y = df_raw_dmso[[ID]])
    df_before_p <- data.frame(Order = meta_p$Order, 
                              y = df_raw_p[[ID]])
    # after correction
    df_after_dmso <- data.frame(Order = meta_dmso$Order, 
                                y = df_dmso[[ID]])
    df_after_p <- data.frame(Order = meta_p$Order, 
                             y = df_p[[ID]])
    # assign drift model type 
    if (fits$best_grad2[ID] != 0) {
      formula <- y ~ x + I(x^2)
      col_line <- pastel_cols[1]
    } else if (fits$best_grad1[ID] != 0) {
      formula <- y ~ x
      col_line <- pastel_cols[2]
    } else {
      formula <- y ~ 1
      col_line <- pastel_cols[3]
    }
    # make sure y limits match for before and after, and dont cutoff points
    ymax <- max(c(df_before_dmso$y,
                  df_before_p$y,
                  df_after_dmso$y,
                  df_after_p$y),
                na.rm = TRUE)
    ymin <- min(c(df_before_dmso$y,
                  df_before_p$y,
                  df_after_dmso$y,
                  df_after_p$y),
                na.rm = TRUE)
    # plot before drift correction for given ID (y)
    plot1 <- ggplot() +
      geom_point(data = df_before_p,
                 aes(x = Order, y = y),
                 colour = "grey70",
                 size = 0.4) +
      geom_point(data = df_before_dmso,
                 aes(x = Order, y = y),
                 colour = "black",
                 size = 0.4) +
      stat_smooth(data = df_before_dmso,
                  aes(x = Order, y = y),
                  method = "lm",
                  formula = formula,
                  se = TRUE,
                  colour = col_line,
                  fill = col_line,
                  alpha = 0.2,
                  linewidth = 0.8) +
      coord_cartesian(ylim = c(ymin, ymax)) +
      scale_x_continuous(breaks = c(0, 100, 200, 300)) +
      theme_pubr() +
      labs(title = "Before Drift Correction",
           x = "Order") +
      theme(axis.text = element_text(size = 6),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_blank(),
            plot.title = element_text(size = 8, hjust = 0.5))
    # plot after drift correction for given ID (y)
    plot2 <- ggplot() +
      geom_point(data = df_after_p,
                 aes(x = Order, y = y),
                 colour = "grey70",
                 size = 0.4) +
      geom_point(data = df_after_dmso,
                 aes(x = Order, y = y),
                 colour = "black",
                 size = 0.4) +
      stat_smooth(data = df_after_dmso,
                  aes(x = Order, y = y),
                  method = "lm",
                  formula = y ~ 1,
                  se = TRUE,
                  colour = pastel_cols[3],
                  fill = pastel_cols[3],
                  alpha = 0.2,
                  linewidth = 0.8) +
      coord_cartesian(ylim = c(ymin, ymax)) +
      scale_x_continuous(breaks = c(0, 100, 200, 300)) +
      theme_pubr() +
      labs(title = "After Drift Correction",
           x = "Order") +
      theme(axis.text = element_text(size = 6),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_blank(),
            plot.title = element_text(size = 8, hjust = 0.5))
    # define a main title (feature name) for combined before and after 
    title <- ggdraw() +
      draw_label(str_wrap(feature_name, width = 40), size = 10)
    # combine before and after together
    plot_grid(title,
              plot_grid(plot1, plot2),
              ncol = 1,
              rel_heights = c(0.3, 1))
  })
  # final_plot for each N
  # final_plot is 4 cols x 3 rows
  # each subplot is before and after for a feature
  final_plot <- plot_grid(plotlist = pdf_list, ncol = 4)
  batch_obj$qc_ids  <- qc_ids
  batch_obj$qc_plot <- final_plot
  return(batch_obj)
}
# run function to generate plate drift QC plots (for dmso only) ####
plots_dmso <- lapply(batches, qc_plots)
# open dmso qc plots in viewer
walk(plots_dmso, ~ print(.x$qc_plot))
# run function to generate plate drift QC plots (for dmso and pos_control) ####
plots_p <- lapply(batches, qc_plots, plot_condition = pos_control)
# open positive control qc plots in viewer
walk(plots_p, ~ print(.x$qc_plot))
# create function to calculate % of drift corrected features per N ####
perc_drift <- function(batch_obj, batch_name) {
  fits <- batch_obj$fits
  total <- nrow(fits)
  poly   <- sum(fits$best_grad2 != 0)
  linear <- sum(fits$best_grad2 == 0 & fits$best_grad1 != 0)
  none   <- sum(fits$best_p == 0)
  data.frame(
    Batch = batch_name,
    Drift = c("Polynomial", "Linear", "None"),
    Count = c(poly, linear, none),
    Percentage = c(poly, linear, none) / total * 100
  )
}
# run function to calculate % of drift corrected features per N ####
drift_df <- imap_dfr(
  batches,
  function(batch_obj, batch_name) {
    perc_drift(batch_obj, batch_name)
  }
)
drift_df$Batch <- factor(
  drift_df$Batch,
  levels = names(batches)
)
drift_df$Drift <- factor(drift_df$Drift,
                         levels = c("Polynomial","Linear","None"))
# plot stacked bar of % of drift corrected features in each N ####
drift_perc_plot <- ggplot(drift_df,
                                aes(x = Batch,
                                    y = Percentage,
                                    fill = Drift)) +
  geom_bar(stat = "identity",
           colour = "black",
           width = 0.8,
           position = position_stack(reverse = TRUE)) +
  geom_label(aes(label = sprintf("%.2f%%", Percentage)),
             position = position_stack(vjust = 0.5, reverse = TRUE),
             colour = "black",
             fill = "white",
             label.size = 0.4,
             size = 4) +
  scale_fill_manual(values = c(
    "Polynomial"   = pastel_cols[1],
    "Linear" = pastel_cols[2],
    "None"     = pastel_cols[3]
  )) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100)) +
  labs(y = "Percentage of Features (%)",
       x = "Batch",
       fill = "Plate drift \nModel of best fit") +
  theme_pubr() +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 15),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.title = element_text(hjust = 0.5, size = 8),
    legend.text = element_text(size = 7),
    # Remove extra padding
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    # Critical for true alignment
    plot.margin = margin(t = 5, r = 5, b = 5, l = 0)
  )
drift_perc_plot
# save plots ####
# save dmso drift qc plots, for however many number of N
iwalk(plots_dmso, function(batch_obj, batch_name) {
  ggsave(
    paste(
      "outputs/figures/",
      batch_obj$file_name,
      "_qc_plot.pdf",
      sep = ""
    ),
    batch_obj$qc_plot,
    width = 13.4,
    height = 6.4,
    units = "in",
    dpi = 300
  )
})
# save dmso + pos control drift qc plots, for however many number of N
iwalk(plots_p, function(batch_obj, batch_name) {
  ggsave(
    paste(
      "outputs/figures/",
      batch_obj$file_name,
      "_qc_p_plot.pdf",
      sep = ""
    ),
    batch_obj$qc_plot,
    width = 13.4,
    height = 6.4,
    units = "in",
    dpi = 300
  )
})
ggsave(
  paste(
    "outputs/figures/", file_name, "_drift_perc_plot.pdf", sep = ""),
  drift_perc_plot,
  width = perc_width,
  height = perc_height,
  units = "in",
  dpi = 300
)
rm(list = ls())
