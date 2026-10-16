# Title: abs zscore feature ranking vis-NN (mean per well) v1
# Step: 6.3.2
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 15-10-2026

# load packages ####
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(purrr)
library(lme4)
library(lmerTest)
library(emmeans)
library(tibble)
library(viridis)
library(colorspace)
# set variables ####
file_name <- "mPaintFDA_N1_N2_N3_N4_N5_N6_N7_N8"
redu_state <- "redu"
integrate_state <- "integrated"
ctrl_cond <- "DMSO_0"
ctrl_cond_col <- "Condition"
plot_cond <- c(seq(from = 0, to = 27, by = 1))
plot_cond_col <- "UMAP_NN"
post_hoc <- "bonferroni"
contrast_all <- FALSE
rank <- 10
avg_profile <- TRUE
chan_cols <- c("TMRM" = "#FFB000", "CellRox" = "#DC267F", "mt-Keima pH7" = "#23CC86", "mt-Keima pH4" = "#FE6100")
# load data ####
# load data as df
if (avg_profile == TRUE) {
  df <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_", integrate_state ,"_redu_avg_data.csv", sep = ""), 
      header = TRUE)
  )
} else {
  df <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_data_", integrate_state ,"_redu.csv", sep = ""), 
      header = TRUE)
  )
}
# keep rownames as WELL_BATCH
rownames(df) <- df$V1
df$V1 <- NULL
# load metadata as meta
if (avg_profile == TRUE) {
  meta <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_avg_dimred_meta.csv", sep = ""), 
      header = TRUE)
  )
} else {
  meta <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_", integrate_state, "_", redu_state, "_dimred_meta.csv", sep = ""), 
      header = TRUE)
  ) 
}
# keep rownames as WELL_BATCH
rownames(meta) <- meta$V1
meta$V1 <- NULL
# align data ####
# ensure df and meta are exactly aligned
meta <- meta[rownames(df), , drop = FALSE]
meta$Batch <- factor(meta$Batch)
# create function for linear mixed effect model, to quantify absolute difference (estimated marginal means) ####
# function to compare one cluster against DMSO #
lmme_cluster_vs_ctrl <- function(
    feature_name,
    cluster_value,
    df,
    meta,
    plot_cond_col = "UMAP_NN",
    ctrl_cond_col = "Condition",
    ctrl_cond = "DMSO_0"
) {
  # align metadata and feature values
  dat <- meta |>
    dplyr::mutate(
      value = df[[feature_name]],
      row_id = rownames(meta)
    )
  # wells belonging to cluster being tested
  cluster_dat <- dat |>
    dplyr::filter(
      .data[[plot_cond_col]] == cluster_value
    ) |>
    dplyr::mutate(
      ModelGroup = paste0(
        plot_cond_col,
        "_",
        cluster_value
      )
    )
  # all DMSO wells
  # exclude DMSO wells already included in cluster_dat
  ctrl_dat <- dat |>
    dplyr::filter(
      .data[[ctrl_cond_col]] == ctrl_cond,
      !row_id %in% cluster_dat$row_id
    ) |>
    dplyr::mutate(
      ModelGroup = ctrl_cond
    )
  # combine cluster and control observations
  model_dat <- dplyr::bind_rows(
    ctrl_dat,
    cluster_dat
  ) |>
    dplyr::filter(
      is.finite(value)
    )
  # DMSO is reference
  model_dat$ModelGroup <- relevel(
    factor(model_dat$ModelGroup),
    ref = ctrl_cond
  )
  model_dat$Batch <- factor(
    model_dat$Batch
  )
  # mixed model
  fit <- suppressWarnings(
    suppressMessages(
      lmerTest::lmer(
        value ~ ModelGroup + (1 | Batch),
        data = model_dat
      )
    )
  )
  # EMM
  emm <- emmeans::emmeans(
    fit,
    ~ ModelGroup
  )
  # cluster vs DMSO
  contr_df <- as.data.frame(
    emmeans::contrast(
      emm,
      method = "trt.vs.ctrl",
      ref = ctrl_cond,
      adjust = "none"
    )
  )
  # add identifiers
  contr_df |>
    dplyr::mutate(
      feature = feature_name,
      cluster = cluster_value,
      n_cluster = nrow(cluster_dat),
      n_control = nrow(ctrl_dat),
      .before = 1
    )
}
# run loop for linear mixed effect model, to quantify absolute difference ####
# identify UMAP clusters
plot_cond <- sort(
  unique(meta[[plot_cond_col]])
)
plot_cond <- plot_cond[
  !is.na(plot_cond)
]
# run loop
# note: computationally intensive for each NN group (approx ~5mins per NN group to calculate for ~1000 features)
abs_df <- purrr::map_dfr(
  plot_cond,
  function(cluster_value) {
    message(
      "Calculating ",
      plot_cond_col,
      " ",
      cluster_value
    )
    purrr::map_dfr(
      colnames(df),
      function(feat) {
        lmme_cluster_vs_ctrl(
          feature_name = feat,
          cluster_value = cluster_value,
          df = df,
          meta = meta,
          plot_cond_col = plot_cond_col,
          ctrl_cond_col = ctrl_cond_col,
          ctrl_cond = ctrl_cond
        )
      }
    )
  }
)
# apply post hoc correction
# adjust p-values only if post_hoc is specified
if (!is.null(post_hoc) && !is.na(post_hoc)) {
  abs_df <- abs_df |>
    dplyr::mutate(
      p.value_raw = p.value,
      p.value = p.adjust(
        p.value_raw,
        method = post_hoc
      )
    )
}
# calculate absolute effect size and arrange
abs_df <- abs_df |>
  dplyr::mutate(
    abs_estimate = abs(estimate)
  ) |>
  dplyr::arrange(
    cluster,
    dplyr::desc(abs_estimate)
  )
# add absolute effect size column
abs_df <- abs_df|> mutate(abs_estimate = abs(estimate))
# sort by absolute effect size column
abs_df <- abs_df|> arrange(desc(abs_estimate)) 
# split abs_df into a list of df for each unique contrast
abs_df <- abs_df |> split(abs_df$contrast)
# create function for plotting ranking by absolute difference ####
abs_plot <- function(
    df,
    rank = 10,
    rank_by = "abs_estimate",
    rank_dir = "desc",
    x_by = "abs_estimate",
    chan_cols = NULL,
    x_lab = "Absolute Estimated Marginal Means\nof Robust Z-score"
) {
  # rank/select features using ORIGINAL values
  if (rank_dir == "desc") {
    plot_df <- df |>
      dplyr::arrange(
        dplyr::desc(.data[[rank_by]])
      ) |>
      dplyr::slice_head(n = rank)
  } else if (rank_dir == "asc") {
    plot_df <- df |>
      dplyr::arrange(
        .data[[rank_by]]
      ) |>
      dplyr::slice_head(n = rank)
  } else {
    stop("rank_dir must be 'asc' or 'desc'")
  }
  # preserve ranking order in plot
  plot_df$feature <- factor(
    plot_df$feature,
    levels = rev(plot_df$feature)
  )
  # determine bar colours
  if (is.null(chan_cols)) {
    plot_df$bar_fill <- "steelblue"
  } else {
    # default features without matching channel
    plot_df$bar_fill <- "grey60"
    # assign channel-specific colours
    for (chan in names(chan_cols)) {
      matches <- stringr::str_detect(
        as.character(plot_df$feature),
        stringr::fixed(chan)
      )
      plot_df$bar_fill[matches] <- chan_cols[[chan]]
    }
  }
  # determine x-axis values
  if (x_by == "p.value") {
    # find minimum NON-ZERO p-value from the full contrast dataframe
    min_pos_p <- min(
      df$p.value[
        df$p.value > 0 &
          is.finite(df$p.value)
      ],
      na.rm = TRUE
    )
    # replace zero only for plotting
    # original p.value remains unchanged
    plot_df$p_plot <- dplyr::if_else(
      plot_df$p.value == 0,
      min_pos_p,
      plot_df$p.value
    )
    # transform for plotting
    plot_df$x_plot <- -log10(
      plot_df$p_plot
    )
  } else {
    plot_df$x_plot <- plot_df[[x_by]]
  }
  # make plot
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = x_plot,
      y = feature
    )
  ) +
    ggplot2::geom_col(
      ggplot2::aes(fill = bar_fill),
      colour = "black",
      width = 0.8
    ) +
    ggplot2::scale_fill_identity()
  # only add SE error bars when plotting effect size
  if (x_by == "abs_estimate") {
    p <- p +
      ggplot2::geom_errorbarh(
        ggplot2::aes(
          xmin = abs_estimate - SE,
          xmax = abs_estimate + SE
        ),
        height = 0.4,
        colour = "black"
      )
    p <- p +
      # significance threshold at p = 0.05
      ggplot2::geom_vline(
        xintercept = 2,
        linetype = "dashed",
        colour = "black"
      ) 
  }
  # p-value axis formatting and significance threshold
  if (x_by == "p.value") {
    p <- p +
      # significance threshold at p = 0.05
      ggplot2::geom_vline(
        xintercept = -log10(0.05),
        linetype = "dashed",
        colour = "black"
      ) +
      # format x-axis labels to 2 decimal places
      ggplot2::scale_x_continuous(
        labels = scales::label_number(
          accuracy = 0.01
        )
      )
  }
  # finish plot
  p +
    ggplot2::labs(
      x = x_lab,
      y = NULL,
      title = unique(plot_df$contrast)
    ) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.margin = ggplot2::margin(5, 5, 5, 10),
      axis.text.y = ggplot2::element_text(size = 9),
      axis.text.x = ggplot2::element_text(size = 9),
      axis.title.x = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 12,
        face = "bold"
      ),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    )
}
# run function for absolute difference bar plot ####
plots <- purrr::map2(
  abs_df,
  names(abs_df),
  function(df, contrast_name) {
    abs_plot(
      df = df,
      rank = rank,
      rank_by = "abs_estimate",
      rank_dir = "desc",
      x_by = "abs_estimate",
      chan_cols = chan_cols,
      x_lab = "Absolute EMM Difference\n(Robust Z-score)"
    ) +
      ggplot2::labs(
        title = contrast_name
      )
  }
)
plots
plots2 <- purrr::map2(
  abs_df,
  names(abs_df),
  function(df, contrast_name) {
    abs_plot(
      df = df,
      rank = rank,
      rank_by = "p.value",
      rank_dir = "asc",
      x_by = "p.value",
      chan_cols = chan_cols,
      x_lab = "-log10(adjusted p-value)\nfor EMM Difference"
    ) +
      ggplot2::labs(
        title = contrast_name
      )
  }
)
plots2
# save data ####
iwalk(
  abs_df,
  function(df, contrast_name) {
    safe_name <- stringr::str_replace_all(contrast_name, "[^[:alnum:]]+", "_")
    write.csv(
      df,
      file = paste0("outputs/data/", file_name, "_", safe_name, "_abs_emm", ".csv")
    )
  }
)
# save plots ####
# align x and y axis of plots
plots <- cowplot::align_plots(
  plotlist = plots,
  align = "hv",
  axis = "tblr"
)
plots2 <- cowplot::align_plots(
  plotlist = plots2,
  align = "hv",
  axis = "tblr"
)
iwalk(
  plots,
  function(plot, plot_name) {
    safe_name <- stringr::str_replace_all(plot_name, "[^[:alnum:]]+", "_")
    ggsave(
      filename = paste0("outputs/figures/", file_name, "_", safe_name, "_abs_emm", ".pdf"),
      plot = plot,
      width = 8,
      height = 3
    )
  }
)
iwalk(
  plots2,
  function(plot, plot_name) {
    safe_name <- stringr::str_replace_all(plot_name, "[^[:alnum:]]+", "_")
    ggsave(
      filename = paste0("outputs/figures/", file_name, "_", safe_name, "_p_emm", ".pdf"),
      plot = plot,
      width = 8,
      height = 3
    )
  }
)
#rm(list = ls())
