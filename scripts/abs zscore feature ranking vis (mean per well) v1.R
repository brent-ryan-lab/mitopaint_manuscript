# Title: abs zscore feature ranking vis (mean per well) v2
# Step: 6.3
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 23-11-2026

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
file_name <- "cPaintSpace2_N1_N2_N3"
redu_state <- "redu"
integrate_state <- "integrated"
ctrl_cond <- "DMSO_0"
plot_cond <- c("Nigericin_3", "Oligomycin_10", "CCCP_20", "Rapamycin_10", "Valinomycin_5",
               "Cyclosporin A_10", "ROT_3", "CQ_20", "BAM15_10", "MitoQ_3", "Nocodazole_10",
               "Cytochalasin D_5")
contrast_all <- FALSE
rank <- 10
avg_profile <- FALSE
post_hoc <- "bonferroni"
chan_cols <- c(
  `Concanavalin A 488` = "yellowgreen",
  `MitoTracker Deep Red` = "#DC267F",
  `Phalloidin 568 WGA 555` = "#FFB000",
  `SYTO 14` = "#23CC86",
  `Hoechst` = "#7E8CF3FF"
)
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
# keep conditions of interest ####
# make factors
meta$Batch <- factor(meta$Batch)
meta$Condition <- factor(meta$Condition)
# keep only ctrl_cond + plot_cond
meta <- meta[meta$Condition %in% c(ctrl_cond, plot_cond), , drop = FALSE]
df <- df[rownames(df) %in% rownames(meta), ]
# create function for linear mixed effect model, to quantify absolute difference ####
lmme_fit <- function(feature_name,
                     df,
                     meta,
                     ctrl_cond = "DMSO_0",
                     plot_cond,
                     contrast_all = FALSE) {
  # pull data
  dat <- cbind(
    meta,
    value = df[[feature_name]]
  )
  # keep a sensible order for the condition factor
  dat$Condition <- factor(dat$Condition, levels = unique(c(ctrl_cond, plot_cond)))
  # fit is lmme of Condition, accounting for random effect from batch
  # suppress warnings that Batch has no effect "boundary (singular) fit: see help('isSingular')"
  fit <- suppressWarnings(
    suppressMessages(
      lmer(value ~ Condition + (1 | Batch), data = dat)
    )
  )
  emm <- emmeans(fit, ~ Condition)
  contr_df <- if (contrast_all) {
    # if contrast_all = T, calculate all unique pairwise comparisons, no redundant reverse directions
    as.data.frame(
      pairs(
        emm,
        adjust = "none"
      )
    )
  } else {
    # if contrast_all = F, only compare each plot_cond to ctrl_cond
    as.data.frame(
      contrast(
        emm,
        method = "trt.vs.ctrl",
        ref = ctrl_cond,
        adjust = "none"
      )
    )
  }
  contr_df |>
    mutate(feature = feature_name, .before = 1)
}
# run loop for linear mixed effect model, to quantify absolute difference ####
# run over every feature
abs_df <- imap_dfr(
  colnames(df),
  function(feat, i) {
    # print a progress message in console (since it takes a while to calculate for all feature)
    message(sprintf("[%d/%d] Calculating %s", i, length(colnames(df)), feat))
    lmme_fit(
      feature_name = feat,
      df = df,
      meta = meta,
      ctrl_cond = ctrl_cond,
      plot_cond = plot_cond,
      contrast_all = contrast_all
    )
  }
)
# apply post hoc correction ####
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
# add absolute effect size column ####
# calculate absolute effect size
abs_df <- abs_df |>
  dplyr::mutate(
    abs_estimate = abs(
      estimate
    )
  ) |>
  dplyr::arrange(
    dplyr::desc(abs_estimate)
  )

# split by contrast
abs_df <- abs_df |>
  split(abs_df$contrast)
# create function for plotting ranking by absolute difference ####
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
# plot ranked by adjusted p-value ####
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
# align plots ####
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
