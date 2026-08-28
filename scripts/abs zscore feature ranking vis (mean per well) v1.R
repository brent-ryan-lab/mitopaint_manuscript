# Title: abs zscore feature ranking vis (mean per well) v1
# Step: 6.3
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 18-08-2026

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
file_name <- "mPaintSpace2_N1_N2_N3"
redu_state <- "redu"
integrate_state <- "integrated"
ctrl_cond <- "DMSO_0"
plot_cond <- c("Nigericin_3", "Oligomycin_10", "CCCP_20", "Rapamycin_10", "Valinomycin_5",
               "Cyclosporin A_10", "ROT_3", "CQ_20", "BAM15_10", "MitoQ_3", "Nocodazole_10",
               "Cytochalasin D_5")
contrast_all <- FALSE
rank <- 10
# load data ####
# load data as df
if (redu_state == "redu") {
  df <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_data_", integrate_state ,"_redu.csv", sep = ""), 
      header = TRUE)
  )
} else {
  df <- as.data.frame(
    fread(
      paste(
        "data/processed/", file_name, "_data_", integrate_state ,".csv", sep = ""), 
      header = TRUE)
  )
}
# keep rownames as WELL_BATCH
rownames(df) <- df$V1
df$V1 <- NULL
# load metadata as meta
meta <- as.data.frame(
  fread(
    paste(
      "data/processed/", file_name, "_meta_", integrate_state, ".csv", sep = ""), 
    header = TRUE)
)
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
        adjust = "bonferroni"
      )
    )
  } else {
    # if contrast_all = F, only compare each plot_cond to ctrl_cond
    as.data.frame(
      contrast(
        emm,
        method = "trt.vs.ctrl",
        ref = ctrl_cond,
        adjust = "bonferroni"
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
# add absolute effect size column
abs_df <- abs_df|> mutate(abs_estimate = abs(estimate))
# sort by absolute effect size column
abs_df <- abs_df|> arrange(desc(abs_estimate)) 
# split abs_df into a list of df for each unique contrast
abs_df <- abs_df |> split(abs_df$contrast)
# create function for plotting ranking by absolute difference ####
abs_plot <- function(df,
                     rank = 10,
                     # optional rank_by and rank_dir variables if wanting to rank by something else (eg. p.value)
                     rank_by = "abs_estimate",
                     rank_dir = "desc",
                     fill_col = "steelblue",
                     x_lab = "Absolute Estimated Marginal Means\nof Robust Z-score") {
  # pull data
  if (rank_dir == "desc") {
    plot_df <- df |>
      arrange(desc(.data[[rank_by]])) |>
      # slice by rank (eg. top 10 features, if rank <- 10)
      slice_head(n = rank)
  } else if (rank_dir == "asc") {
    plot_df <- df |>
      arrange(.data[[rank_by]]) |>
      # slice by rank (eg. top 10 features, if rank <- 10)
      slice_head(n = rank)
  } else {
    stop("rank_dir must be 'asc' or 'desc'")
  }
  plot_df$feature <- factor(plot_df$feature, levels = rev(plot_df$feature))
  # plot
  ggplot(plot_df, aes(x = abs_estimate, y = feature)) +
    geom_col(
      fill = fill_col,
      colour = "black",
      width = 0.8
    ) +
    # error bar is SE of estimate marginal means
    geom_errorbarh(
      aes(
        xmin = abs_estimate - SE,
        xmax = abs_estimate + SE
      ),
      height = 0.4,
      colour = "black"
    ) +
    labs(
      x = x_lab,
      # y is feature name (if name is too long it wraps to next line)
      y = stringr::str_wrap(NULL, width = 35),
      title = unique(plot_df$contrast)
    ) +
    # tidy theme
    theme_pubr() +
    theme(
      # fix dimensions of plot (square)
      aspect.ratio = 1,
      # keep left width wide to avoid squishing plot from long feature names
      plot.margin = margin(5, 5, 5, 10),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9),
      axis.title.x = element_text(size = 10),
      plot.title = element_text(
        hjust = 0.5,
        size = 12,
        face = "bold"),
      panel.grid = element_blank(),
      legend.position = "none"
    )
}
# run function for absolute difference bar plot ####
plots <- map2(
  abs_df,
  names(abs_df),
  function(df, contrast_name) {
      abs_plot(
      df = df,
      rank = rank,
      rank_by = "abs_estimate",
      rank_dir = "desc"
    ) +
      labs(title = contrast_name)
  }
)
plots
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
rm(list = ls())
