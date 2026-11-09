# Title: hits- mahalanobis distance mitopaint vis (mean per well) v2
# Step: 7.2
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 09-11-26

# load packages ####
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
# set file variables ####
file_name <- "mPaintFDA_N1_N2_N3_N4_N5_N6_N7_N8"
integrate_state <- "integrated"
redu_state <- "redu"
pc_use <- 10
grid_width <- 2
plot_width <- 8
plot_height <- 4
avg_profile <- TRUE
# create function to load data ####
load_data <- function(md_chi_file,
                      md_lmme_file,
                      integrate_state) {
  # load mahalanobis distance (md)
  md_chi_mean <- as.data.frame(
    fread(
      paste("data/processed/", file_name, "_mahal_chi_mean_PC", pc_use, ".csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames as WELL_BATCH
  rownames(md_chi_mean) <- md_chi_mean$V1
  md_chi_mean$V1 <- NULL
  # load md_lmme_p
  md_lmme_p <- as.data.frame(
    fread(
      paste("data/processed/", file_name, "_mahal_lmme_p_PC", pc_use,".csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames
  rownames(md_lmme_p) <- md_lmme_p$V1
  md_lmme_p$V1 <- NULL  
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
  # return list
  return(list(
    meta = meta,
    md_chi_mean = md_chi_mean,
    md_lmme_p = md_lmme_p
  ))
}
# load data ####
data <- load_data(md_chi_file,
                  md_lmme_file,
                  integrate_state)
# plot all hits by lmme p ####
# filter hits to those with p < 0.05
plot_df <- data$md_lmme_p |>
  # replace exact 0 p.value with slightly smaller than min value
  mutate(
    p.value = if_else(
      p.value == 0,
      0.9 * min(p.value[p.value > 0], na.rm = TRUE),
      p.value
    )
  ) |>
  filter(p.value < 0.05) |>
  arrange(p.value) |>
  mutate(
    Condition = factor(Condition, levels = rev(unique(Condition)))
  )
# initialize plots list
plots <- list()
# hits_lmme is dot plot of compounds with significant p from lmme test
plots$hits_lmme <- ggplot(plot_df, 
                         # x = -log(p.value), higher x is stronger hit
                         # y = condition, in desc order (top hits at top)
                         aes(x = -log(p.value), y = Condition)) +
  geom_point() +
  # red dashed line at p = 0.05
  geom_vline(xintercept = -log(0.05), linetype = "dashed", colour = "red") +
  labs(x = "-log(p.value)",
       y = "Condition",
       title = "Mahalanobis distance hits\n(linear mixed effect model)") +
  theme_pubr() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    panel.grid = element_blank()
  )
plots$hits_lmme
# plot all hits by chi cutoff ####
# save cutoff threshold for plot annotation
md_cutoff <- qchisq(0.95, df = pc_use)
# save DMSO avg md for plot annotation
dmso_md_mean <- data$md_chi_mean |>
  dplyr::filter(Condition == "DMSO_0") |>
  dplyr::pull(md2_mean)
# filter hits to those with md > cutoff
plot_df <- data$md_chi_mean |>
  filter(md2_mean > md_cutoff) |>
  arrange(desc(md2_mean)) |>
  mutate(
    Condition = factor(Condition, levels = rev(unique(Condition)))
  )
# hits_md is dot plot of compounds with significant p from chi squared test
plots$hits_chi <- ggplot(
  plot_df,
  aes(
    x = log(md2_mean),
    y = Condition
  )
) +
  # SEM bars
  # rows with NA SEM are automatically omitted
  geom_errorbarh(
    aes(
      xmin = log(md2_mean - md2_sem),
      xmax = log(md2_mean + md2_sem)
    ),
    height = 0.3,
    colour = "black",
    na.rm = TRUE
  ) +
  # mean points
  geom_point() +
  # red dashed line at chi-squared significance threshold
  geom_vline(
    xintercept = log(md_cutoff),
    linetype = "dashed",
    colour = "red"
  ) +
  # grey dashed line at mean DMSO Mahalanobis distance
  geom_vline(
    xintercept = log(dmso_md_mean),
    linetype = "dashed",
    colour = "grey60"
  ) +
  labs(
    x = "log(md2_mean)",
    y = "Condition",
    title = "Mahalanobis distance hits\n(Chi squared test)"
  ) +
  theme_pubr() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(
      hjust = 0.5,
      size = 10,
      face = "bold"
    ),
    panel.grid = element_blank()
  )
plots$hits_chi
# create function to plot grid of lmme_p by compound ####
# identify compounds with more than one unique concentration
multi_conc_compounds <- data$md_lmme_p |>
  dplyr::filter(
    !is.na(Compound),
    !is.na(Concentration)
  ) |>
  dplyr::group_by(Compound) |>
  dplyr::summarise(
    n_concentrations = dplyr::n_distinct(Concentration),
    .groups = "drop"
  ) |>
  dplyr::filter(
    n_concentrations > 1
  ) |>
  dplyr::pull(Compound)
# subset data to compounds with multiple concentrations
multi_conc_df <- data$md_lmme_p |>
  dplyr::filter(
    Compound %in% multi_conc_compounds
  )
# replace zeros once for determining limits
p_vals <- multi_conc_df$p.value
min_pos_p <- min(
  p_vals[
    p_vals > 0 &
      is.finite(p_vals)
  ],
  na.rm = TRUE
)
zero_replacement <- 0.9 * min_pos_p
# shared y limits based only on multi-concentration compounds
y_vals <- -log(
  ifelse(
    p_vals == 0,
    zero_replacement,
    p_vals
  )
)
y_vals <- y_vals[
  is.finite(y_vals)
]
# consistent y limits for all plots
y_limits <- range(
  y_vals,
  na.rm = TRUE
) * c(0.95, 1.05)

# create function to plot one compound
plot_compound_hits <- function(
    df,
    compound_name,
    y_limits,
    zero_replacement
) {
  plot_df <- df |>
    dplyr::filter(
      Compound == compound_name
    ) |>
    dplyr::mutate(
      p.value = dplyr::if_else(
        p.value == 0,
        zero_replacement,
        p.value
      ),
      neg_log_p = -log(p.value),
      point_col = dplyr::if_else(
        neg_log_p < -log(0.05),
        "grey70",
        "black"
      )
    ) |>
    dplyr::arrange(
      Concentration
    )
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = Concentration,
      y = neg_log_p
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        colour = point_col
      )
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_hline(
      yintercept = -log(0.05),
      linetype = "dashed",
      colour = "red"
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_continuous(
      limits = y_limits
    ) +
    ggplot2::labs(
      x = "log10([uM])",
      y = "-log(p.value)",
      title = compound_name
    ) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 8,
        angle = 45,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(size = 8),
      axis.title.x = ggplot2::element_text(size = 8),
      axis.title.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 10,
        face = "bold"
      ),
      panel.grid = ggplot2::element_blank()
    )
}
# run function to plot grid of lmme_p by compound ####
# order ONLY multi-concentration compounds by smallest p-value
compound_order <- multi_conc_df |>
  dplyr::group_by(Compound) |>
  dplyr::summarise(
    min_p = min(
      p.value,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    min_p
  ) |>
  dplyr::pull(
    Compound
  )
# make one plot per multi-concentration compound
compound_plots <- purrr::map(
  compound_order,
  ~ plot_compound_hits(
    multi_conc_df,
    .x,
    y_limits = y_limits,
    zero_replacement = zero_replacement
  )
)
# name plots by compound
names(compound_plots) <- compound_order
# add individual plots to existing plots list
plots <- c(
  plots,
  compound_plots
)
# align plots
aligned_plots <- cowplot::align_plots(
  plotlist = compound_plots,
  align = "hv",
  axis = "tblr"
)
# make grid
plots$p_grid <- cowplot::plot_grid(
  plotlist = aligned_plots,
  ncol = grid_width
)
plots$p_grid
# create function to plot grid of md by compound ####
# identify conditions for compounds with more than one unique concentration
multi_conc_conditions <- data$md_lmme_p |>
  dplyr::filter(
    !is.na(Compound),
    !is.na(Concentration)
  ) |>
  dplyr::group_by(Compound) |>
  dplyr::filter(
    dplyr::n_distinct(Concentration) > 1
  ) |>
  dplyr::ungroup() |>
  dplyr::distinct(Condition) |>
  dplyr::pull(Condition)
# lookup table for Condition -> Compound + Concentration
condition_lookup <- data$md_lmme_p |>
  dplyr::distinct(
    Condition,
    Compound,
    Concentration
  )
# keep only compounds with multiple concentrations
multi_conc_md_df <- data$md_chi_mean |>
  dplyr::filter(
    Condition %in% multi_conc_conditions
  ) |>
  dplyr::left_join(
    condition_lookup,
    by = "Condition"
  )
# create function to plot one compound
plot_md_compound_hits <- function(df, compound_name) {
  plot_df <- df |>
    dplyr::filter(
      Compound == compound_name
    ) |>
    dplyr::arrange(
      Concentration
    ) |>
    dplyr::mutate(
      # grey if below cutoff, black if above cutoff
      point_col = dplyr::if_else(
        md2_mean < md_cutoff,
        "grey60",
        "black"
      )
    )
  # determine y-axis maximum for this compound
  # include SEM so error bars are not cut off
  max_md <- max(
    plot_df$md2_mean +
      dplyr::coalesce(plot_df$md2_sem, 0),
    na.rm = TRUE
  )
  y_max <- dplyr::case_when(
    max_md <= 250  ~ 250,
    max_md <= 500  ~ 500,
    max_md <= 1000 ~ 1000,
    max_md <= 2000 ~ 2000,
    max_md <= 4000 ~ 4000,
    TRUE           ~ 5000
  )
  # plot
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = Concentration,
      y = md2_mean
    )
  ) +
    # mean MD2
    ggplot2::geom_point(
      ggplot2::aes(
        colour = point_col
      ),
      size = 1.5
    ) +
    # SEM
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = md2_mean - md2_sem,
        ymax = md2_mean + md2_sem,
        colour = point_col
      ),
      width = 0.1,
      linewidth = 0.4,
      na.rm = TRUE
    ) +
    # use actual colour values stored in point_col
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_continuous(
      limits = c(0, y_max)
    ) +
    # chi-squared cutoff
    ggplot2::geom_hline(
      yintercept = md_cutoff,
      linetype = "dashed",
      colour = "red"
    ) +
    ggplot2::labs(
      x = "log10([uM])",
      y = "Mean Mahalanobis d²",
      title = compound_name
    ) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 8,
        angle = 45,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(size = 8),
      axis.title.x = ggplot2::element_text(size = 8),
      axis.title.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 10,
        face = "bold"
      ),
      panel.grid = ggplot2::element_blank()
    )
}
# run function to plot grid of md by compound ####
compound_order <- multi_conc_md_df |>
  dplyr::group_by(Compound) |>
  dplyr::summarise(
    max_md = max(md2_mean, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(max_md)
  ) |>
  dplyr::pull(Compound)
# use same compound order as p-value grid
compound_md_plots <- purrr::map(
  compound_order,
  ~ plot_md_compound_hits(
    multi_conc_md_df,
    .x
  )
)
# name plots by compound
names(compound_md_plots) <- compound_order
# add individual plots to existing plots list
plots <- c(
  plots,
  compound_md_plots
)
# align plots
aligned_md_plots <- cowplot::align_plots(
  plotlist = compound_md_plots,
  align = "hv",
  axis = "tblr"
)
# make grid
plots$md_grid <- cowplot::plot_grid(
  plotlist = aligned_md_plots,
  ncol = grid_width
)
plots$md_grid
# save data ####
write.csv(plot_df,
          paste(
            "outputs/data/", file_name, "_", integrate_state, "_", redu_state, "_md_lmme_hits.csv", sep = "")
)
# save plots ####
ggsave(
  filename = paste0("outputs/figures/", file_name, "_md_p_dr_grid.pdf"),
  plot = plots$p_grid,
  width = plot_width,
  height = plot_height
)
ggsave(
  filename = paste0("outputs/figures/", file_name, "_md_dr_grid.pdf"),
  plot = plots$md_grid,
  width = plot_width,
  height = plot_height
)
ggsave(
  filename = paste0("outputs/figures/", file_name, "_lmme_hits.pdf"),
  plot = plots$hits_lmme,
  width = 3.5,
  height = 6
)
ggsave(
  filename = paste0("outputs/figures/", file_name, "_chi_hits.pdf"),
  plot = plots$hits_chi,
  width = 3.5,
  height = 20
)
rm(list = ls())
