# Title: hits- mahalanobis distance mitopaint vis (mean per well) v1
# Step: 7.2
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 07-09-26

# load packages ####
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
# set file variables ####
file_name <- "mPaintDR2_N2_N3_N4"
md_file_name <- "mPaintDR2_N2_N3_N4_mahal_lmme_p_PC10"
integrate_state <- "integrated"
redu_state <- "redu"
pc_use <- 10
grid_width <- 2
plot_width <- 5
plot_height <- 2.5
# create function to load data ####
load_data <- function(file_name,
                      md_file_name,
                      classic_file_name,
                      integrate_state) {
  # load mahalanobis distance (md)
  md <- as.data.frame(
    fread(
      paste("data/processed/", md_file_name, ".csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames as WELL_BATCH
  rownames(md) <- md$V1
  md$V1 <- NULL
  # load md_lmme_p
  md_lmmme_p <- as.data.frame(
    fread(
      paste("data/processed/", file_name, "_mahal_lmme_p_", "PC", pc_use,".csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames
  rownames(md_lmmme_p) <- md_lmmme_p$V1
  md_lmmme_p$V1 <- NULL  
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
    md = md,
    meta = meta,
    md_lmmme_p = md_lmmme_p
  ))
}
# load data ####
data <- load_data(file_name,
                  md_file_name,
                  classic_file_name,
                  integrate_state)
# plot all hits ####
# filter hits to those with p < 0.05
plot_df <- data$md_lmmme_p |>
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
# all_hits is dot plot of compounds with significant p
plots$all_hits <- ggplot(plot_df, 
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
plots$all_hits
# create function to plot grid of lmme_p by compound ####
# replace zeros once for determining limits
p_vals <- data$md_lmmme_p$p.value
min_pos_p <- min(p_vals[p_vals > 0 & is.finite(p_vals)], na.rm = TRUE)
zero_replacement <- 0.9 * min_pos_p
# shared y limits based on all p-values in md_lmmme_p
y_vals <- -log(ifelse(
  p_vals == 0,
  0.9 * min_pos_p,
  p_vals
)[is.finite(p_vals)])
# consistent y limits for all plots based on min and max p val
y_limits <- range(y_vals, na.rm = TRUE) * c(0.95, 1.05)
# create function to plot one compound
plot_compound_hits <- function(df, compound_name, y_limits, zero_replacement) {
  plot_df <- df |>
    dplyr::filter(Compound == compound_name) |>
    dplyr::mutate(
      p.value = dplyr::if_else(
        p.value == 0,
        zero_replacement,
        p.value
      ),
      neg_log_p = -log(p.value),
      point_col = dplyr::if_else(neg_log_p < -log(0.05), "grey70", "black")
    ) |>
    dplyr::arrange(Concentration)
  
  ggplot2::ggplot(plot_df, ggplot2::aes(x = Concentration, y = neg_log_p)) +
    ggplot2::geom_point(ggplot2::aes(colour = point_col)) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_hline(
      yintercept = -log(0.05),
      linetype = "dashed",
      colour = "red"
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_continuous(limits = y_limits) +
    ggplot2::labs(
      x = "log10([uM])",
      y = "-log10(p.value)",
      title = compound_name
    ) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8),
      axis.title.x = ggplot2::element_text(size = 8),
      axis.title.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold"),
      panel.grid = ggplot2::element_blank()
    )
}
# run function to plot grid of lmme_p by compound ####
# make one plot per compound #
# order plots on grid by smallest to largest p (biggest to smallest -log(p))
compound_order <- data$md_lmmme_p |>
  group_by(Compound) |>
  summarise(min_p = min(p.value, na.rm = TRUE), .groups = "drop") |>
  arrange(min_p) |>
  pull(Compound)
# loop plot all compounds
compound_plots <- purrr::map(
  compound_order,
  ~ plot_compound_hits(
    data$md_lmmme_p,
    .x,
    y_limits = y_limits,
    zero_replacement = zero_replacement
  )
)
# name plots by compound name
names(compound_plots) <- compound_order
plots <- c(plots, compound_plots)
# add them to the existing plots list #
# match size/ scaling of all plots 
aligned_plots <- align_plots(
  plotlist = compound_plots,
  align = "hv",
  axis = "tblr"
)
# make grid
plots$p_grid <- plot_grid(
  plotlist = aligned_plots,
  ncol = grid_width
)
plots$p_grid
# create function to plot grid of md by compound ####
# create function to plot one compound
plot_md_compound_hits <- function(df, compound_name) {
  plot_df <- df |>
    dplyr::filter(Compound == compound_name) |>
    dplyr::arrange(Concentration)
  # determine y-axis limits for this compound
  max_md <- max(plot_df$estimate, na.rm = TRUE)
  y_max <- dplyr::case_when(
    max_md <= 250  ~ 250,
    max_md <= 500  ~ 500,
    max_md <= 1000 ~ 1000,
    TRUE           ~ 2000
  )
  y_min <- min(plot_df$estimate - plot_df$SE, na.rm = TRUE)
  ggplot2::ggplot(plot_df, ggplot2::aes(x = Concentration, y = estimate)) +
    ggplot2::geom_point(size = 1.5, colour = "black") +
    ggplot2::geom_errorbar(
      aes(
        ymin = estimate - SE,
        ymax = estimate + SE
      ),
      width = 0.1,
      linewidth = 0.4
    ) +
    ggplot2::scale_x_log10() +
    scale_y_continuous(
      limits = c(y_min, y_max)
    ) +
    ggplot2::labs(
      x = "log10([uM])",
      y = "md",
      title = compound_name
    ) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8),
      axis.title.x = ggplot2::element_text(size = 8),
      axis.title.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold"),
      panel.grid = ggplot2::element_blank()
    )
}
# run function to plot grid of md by compound ####
# make one plot per compound #
# order plots on grid by smallest to largest md
# loop plot all compounds
compound_md_plots <- purrr::map(
  compound_order,
  ~ plot_md_compound_hits(
    data$md,
    .x
  )
)
# name plots by compound name
names(compound_md_plots) <- compound_order
plots <- c(plots, compound_md_plots)
# add them to the existing plots list #
# match size/ scaling of all plots 
aligned_md_plots <- align_plots(
  plotlist = compound_md_plots,
  align = "hv",
  axis = "tblr"
)
# make grid
plots$md_grid <- plot_grid(
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
  filename = paste0("outputs/figures/", file_name, "_md_hits_list.pdf"),
  plot = plots$all_hits,
  width = 3.5,
  height = 6
)
rm(list = ls())
