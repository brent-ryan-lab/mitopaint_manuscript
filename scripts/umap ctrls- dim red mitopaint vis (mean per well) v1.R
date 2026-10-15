# Title: umap ctrls- dim red mitopaint vis (mean per well) v1
# Step: 6.2.3
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 15-10-2026

# load packages ####
library(data.table)
library(ggplot2)
library(colorspace)
library(ggpubr)
library(tidyverse)
library(cowplot)
# set variables ####
file_name <- "mPaintFDA_N1_N2_N3_N4_N5_N6_N7_N8"
redu_state <- "redu"
integrate_state <- "integrated"
avg_profile <- TRUE
pos_ctrls <- c("SodiumArsenite", "ROT", "Oligomycin", "Nocodazole", "MLN4924", 
               "MitoQ", "CytochalasinD", "Chloroquine", "CCCP", "AntimycinA")
plot_width <- 5
plot_height <- 4
grid_plot_width <- 3
grid_width <- 5
grid_height <- 7.5
# create function to load data ####
load_data <- function(file_name,
                      integrate_state,
                      redu_state) {
  if (avg_profile == TRUE) {
    # load umap embeddings
    umap_embeddings <- as.data.frame(
      fread(paste0("data/processed/",file_name,"_",integrate_state,
                   "_",redu_state,"_avg_umap_embeddings.csv"),
            header = TRUE))
    rownames(umap_embeddings) <- umap_embeddings$V1
    umap_embeddings$V1 <- NULL 
    # load metadata
    meta <- as.data.frame(
      fread(paste0("data/processed/",file_name,"_",integrate_state,
                   "_",redu_state,"_avg_dimred_meta.csv"),
            header = TRUE))
    rownames(meta) <- meta$V1
    meta$V1 <- NULL
    # keep rows aligned between meta and umap embeddings
    meta <- meta[rownames(umap_embeddings), , drop = FALSE]
  } else {
    # load umap embeddings
    umap_embeddings <- as.data.frame(
      fread(paste0("data/processed/",file_name,"_",integrate_state,
                   "_",redu_state,"_umap_embeddings.csv"),
            header = TRUE))
    rownames(umap_embeddings) <- umap_embeddings$V1
    umap_embeddings$V1 <- NULL
    # load metadata
    meta <- as.data.frame(
      fread(paste0("data/processed/",file_name,"_",integrate_state,
                   "_",redu_state,"_dimred_meta.csv"),
            header = TRUE))
    rownames(meta) <- meta$V1
    meta$V1 <- NULL
    # keep rows aligned between meta and umap embeddings
    meta <- meta[rownames(umap_embeddings), , drop = FALSE]
  }
  return(list(
    umap_embeddings = umap_embeddings,
    meta = meta))
}
# load data ####
data <- load_data(
  file_name = file_name,
  integrate_state = integrate_state,
  redu_state = redu_state
)
# plot umap - color each ctrl ####
# add Condition_label column to meta
data$meta <- data$meta |>
  dplyr::mutate(
    Condition_label = dplyr::case_when(
      Compound %in% pos_ctrls ~ as.character(Compound),
      Compound == "DMSO" ~ "DMSO",
      TRUE ~ "Library"
    )
  )
plot_df <- cbind(data$meta, data$umap_embeddings)
# arrange rows in plotting order
plot_df <- plot_df |>
  dplyr::mutate(
    plot_order = dplyr::case_when(
      Condition_label == "DMSO" ~ 1,
      Condition_label == "Library" ~ 2,
      TRUE ~ 3
    )
  ) |>
  dplyr::arrange(plot_order)
plots <- list()
plots$umap_pos_ctrls <- 
  ggplot(plot_df,aes(x = umap_1,y = umap_2,
                   colour = Condition_label)) +
    geom_point(size = 2) +
  labs(
    title = "UMAP by Condition",
    x = "UMAP_1",
    y = "UMAP_2",
    colour = "Condition") +
  scale_colour_manual(
    values = c(
      DMSO = "grey40",
      Library = "grey60",
      setNames(
        scales::hue_pal()(length(pos_ctrls)),
        pos_ctrls
      )
    ),
    drop = FALSE
  ) +
  theme_pubr() +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(
      hjust = 0.5,
      size = 9,
      face = "bold"
    ),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.position = "right",
    panel.grid = element_blank()
  )
plots$umap_pos_ctrls
# plot umap - individual highlight ctrl ####
# assign consistent colour to each positive control
pos_ctrl_cols <- setNames(
  scales::hue_pal()(length(pos_ctrls)),
  pos_ctrls
)
# function to highlight one positive control #
plot_umap_pos_ctrl <- function(
    plot_df,
    pos_ctrl,
    pos_ctrl_cols
) {
  # all points except highlighted positive control
  bg_df <- plot_df |>
    dplyr::filter(Condition_label != pos_ctrl)
  # highlighted positive control
  fg_df <- plot_df |>
    dplyr::filter(Condition_label == pos_ctrl)
  # colour corresponding to this positive control
  ctrl_col <- pos_ctrl_cols[[pos_ctrl]]
  # plot
  ggplot2::ggplot() +
    # all other points plotted first in grey
    ggplot2::geom_point(
      data = bg_df,
      ggplot2::aes(
        x = umap_1,
        y = umap_2
      ),
      colour = "grey60",
      size = 2
    ) +
    # highlighted positive control plotted on top
    ggplot2::geom_point(
      data = fg_df,
      ggplot2::aes(
        x = umap_1,
        y = umap_2
      ),
      colour = ctrl_col,
      size = 2
    ) +
    ggplot2::labs(
      title = pos_ctrl,
      x = "UMAP_1",
      y = "UMAP_2"
    ) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      aspect.ratio = 1,
      # title uses same colour as highlighted points
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 9,
        face = "bold",
        colour = ctrl_col
      ),
      axis.text = ggplot2::element_text(size = 7),
      axis.title = ggplot2::element_text(size = 7),
      legend.position = "none",
      panel.grid = ggplot2::element_blank()
    )
}
pos_ctrl_plots <- purrr::map(
  pos_ctrls,
  ~ plot_umap_pos_ctrl(
    plot_df = plot_df,
    pos_ctrl = .x,
    pos_ctrl_cols = pos_ctrl_cols
  )
)
# name each plot
names(pos_ctrl_plots) <- paste0(
  "umap_",
  pos_ctrls
)
pos_ctrl_plots
plots <- c(
  plots,
  pos_ctrl_plots
)
# plot umap grid - individual highlight ctrl ####
# make grid of positive control highlight plots #
plots$umap_pos_ctrl_grid <- cowplot::plot_grid(
  plotlist = pos_ctrl_plots,
  ncol = grid_plot_width,
  align = "hv",
  axis = "tblr"
)
plots$umap_pos_ctrl_grid
# save plots ####
iwalk(
  plots,
  function(plot, plot_name) {
    ggsave(
      filename = paste0(
        "outputs/figures/umap/",
        file_name,
        "_",
        plot_name,
        ".pdf"
      ),
      plot = plot,
      width = ifelse(plot_name == "umap_pos_ctrl_grid", grid_width, plot_width),
      height = ifelse(plot_name == "umap_pos_ctrl_grid", grid_height, plot_height),
      units = "in",
      dpi = 300
    )
  }
)
