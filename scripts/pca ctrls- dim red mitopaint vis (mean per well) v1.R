# Title: pca ctrls- dim red mitopaint vis (mean per well) v1
# Step: 6.1.6
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
pos_ctrls <- c(
  "CCCP", "ROT", "Nocodazole", "47", "923", "954", "847", "468", "441","647",
  "147","481","859","118","184","683","1091","1010","344","313","825","1099",
  "1008","1067","516","489","AntimycinA","MLN4924","Chloroquine")
plot_width <- 6
plot_height <- 5
grid_plot_width <- 10
grid_width <- 12
grid_height <- 4
# create function to load data ####
load_data <- function(file_name,
                      integrate_state,
                      redu_state) {
  if (avg_profile == TRUE) {
    # load pca embeddings
    pca_embeddings <- as.data.frame(
      fread(paste0("data/processed/",file_name,"_",integrate_state,
                   "_",redu_state,"_pca_embeddings.csv"),
            header = TRUE))
    rownames(pca_embeddings) <- pca_embeddings$V1
    pca_embeddings$V1 <- NULL 
    # load pca var
    pca_var <- as.data.frame(
      fread(paste0("data/processed/",file_name,"_",integrate_state,
                   "_",redu_state,"_pca_var.csv"),
            header = TRUE))
    rownames(pca_var) <- pca_var$V1
    pca_var$V1 <- NULL
    # load metadata
    meta <- as.data.frame(
      fread(paste0("data/processed/",file_name,"_",integrate_state,
                   "_",redu_state,"_avg_dimred_meta.csv"),
            header = TRUE))
    rownames(meta) <- meta$V1
    meta$V1 <- NULL
    # keep rows aligned between meta and pca embeddings
    meta <- meta[rownames(pca_embeddings), , drop = FALSE]
  } else {
    # load pca embeddings
    pca_embeddings <- as.data.frame(
      fread(paste0("data/processed/",file_name,"_",integrate_state,
                   "_",redu_state,"_pca_embeddings.csv"),
            header = TRUE))
    rownames(pca_embeddings) <- pca_embeddings$V1
    pca_embeddings$V1 <- NULL
    # load pca var
    pca_var <- as.data.frame(
      fread(paste0("data/processed/",file_name,"_",integrate_state,
                   "_",redu_state,"_pca_var.csv"),
            header = TRUE))
    rownames(pca_var) <- pca_var$V1
    pca_var$V1 <- NULL
    # load metadata
    meta <- as.data.frame(
      fread(paste0("data/processed/",file_name,"_",integrate_state,
                   "_",redu_state,"_dimred_meta.csv"),
            header = TRUE))
    rownames(meta) <- meta$V1
    meta$V1 <- NULL
    # keep rows aligned between meta and pca embeddings
    meta <- meta[rownames(pca_embeddings), , drop = FALSE]
  }
  return(list(
    pca_embeddings = pca_embeddings,
    pca_var = pca_var,
    meta = meta))
}
# load data ####
data <- load_data(
  file_name = file_name,
  integrate_state = integrate_state,
  redu_state = redu_state
)
x_lab <- paste0("PC_1 (",round(data$pca_var$Percent_Variance[1],2),"%)")
y_lab <-paste0("PC_2 (",round(data$pca_var$Percent_Variance[2],2),"%)")
# make plot_df ####

# add Condition_label column to meta
data$meta <- data$meta |>
  dplyr::mutate(
    Condition_label = dplyr::case_when(
      Compound %in% pos_ctrls ~ as.character(Compound),
      Compound == "DMSO" ~ "DMSO",
      TRUE ~ "Library"
    )
  )

plot_df <- cbind(
  data$meta,
  data$pca_embeddings
)

# identify PCA_NN levels represented by positive controls
pos_ctrl_pca_nn <- plot_df |>
  dplyr::filter(
    Compound %in% pos_ctrls
  ) |>
  dplyr::pull(PCA_NN) |>
  unique() |>
  sort()

# assign one colour to each PCA_NN
pca_nn_cols <- setNames(
  scales::hue_pal()(length(pos_ctrl_pca_nn)),
  paste0("PCA_NN ", pos_ctrl_pca_nn)
)

# create plotting colour group
plot_df <- plot_df |>
  dplyr::mutate(
    
    # DMSO and Library remain grey;
    # positive controls are grouped by PCA_NN
    PCA_NN_highlight = dplyr::case_when(
      Compound == "DMSO" ~ "DMSO",
      
      Compound %in% pos_ctrls ~ paste0(
        "PCA_NN ",
        PCA_NN
      ),
      
      TRUE ~ "Library"
    ),
    
    # plotting order:
    # DMSO first -> Library -> positive controls last
    plot_order = dplyr::case_when(
      Compound == "DMSO" ~ 1,
      !Compound %in% pos_ctrls ~ 2,
      TRUE ~ 3
    )
  ) |>
  dplyr::arrange(
    plot_order
  )
# plot pca - color each ctrl ####
plots <- list()

plots$pca_pos_ctrls <- 
  ggplot(
    plot_df,
    aes(
      x = PC_1,
      y = PC_2,
      colour = Condition_label
    )
  ) +
  geom_point(size = 2) +
  labs(
    title = "PCA by Condition",
    x = x_lab,
    y = y_lab,
    colour = "Condition"
  ) +
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

plots$pca_pos_ctrls
# plot pca - colour positive controls by PCA_NN ####
plots$pca_ctrls_nn <- ggplot(
  plot_df,
  aes(
    x = PC_1,
    y = PC_2,
    colour = PCA_NN_highlight
  )
) +
  geom_point(
    size = 2
  ) +
  
  scale_colour_manual(
    values = c(
      "DMSO" = "grey40",
      "Library" = "grey60",
      pca_nn_cols
    )
  ) +
  
  labs(
    title = "PCA PCA_NN highlight",
    x = x_lab,
    y = y_lab,
    colour = "PCA_NN highlight"
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

plots$pca_ctrls_nn
# plot pca - individual highlight ctrl ####
# assign consistent colour to each positive control
pos_ctrl_cols <- setNames(
  scales::hue_pal()(length(pos_ctrls)),
  pos_ctrls
)
# function to highlight one positive control #
plot_pca_pos_ctrl <- function(
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
        x = PC_1,
        y = PC_2
      ),
      colour = "grey60",
      size = 2
    ) +
    # highlighted positive control plotted on top
    ggplot2::geom_point(
      data = fg_df,
      ggplot2::aes(
        x = PC_1,
        y = PC_2
      ),
      colour = ctrl_col,
      size = 2
    ) +
    ggplot2::labs(
      title = pos_ctrl,
      x = x_lab,
      y = y_lab
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
  ~ plot_pca_pos_ctrl(
    plot_df = plot_df,
    pos_ctrl = .x,
    pos_ctrl_cols = pos_ctrl_cols
  )
)
# name each plot
names(pos_ctrl_plots) <- paste0(
  "pca_",
  pos_ctrls
)
pos_ctrl_plots
plots <- c(
  plots,
  pos_ctrl_plots
)
# plot pca grid - individual highlight ctrl ####
# make grid of positive control highlight plots #
plots$pca_pos_ctrl_grid <- cowplot::plot_grid(
  plotlist = pos_ctrl_plots,
  ncol = grid_plot_width,
  align = "hv",
  axis = "tblr"
)
plots$pca_pos_ctrl_grid

# save plots ####
# plots to save
save_plots <- c(
  "pca_pos_ctrls",
  "pca_ctrls_nn",
  "pca_pos_ctrl_grid"
)

# save selected plots only
purrr::iwalk(
  plots[save_plots],
  function(plot, plot_name) {
    
    ggsave(
      filename = paste0(
        "outputs/figures/pca/",
        file_name,
        "_",
        plot_name,
        ".pdf"
      ),
      plot = plot,
      width = ifelse(
        plot_name == "pca_pos_ctrl_grid",
        grid_width,
        plot_width
      ),
      height = ifelse(
        plot_name == "pca_pos_ctrl_grid",
        grid_height,
        plot_height
      ),
      units = "in",
      dpi = 300
    )
  }
)
