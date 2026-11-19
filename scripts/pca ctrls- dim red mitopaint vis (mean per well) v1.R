# Title: pca ctrls- dim red mitopaint vis (mean per well) v2
# Step: 6.1.6
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 19-11-2026

# load packages ####
library(data.table)
library(ggplot2)
library(colorspace)
library(ggpubr)
library(tidyverse)
library(cowplot)
library(ggrepel)
# set variables ####
file_name <- "mPaintFDA_N1_N2_N3_N4_N5_N6_N7_N8"
redu_state <- "redu"
integrate_state <- "integrated"
avg_profile <- TRUE
pos_ctrls <- c("184", "1377", "524", "246", "325", "1010", "1086")
plot_width <- 5
plot_height <- 4
grid_plot_width <- 3
grid_width <- 5
grid_height <- 7.5
ctrls_lab <- c(
  "184" = "Olanzapine",
  "1377" = "Clozapine",
  "524" = "Valproate",
  "246" = "Mirtazapine", 
  "325" = "Aripiprazole", 
  "1010" = "Sertraline", 
  "1086" = "Escitalopram")
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
# legend labels
if (is.null(ctrls_lab)) {
  
  # use original IDs
  ctrl_legend_labels <- pos_ctrls
  
} else {
  
  # use renamed labels
  ctrl_legend_labels <- unname(ctrls_lab[pos_ctrls])
}

plots$pca_pos_ctrls <- ggplot(
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
    breaks = c(
      "DMSO",
      "Library",
      pos_ctrls
    ),
    labels = c(
      "DMSO",
      "Library",
      ctrl_legend_labels
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

# function to highlight one positive control
plot_pca_pos_ctrl <- function(
    plot_df,
    pos_ctrl,
    pos_ctrl_cols,
    ctrls_lab = NULL
) {
  
  # determine display title
  # use renamed control label if supplied,
  # otherwise use original control ID
  ctrl_title <- if (is.null(ctrls_lab)) {
    pos_ctrl
  } else {
    ctrls_lab[[pos_ctrl]]
  }
  
  # all points except highlighted positive control
  bg_df <- plot_df |>
    dplyr::filter(
      Condition_label != pos_ctrl
    )
  
  # highlighted positive control
  fg_df <- plot_df |>
    dplyr::filter(
      Condition_label == pos_ctrl
    )
  
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
    
    # labels
    ggplot2::labs(
      title = ctrl_title,
      x = x_lab,
      y = y_lab
    ) +
    
    # tidy theme
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
      
      axis.text = ggplot2::element_text(
        size = 7
      ),
      
      axis.title = ggplot2::element_text(
        size = 7
      ),
      
      legend.position = "none",
      
      panel.grid = ggplot2::element_blank()
    )
}


# make one plot for each positive control
pos_ctrl_plots <- purrr::map(
  pos_ctrls,
  ~ plot_pca_pos_ctrl(
    plot_df = plot_df,
    pos_ctrl = .x,
    pos_ctrl_cols = pos_ctrl_cols,
    ctrls_lab = ctrls_lab
  )
)


# name each plot using original control ID
names(pos_ctrl_plots) <- paste0(
  "pca_",
  pos_ctrls
)

# inspect individual plots
pos_ctrl_plots


# add individual plots to main plots list
plots <- c(
  plots,
  pos_ctrl_plots
)


# plot pca grid - individual highlight ctrl ####

# make grid of positive-control highlight plots
plots$pca_pos_ctrl_grid <- cowplot::plot_grid(
  plotlist = pos_ctrl_plots,
  ncol = grid_plot_width,
  align = "hv",
  axis = "tblr"
)

plots$pca_pos_ctrl_grid
# plot pca labs - color each ctrl and annotate with geom_label ####
# positive-control display labels
if (is.null(ctrls_lab)) {
  ctrl_display_labels <- setNames(
    pos_ctrls,
    pos_ctrls
  )
} else {
  ctrl_display_labels <- ctrls_lab[pos_ctrls]
}
# calculate median PCA position for each positive control
ctrl_labels <- plot_df |>
  dplyr::filter(
    Compound %in% pos_ctrls
  ) |>
  dplyr::group_by(
    Compound
  ) |>
  dplyr::summarise(
    PC_1 = median(PC_1, na.rm = TRUE),
    PC_2 = median(PC_2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    label = ctrl_display_labels[
      as.character(Compound)
    ]
  )
# only positive-control points act as obstacles for labels
label_obstacles <- plot_df |>
  dplyr::filter(
    Compound %in% pos_ctrls
  ) |>
  dplyr::transmute(
    PC_1,
    PC_2,
    label = ""
  )

# add actual positive-control label anchors
label_repel_df <- dplyr::bind_rows(
  label_obstacles,
  ctrl_labels |>
    dplyr::select(
      PC_1,
      PC_2,
      label
    )
)
# plot
plots$pca_pos_ctrls_label <- ggplot(
  plot_df,
  aes(
    x = PC_1,
    y = PC_2,
    colour = Condition_label
  )
) +
  
  # PCA points
  geom_point(
    size = 2
  ) +
  
  # positive-control labels
  ggrepel::geom_label_repel(
    data = label_repel_df,
    aes(
      x = PC_1,
      y = PC_2,
      label = label
    ),
    inherit.aes = FALSE,
    
    # label appearance
    colour = "black",
    fill = "white",
    linewidth = 0.35,
    size = 3,
    fontface = "bold",
    
    # keep labels away from points / other labels
    box.padding = 0.5,
    point.padding = 0.5,
    force = 5,
    force_pull = 0.5,
    
    # keep segment lines
    segment.colour = "black",
    segment.linewidth = 0.35,
    min.segment.length = 0,
    
    max.overlaps = Inf,
    show.legend = FALSE
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
    breaks = c(
      "DMSO",
      "Library",
      pos_ctrls
    ),
    labels = c(
      "DMSO",
      "Library",
      unname(ctrl_display_labels)
    ),
    drop = FALSE
  ) +
  
  labs(
    title = "PCA by Condition",
    x = x_lab,
    y = y_lab,
    colour = "Condition"
  ) +
  
  ggpubr::theme_pubr() +
  
  theme(
    aspect.ratio = 1,
    plot.title = element_text(
      hjust = 0.5,
      size = 9,
      face = "bold"
    ),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    # legend remains visible
    legend.position = "none",
    panel.grid = element_blank()
  )
plots$pca_pos_ctrls_label
# save plots ####
# plots to save
save_plots <- c(
  "pca_pos_ctrls",
  "pca_ctrls_nn",
  "pca_pos_ctrl_grid",
  "pca_pos_ctrls_label"
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
