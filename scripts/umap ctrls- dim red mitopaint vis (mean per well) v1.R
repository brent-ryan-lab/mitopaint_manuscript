# Title: umap ctrls- dim red mitopaint vis (mean per well) v3
# Step: 6.2.3
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
pos_ctrls <- c(
  "184", "1377", "524", "246", "325", "1010", "1086")
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
  "1086" = "Escitalopram"
)
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
  data$umap_embeddings
)

# identify UMAP_NN levels represented by positive controls
pos_ctrl_umap_nn <- plot_df |>
  dplyr::filter(
    Compound %in% pos_ctrls
  ) |>
  dplyr::pull(UMAP_NN) |>
  unique() |>
  sort()

# assign one colour to each UMAP_NN
umap_nn_cols <- setNames(
  scales::hue_pal()(length(pos_ctrl_umap_nn)),
  paste0("UMAP_NN ", pos_ctrl_umap_nn)
)

# create plotting colour group
plot_df <- plot_df |>
  dplyr::mutate(
    
    UMAP_NN_highlight = dplyr::case_when(
      Compound == "DMSO" ~ "DMSO",
      
      Compound %in% pos_ctrls ~ paste0(
        "UMAP_NN ",
        UMAP_NN
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
# plot umap - color each ctrl ####
plots <- list()

# legend labels
if (is.null(ctrls_lab)) {
  
  # use original IDs
  ctrl_legend_labels <- pos_ctrls
  
} else {
  
  # use renamed labels
  ctrl_legend_labels <- unname(
    ctrls_lab[pos_ctrls]
  )
}

plots$umap_pos_ctrls <- ggplot(
  plot_df,
  aes(
    x = umap_1,
    y = umap_2,
    colour = Condition_label
  )
) +
  geom_point(
    size = 2
  ) +
  
  labs(
    title = "UMAP by Condition",
    x = "UMAP_1",
    y = "UMAP_2",
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

plots$umap_pos_ctrls
# plot umap - color each by UMAP_NN ####
plots$umap_ctrls_nn <- ggplot(
  plot_df,
  aes(
    x = umap_1,
    y = umap_2,
    colour = UMAP_NN_highlight
  )
) +
  geom_point(
    size = 2
  ) +
  
  scale_colour_manual(
    values = c(
      "DMSO" = "grey40",
      "Library" = "grey60",
      umap_nn_cols
    )
  ) +
  
  labs(
    title = "UMAP UMAP_NN highlight",
    x = "UMAP_1",
    y = "UMAP_2",
    colour = "UMAP_NN highlight"
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

plots$umap_ctrls_nn
# plot umap - individual highlight ctrl ####

# assign consistent colour to each positive control
pos_ctrl_cols <- setNames(
  scales::hue_pal()(length(pos_ctrls)),
  pos_ctrls
)

# function to highlight one positive control
plot_umap_pos_ctrl <- function(
    plot_df,
    pos_ctrl,
    pos_ctrl_cols,
    ctrls_lab = NULL
) {
  
  # determine display title
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
      title = ctrl_title,
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
  ~ plot_umap_pos_ctrl(
    plot_df = plot_df,
    pos_ctrl = .x,
    pos_ctrl_cols = pos_ctrl_cols,
    ctrls_lab = ctrls_lab
  )
)

# name each plot using original control ID
names(pos_ctrl_plots) <- paste0(
  "umap_",
  pos_ctrls
)

pos_ctrl_plots

# add individual plots to main plots list
plots <- c(
  plots,
  pos_ctrl_plots
)
# plot umap grid - individual highlight ctrl ####
plots$umap_pos_ctrl_grid <- cowplot::plot_grid(
  plotlist = pos_ctrl_plots,
  ncol = grid_plot_width,
  align = "hv",
  axis = "tblr"
)

plots$umap_pos_ctrl_grid
# plot umap labs - color each ctrl and annotate with geom_label ####

# positive-control display labels
if (is.null(ctrls_lab)) {
  
  ctrl_display_labels <- setNames(
    pos_ctrls,
    pos_ctrls
  )
  
} else {
  
  ctrl_display_labels <- ctrls_lab[
    pos_ctrls
  ]
}

# calculate median UMAP position for each positive control
ctrl_labels <- plot_df |>
  dplyr::filter(
    Compound %in% pos_ctrls
  ) |>
  dplyr::group_by(
    Compound
  ) |>
  dplyr::summarise(
    umap_1 = median(
      umap_1,
      na.rm = TRUE
    ),
    umap_2 = median(
      umap_2,
      na.rm = TRUE
    ),
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
    umap_1,
    umap_2,
    label = ""
  )

# add actual positive-control label anchors
label_repel_df <- dplyr::bind_rows(
  label_obstacles,
  ctrl_labels |>
    dplyr::select(
      umap_1,
      umap_2,
      label
    )
)

# plot
plots$umap_pos_ctrls_label <- ggplot(
  plot_df,
  aes(
    x = umap_1,
    y = umap_2,
    colour = Condition_label
  )
) +
  
  # UMAP points
  geom_point(
    size = 2
  ) +
  
  # positive-control labels
  ggrepel::geom_label_repel(
    data = label_repel_df,
    aes(
      x = umap_1,
      y = umap_2,
      label = label
    ),
    inherit.aes = FALSE,
    
    # label appearance
    colour = "black",
    fill = "white",
    linewidth = 0.35,
    size = 3,
    fontface = "bold",
    
    # label repulsion
    box.padding = 0.5,
    point.padding = 0.3,
    force = 1,
    force_pull = 1,
    
    # segment line
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
    title = "UMAP by Condition",
    x = "UMAP_1",
    y = "UMAP_2",
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
    axis.text = element_text(
      size = 7
    ),
    axis.title = element_text(
      size = 7
    ),
    
    # hide legend on labelled plot
    legend.position = "none",
    
    panel.grid = element_blank()
  )

plots$umap_pos_ctrls_label
# save plots ####

# plots to save
save_plots <- c(
  "umap_pos_ctrls",
  "umap_ctrls_nn",
  "umap_pos_ctrl_grid",
  "umap_pos_ctrls_label"
)

purrr::iwalk(
  plots[save_plots],
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
      width = ifelse(
        plot_name == "umap_pos_ctrl_grid",
        grid_width,
        plot_width
      ),
      height = ifelse(
        plot_name == "umap_pos_ctrl_grid",
        grid_height,
        plot_height
      ),
      units = "in",
      dpi = 300
    )
  }
)
