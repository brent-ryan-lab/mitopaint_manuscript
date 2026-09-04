# Title: profile heatmap mitopaint vis (mean per well) v1
# Step: 8.1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 03-09-2026

# load packages ####
library(data.table)
library(ComplexHeatmap)
library(tidyverse)
# set variables ####
file_name <- "mPaintSpace2_N1_N2_N3"
redu_state <- "redu"
integrate_state <- "integrated"
plot_cond <- c("DMSO_0", "CCCP_20", "ROT_3", "Nigericin_3", "Oligomycin_10",
               "Rapamycin_10", "Valinomycin_5", "Cyclosporin A_10", "CQ_10",
               "BAM15_10", "MitoQ_3", "Nocodazole_10", "Cytochalasin D_5")
plot_ann_meta <- c("Compound", "Batch")
plot_ann_cols <- list(
  Compound = c(
    "DMSO" = lighten("#440154FF", amount = 0.3),
    "CCCP" = lighten("#238A8DFF", amount = 0.3),
    "ROT"  = lighten("#FDE725FF", amount = 0.3),
    "Nigericin" = scales::hue_pal()(10)[1],
    "Oligomycin" = scales::hue_pal()(10)[2],
    "Rapamycin" = scales::hue_pal()(10)[3],
    "Valinomycin"  = scales::hue_pal()(10)[4],
    "Cyclosporin A" = scales::hue_pal()(10)[5],
    "CQ" = scales::hue_pal()(10)[6], 
    "BAM15" = scales::hue_pal()(10)[7],  
    "MitoQ" = scales::hue_pal()(10)[8],   
    "Nocodazole" = scales::hue_pal()(10)[9],  
    "Cytochalasin D" = scales::hue_pal()(10)[10]
  ),
  Batch = c(
    "N1A" = lighten(viridis(6), amount = 0.3)[1],
    "N1B" = lighten(viridis(6), amount = 0.3)[2],
    "N2A" = lighten(viridis(6), amount = 0.3)[3],
    "N2B" = lighten(viridis(6), amount = 0.3)[4],
    "N3A" = lighten(viridis(6), amount = 0.3)[5],
    "N3B" = lighten(viridis(6), amount = 0.3)[6]
  )
)
plot_height <- 4.5
plot_width <- 8
visible_legend <- TRUE
# load data ####
# create function to load data ####
load_data <- function(file_name,
                      integrate_state,
                      redu_state) {
  # load data
  df <- as.data.frame(
    fread(paste0("data/processed/",file_name,"_data_",integrate_state,
                 "_",redu_state,".csv"),
          header = TRUE))
  rownames(df) <- df$V1
  df$V1 <- NULL
  # load metadata
  meta <- as.data.frame(
    fread(paste0("data/processed/",file_name,"_",integrate_state,
                 "_",redu_state,"_dimred_meta.csv"),
          header = TRUE))
  rownames(meta) <- meta$V1
  meta$V1 <- NULL
  # keep rows aligned between meta and umap embeddings
  meta <- meta[rownames(df), , drop = FALSE]
  return(list(
    df = df,
    meta = meta))
}
# load data ####
data <- load_data(
  file_name = file_name,
  integrate_state = integrate_state,
  redu_state = redu_state
)
# average data ####
# average technical replicates by Batch and Condition #
# average technical replicates by Batch and Condition ####
avg <- data$meta |>
  dplyr::select(
    Compound,
    Concentration,
    Condition,
    Batch
  ) |>
  cbind(data$df) |>
  dplyr::group_by(
    Batch,
    Condition,
    Compound,
    Concentration
  ) |>
  dplyr::summarise(
    dplyr::across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    ID = paste(Condition, Batch, sep = "_")
  )
# split into averaged data and metadata
avg <- list(
  df = avg |>
    tibble::column_to_rownames("ID") |>
    dplyr::select(
      -Batch,
      -Condition,
      -Compound,
      -Concentration
    ),
  meta = avg |>
    dplyr::select(
      ID,
      Compound,
      Concentration,
      Condition,
      Batch
    ) |>
    tibble::column_to_rownames("ID")
)
# organise data to plot
# meta
plot_meta <- avg$meta |>
  dplyr::filter(Condition %in% plot_cond) |>
  dplyr::mutate(
    Condition = factor(Condition, levels = plot_cond)
  ) |>
  dplyr::arrange(Condition, Batch)
# matrix
plot_mat <- t(as.matrix(avg$df[rownames(plot_meta), ]))
# annotation
plot_ann <- plot_meta[, plot_ann_meta]
plot_ann <- HeatmapAnnotation(
  df = plot_ann,
  col = plot_ann_cols,
  show_legend = visible_legend,
  annotation_name_gp = gpar(fontsize = 9) 
)
# plot
plot <-  Heatmap(
  plot_mat,
  name = "Robust Z",
  col = colorRamp2(c(-30, 0, 30),
                   c("blue", "white", "red")),
  top_annotation = plot_ann,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_title = "Features",
  rect_gp = gpar(col = NA),
  show_row_dend = FALSE
)
draw(plot)
# save plot ####
pdf(paste("outputs/figures/",
          file_name,
          "_",
          integrate_state,
          "_",
          redu_state,
          "_avg_fingerprint.pdf", sep = ""),
    width = plot_width ,
    height = plot_height,
    useDingbats = FALSE)
draw(plot)
dev.off()
