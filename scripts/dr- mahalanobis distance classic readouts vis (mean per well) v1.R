# Title: # dr- mahalanobis distance classic readouts vis (mean per well) v1
# Step: 7.1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 10-08-26

# load packages ####
library(data.table)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(dplyr)
library(tidyr)
library(purrr)
library(lme4)
library(lmerTest)
library(emmeans)
library(cowplot)
library(grid)
library(stringr)
# set file variables ####
file_name <- "mPaintSpace2_N1_N2_N3"
md_file_name <- "mPaintSpace2_N1_N2_N3_mahal_lmme_p_PC10"
classic_file_name <- "mPaintSpace2_N1_N2_N3_Classic"
integrate_state <- "integrated"
plot_feats <- c("Intensity Cytoplasm CellRox Deep Red test Mean",
                "Intensity Cytoplasm TMRM test Mean",
                "mkeima ph7 mitochondria Ratio Width to Length",
                "Number of Mitophagy Spots Selected- per Cell")
y_lab <- c("Cytoplasm ROS Intensity",
           "Cytoplasm MMP Intensity",
           "Mitochondria Width:Length",
           "Mitophagy Spots")
plot_cond <- c("CCCP", "ROT", "Nigericin", "Oligomycin", "Rapamycin", "Valinomycin", "Cyclosporin A", "CQ", "BAM15", "MitoQ", "Nocodazole",  "Cytochalasin D")
pastel_cols <- c(lighten(c("#238A8DFF", "#FDE725FF"), amount = 0.3), scales::hue_pal()(10))
# create function to load data ####
y_lab_lookup <- setNames(
  y_lab,
  plot_feats
)
load_data <- function(file_name,
                      md_file_name,
                      classic_file_name,
                      integrate_state) {
  # load md
  md <- as.data.frame(
    fread(
      paste("data/processed/", md_file_name, ".csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames as WELL_BATCH
  rownames(md) <- md$V1
  md$V1 <- NULL
  # load classic readout
  classic <- as.data.frame(
    fread(
      paste("data/processed/", classic_file_name, ".csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames as WELL_BATCH
  rownames(classic) <- classic$V1
  classic$V1 <- NULL
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
  classic <- classic[rownames(classic) %in% rownames(meta),]
  # remove meta cols from classic
  classic <- classic[,!colnames(classic) %in% colnames(meta)]
  # return list
  return(list(
    md = md,
    classic = classic,
    meta = meta
  ))
}
# load data ####
data <- load_data(file_name,
                  md_file_name,
                  classic_file_name,
                  integrate_state)
# create function to dmso norm classic ####
dmso_norm <- function (data, control_condition = "DMSO_0") {
  # pull data tables
  classic <- data$classic
  meta <- data$meta
  # classic_norm populated with classic data as placeholder values
  classic_norm <- classic
  batches <- unique(meta$Batch)
  # for loop to subset data by batch, repeats individually for each batch
  for (b in batches) {
    batch_idx <- meta$Batch == b
    dmso_idx  <- batch_idx & meta$Condition == control_condition
    # for loop to subset data by classic readout, repeats individually for each readout
    for (feat in colnames(classic)) {
      dmso_mean <- mean(classic[dmso_idx, feat], na.rm = TRUE)
      if (is.na(dmso_mean) || dmso_mean == 0) {
        classic_norm[batch_idx, feat] <- NA_real_
      } else {
        classic_norm[batch_idx, feat] <- classic[batch_idx, feat] / dmso_mean
      }
    }
  }
  return(classic_norm)
}
# run function to dmso norm classic ####
data$classic_norm <- dmso_norm(data)
# stats analyse lmme model of dmso norm classic ####
classic_norm_df <- cbind(data$meta, data$classic_norm)
classic_norm_df$Batch <- factor(classic_norm_df$Batch)
classic_norm_df$Condition <- factor(classic_norm_df$Condition)
condition_lookup <- classic_norm_df |>
  distinct(Condition, Compound, Concentration)
classic_lmme_long <- map_dfr(
  colnames(data$classic_norm),
  function(feat) {
    fit <- lmer(
      as.formula(paste0("`", feat, "` ~ Condition + (1|Batch)")),
      data = classic_norm_df
    )
    contr_df <- emmeans(fit, ~ Condition) |>
      contrast(
        method = "trt.vs.ctrl",
        ref = "DMSO_0",
        adjust = "bonferroni"
      ) |>
      as.data.frame()
    contr_df |>
      separate(
        contrast,
        into = c("Condition", "ControlCondition"),
        sep = " - ",
        remove = FALSE
      ) |>
      left_join(condition_lookup, by = "Condition") |>
      mutate(feature = feat, .before = 1) |>
      select(
        feature,
        Condition,
        Compound,
        Concentration,
        ControlCondition,
        contrast,
        estimate,
        SE,
        df,
        t.ratio,
        p.value
      )
  }
)
data$classic_lmme_results <- classic_lmme_long |>
  select(
    Condition,
    Compound,
    Concentration,
    contrast,
    feature,
    p.value
  ) |>
  distinct() |>
  pivot_wider(
    names_from = feature,
    values_from = p.value
  )
rm(classic_lmme_long, condition_lookup, classic_norm_df)


# create function to dot plot pval ####
plot_dot <- function(results_df,
                             feature_col,
                             compound_name,
                             plot_title,
                             pastel_col) {
  # create plot df for specific compound being plotted
  plot_df <- results_df |>
    filter(Compound == compound_name) |>
    # if p is zero, change to pos val to avoid Inf when -log
    mutate(
      p_floor = if (any(.data[[feature_col]] > 0, na.rm = TRUE)) {
        min(.data[[feature_col]][.data[[feature_col]] > 0], na.rm = TRUE) / 10
      } else {
        .Machine$double.xmin
      },
      p_plot = if_else(
        .data[[feature_col]] <= 0,
        p_floor,
        .data[[feature_col]]
      ),
      x = log10(Concentration),
      y = -log10(p_plot)
    )
  # set y_max to avoid cutting off points
  y_max <- max(plot_df$y, na.rm = TRUE) * 1.1
  # plot data
  ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(
      fill = pastel_col,colour = "black",shape = 21,size = 3,stroke = 0.6) +
    # dashed line at p = 0.05 sig
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey70") +
    labs(
      x = "log10([uM])",
      y = "-log10(p.value)",
      title = plot_title
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(limits = c(0, y_max)) +
    theme_pubr() +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8, margin = margin(r = 6)),
      plot.margin = margin(5, 5, 5, 5),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "none",
      aspect.ratio = 1) 
}
# run function to dot plot pval ####
plots <- list()
# plot classic lmme results
for (feat in plot_feats) {
  feat_label <- y_lab_lookup[feat]
  for (comp in plot_cond) {
    plots[[paste0(comp, "_", feat)]] <- plot_dot(
      results_df = data$classic_lmme_results,
      feature_col = feat,
      compound_name = comp,
      plot_title = paste0(comp, "\n", feat_label),
      pastel_col = ifelse(comp == "CCCP", pastel_cols[1], pastel_cols[2])
    )
  }
}
# plot mpaint lmme mahalanobis distance results
plots$md_cccp <- plot_dot(
  results_df = data$md,
  feature_col = "p.value",
  compound_name = "CCCP",
  plot_title = "CCCP\nMahalanobis distance",
  pastel_col = pastel_cols[1]
)
plots$md_rot <- plot_dot(
  results_df = data$md,
  feature_col = "p.value",
  compound_name = "ROT",
  plot_title = "ROT\nMahalanobis distance",
  pastel_col = pastel_cols[2]
)
# align plots vertically and horizontally using cowplot
aligned_plots <- align_plots(
  plotlist = plots,
  align = "v",
  axis = "l"
)
for (p in aligned_plots) {
  grid.newpage()
  grid.draw(p)
}
# save data ####
write.csv(
  data$classic_lmme_results,
  paste("data/processed/", classic_file_name, "_lmme_p.csv", sep = "")
)
# save plots ####
iwalk(
  aligned_plots,
  function(plot, plot_name) {
    safe_name <- str_replace_all(plot_name, "[^[:alnum:]]+", "_")
    ggsave(
      filename = paste0("outputs/figures/", file_name, "_", safe_name, ".pdf"),
      plot = plot,
      width = 3,
      height = 3.3
    )
  }
)
rm(list = ls())
