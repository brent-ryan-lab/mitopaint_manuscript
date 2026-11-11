# Title: classic readouts mitopaint data (mean per well) v2
# Step: 9
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 10-11-2026

# load packages ####
library(data.table)
library(colorspace)
library(stringr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(purrr)
library(lme4)
library(lmerTest)
library(emmeans)
# set file variables ####
# make sure order of plot_cond and plot_lab are the same
# make sure order of plot_feats and y_lab are the same
batches_info <- list(
  N1 = list(
    file_name = "SF240627_mPaintDR2_N2_Classic",
    dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9"),
    batch_name = "N1"
  ),
  N2 = list(
    file_name = "SF240704_mPaintDR2_N3_Classic",
    dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9"),
    batch_name = "N2"
  ),
  N3 = list(
    file_name = "SF240711_mPaintDR2_N4_Classic",
    dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9"),
    batch_name = "N3"
  )
)
nuc_count <- "Non-border cells Selected - Number of Objects"
meta_cols <- c("Row",
               "Column",
               "Compound",	
               "Concentration")
file_name <- "mPaintDR2_Classic_N2_N3_N4"
plot_cond <- c("DMSO_0", "CCCP_30", "ROT_10")
plot_lab <- c("DMSO", "CCCP", "ROT")
plot_feats <- c("Intensity Cytoplasm CellRox Mean",
                "Intensity Cytoplasm TMRM Mean",
                "Mitochondria Selected Ratio Width to Length",
                "Number of Selected Spots/ Selected Cell")
y_lab <- c("Cytoplasm ROS Intensity (a.u.)",
           "Cytoplasm MMP Intensity (a.u.)",
           "Mitochondria Width:Length (a.u.)",
           "Mitophagy Spots (a.u.)")
x_lab <- "Compound"
pastel_cols <- lighten(c("#440154FF", "#238A8DFF", "#FDE725FF"), amount = 0.3)
point_size <- 3
size_annot <- 6
size_axis <- 12
plot_width <- 3
plot_height <- 4
min_nuc <- 100
# pair up plot_feats and y_lab ####
# pairs link the shorthand (y_lab) for each feature to the full feature name
y_lab_lookup <- setNames(
  y_lab,
  plot_feats
)
# create function to load data ####
load_data <- function(file_name, batch_name, dmso_wells) {
  # load classic data as df
  df <- as.data.frame(
    fread(
      paste("data/raw/", file_name, ".csv", sep = ""),
      header = TRUE, skip = "Row"
    )
  )
  # keep rownames as WELL_BATCH
  rownames(df) <- df$V1
  df$V1 <- NULL
  # separate metadata 
  meta <- df %>%
    select(any_of(meta_cols))
  df <- df %>%
    select(-any_of(meta_cols))
  # populate additional metadata
  meta$Well  <- paste(meta$Column, meta$Row, sep = "_")
  meta$Batch <- batch_name
  meta$ID    <- paste(meta$Well, meta$Batch, sep = "_")
  meta$Condition <- paste(meta$Compound, meta$Concentration, sep = "_")
  rownames(meta) <- meta$ID
  rownames(df)   <- meta$ID
  # clean column names
  names(df) <- names(df) |>
    str_remove("^Non-border cells Selected - ") |>
    str_remove(" - Mean per Well$") |>
    str_replace("mKeima ph4(?:/TMRM)?", "mt-Keima pH4") |>
    str_replace("mKeima ph7", "mt-Keima pH7") |>
    str_trim(side = "right")
  # convert all features to numeric
  df[] <- lapply(df, as.numeric)
  # remove any low cell count
  df <- df[df$`Number of Objects` > min_nuc, , drop = FALSE]
  # remove any columns with NA
  df <- df[, colSums(is.na(df)) == 0, drop = FALSE]
  # remove any rows with NA
  df <- df[rowSums(is.na(df)) == 0, , drop = FALSE]
  # keep meta and df aligned
  meta <- meta[rownames(df), , drop = FALSE]
  # return data and metadata as a list
  return(list(data = df,
              meta = meta,
              file_name = file_name,
              batch_name = batch_name,
              dmso_wells = dmso_wells)
  )
}
# run function to load data ####
# loop through loading all Ns given in batch_info
batches <- purrr::map(
  batches_info,
  function(batch_info) {
    load_data(
      file_name = batch_info$file_name,
      batch_name = batch_info$batch_name,
      dmso_wells = batch_info$dmso_wells
    )
  }
)
# count how many batches ####
# number of independent batches
n_batches <- length(batches)
# whether statistical testing should be performed
run_stats <- n_batches > 1
# create function to fold change normalise data to dmso ####
norm_data <- function(df, meta, dmso_wells) {
  # df_norm is where norm values will be stored
  # to start, it is just populated with df data but will soon be replaced
  df_norm <- df
  # subset dmso wells
  meta_dmso <- meta[meta$Well %in% dmso_wells, , drop = FALSE]
  df_dmso <- df[rownames(meta_dmso), , drop = FALSE]
  # for() loop normalises to dmso mean for each feature col
  for(col in colnames(df)) {
    # calculate mean of DMSO wells
    dmso_mean <- mean(df_dmso[[col]], na.rm = TRUE)
    # norm = feature value/mean of dmso
    # if statement is to avoid divide by 0 problems, if dmso_mean = 0
    if (dmso_mean == 0 || is.na(dmso_mean)) {
      df_norm[[col]] <- NA_real_
    } else {
      df_norm[[col]] <- df[[col]] / dmso_mean
    }
  }
  return(df_norm)
}
# run function to fold change normalise data to dmso ####
# loop through fold change normalising all Ns given in batches
batches <- purrr::map(
  batches,
  function(batch_obj) {
    batch_obj$df_norm <- norm_data(
      df = batch_obj$data,
      meta = batch_obj$meta,
      dmso_wells = batch_obj$dmso_wells
    )
    batch_obj
  }
)
# combine Ns in all_norm_filt ####
# creates single dataframe col binding data and meta 
# then row binding all_norm_filt N
# all_norm_filt is filtered to features and conditions of interest (plot_feats, plot_cond)
all_norm_filt <- purrr::map(
  batches,
  function(batch_obj) {
    meta <- batch_obj$meta
    data <- batch_obj$df_norm[, plot_feats, drop = FALSE] |>
      rownames_to_column("ID")
    left_join(meta, data, by = "ID")
  }
) |>
  bind_rows() |>
  filter(Condition %in% plot_cond)
all_norm_filt <- all_norm_filt |>
  column_to_rownames("ID")
# all_norm is not filtered for conditions
all_norm <- purrr::map(
  batches,
  function(batch_obj) {
    meta <- batch_obj$meta
    data <- batch_obj$df_norm[, plot_feats, drop = FALSE] |>
      rownames_to_column("ID")
    left_join(meta, data, by = "ID")
  }
) |>
  bind_rows() 
all_norm <- all_norm |>
  column_to_rownames("ID")
# all is not filtered for conditions or fold change normalised
all <- purrr::map(
  batches,
  function(batch_obj) {
    meta <- batch_obj$meta
    data <- batch_obj$data |>
      rownames_to_column("ID")
    left_join(meta, data, by = "ID")
  }
) |>
  bind_rows() 
all <- all |>
  column_to_rownames("ID")
# calculate batch level means in all_mean ####
# grouping by Batch and Condition, summarises the mean
# group means are the coloured dots in the final plot
all_mean <- all_norm_filt |>
  group_by(Batch, Condition) |>
  summarise(
    across(
      all_of(plot_feats),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )
# calculate mixed effect models for each plot_feats ####
# factor reorder makes sure levels are in the order of plot_cond
# control (DMSO_0) must be the first condition given in plot_cond
all_norm_filt$Condition <- factor(
  all_norm_filt$Condition,
  levels = plot_cond
)
all_mean$Condition <- factor(
  all_mean$Condition,
  levels = plot_cond
)
# model assumes fixed effect of condition, and random effect of batch
# this asks: what is the effect of Condition
# while accounting for variability between Batches?
# 1|Batch, assumes that each batch has its own baseline/ vertical offset
# calculate statistics only with multiple independent batches ####
if (run_stats) {
  message(
    n_batches,
    " batches detected: running linear mixed-effects models."
  )
  # model assumes fixed effect of condition
  # and random intercept for batch
  models <- lapply(
    plot_feats,
    function(feat) {
      lmerTest::lmer(
        as.formula(
          paste0(
            "`",
            feat,
            "` ~ Condition + (1|Batch)"
          )
        ),
        data = all_norm_filt
      )
    }
  )
  # name models by feature
  names(models) <- plot_feats
  # collate p-values
  emm_df <- purrr::imap_dfr(
    models,
    function(model, feat) {
      emmeans::emmeans(
        model,
        pairwise ~ Condition,
        adjust = "bonferroni"
      )$contrasts |>
        as.data.frame() |>
        dplyr::mutate(
          feature = feat,
          .before = 1
        )
    }
  )
} else {
  message(
    "Only one batch detected: skipping statistical testing."
  )
  # placeholders so downstream code still works
  models <- NULL
  emm_df <- NULL
}
# save summary of values plotted ####
# summary stats are the bars and SEM in the final plot
# save summary of values plotted ####
# statistics are calculated across independent batch means
summary_df <- purrr::map_dfr(
  plot_feats,
  function(feature) {
    all_mean |>
      dplyr::group_by(Condition) |>
      dplyr::summarise(
        # number of independent batches
        n = dplyr::n(),
        # mean across batch means
        mean = mean(
          .data[[feature]],
          na.rm = TRUE
        ),
        # SD across batches
        sd = if (dplyr::n() > 1) {
          sd(
            .data[[feature]],
            na.rm = TRUE
          )
        } else {
          NA_real_
        },
        .groups = "drop"
      ) |>
      dplyr::mutate(
        # SEM only meaningful with >1 independent batch
        sem = dplyr::if_else(
          n > 1,
          sd / sqrt(n),
          NA_real_
        ),
        feature = feature,
        .before = 1
      )
  }
)
# create a function to make barplot of plot_feat ####
plot_feature <- function(
    feature,
    model = NULL,
    all_norm_filt,
    all_mean,
    summary_df,
    pastel_cols,
    run_stats = TRUE
) {
  # batch means
  batch_means <- all_mean |>
    dplyr::transmute(
      Batch,
      Condition,
      batch_mean = .data[[feature]]
    )
  # summary for bars
  plot_summary <- summary_df |>
    dplyr::filter(
      .data$feature == .env$feature
    )
  # determine y maximum
  y_max <- max(
    c(
      all_norm_filt[[feature]],
      plot_summary$mean +
        dplyr::coalesce(plot_summary$sem, 0)
    ),
    na.rm = TRUE
  )
  # base plot
  p <- ggplot(
    all_norm_filt,
    aes(
      x = Condition,
      y = .data[[feature]],
      fill = Condition
    )
  ) +
    # bar = mean across independent batch means
    geom_col(
      data = plot_summary,
      aes(
        x = Condition,
        y = mean,
        fill = Condition
      ),
      inherit.aes = FALSE,
      colour = "black",
      width = 0.8
    ) +
    # faint points = technical replicate wells
    geom_point(
      aes(
        group = Batch
      ),
      position = position_jitter(
        width = 0.2
      ),
      shape = 16,
      colour = "black",
      size = point_size,
      alpha = 0.1
    ) +
    # solid points = mean for each independent batch
    geom_point(
      data = batch_means,
      aes(
        x = Condition,
        y = batch_mean,
        fill = Condition
      ),
      inherit.aes = FALSE,
      position = if (n_batches > 1) {
        position_jitter(width = 0.1, height = 0)
      } else {
        position_identity()
      },
      shape = 21,
      size = point_size,
      colour = "black",
      stroke = 0.6
    )
  # SEM only if there are multiple batches
  if (n_batches > 1) {
    p <- p +
      geom_errorbar(
        data = plot_summary,
        aes(
          x = Condition,
          ymin = mean - sem,
          ymax = mean + sem
        ),
        inherit.aes = FALSE,
        width = 0.5,
        colour = "black",
        na.rm = TRUE
      )
  }
  # formatting
  p <- p +
    scale_fill_manual(
      values = pastel_cols
    ) +
    scale_x_discrete(
      labels = setNames(
        plot_lab,
        plot_cond
      )
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(
        0,
        if (run_stats) {
          y_max * 1.25
        } else {
          y_max * 1.05
        }
      )
    ) +
    labs(
      y = y_lab_lookup[feature],
      x = x_lab
    ) +
    theme_pubr() +
    theme(
      axis.text = element_text(
        size = size_axis
      ),
      axis.title = element_text(
        size = size_axis
      ),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  # statistics / annotations
  # only with multiple independent batches
  if (run_stats && !is.null(model)) {
    # estimated marginal means
    emm <- emmeans::emmeans(
      model,
      ~ Condition
    )
    # compare each condition against control
    contrast_df <- emmeans::contrast(
      emm,
      method = "trt.vs.ctrl",
      ref = 1,
      adjust = "bonferroni"
    ) |>
      summary() |>
      as.data.frame()
    # significance labels
    contrast_df$label <- cut(
      contrast_df$p.value,
      breaks = c(
        -Inf,
        0.001,
        0.01,
        0.05,
        Inf
      ),
      labels = c(
        "***",
        "**",
        "*",
        "ns"
      )
    )
    # >3 conditions:
    # stars directly above corresponding bars
    if (length(plot_cond) > 3) {
      contrast_df$bar_x <- seq_along(
        plot_cond
      )[-1]
      contrast_df$star_y <- y_max * 1.08
      p <- p +
        geom_text(
          data = subset(
            contrast_df,
            label != "ns"
          ),
          aes(
            x = bar_x,
            y = star_y,
            label = label
          ),
          inherit.aes = FALSE,
          size = size_annot,
          vjust = 0
        )
    } else {
      # <=3 conditions:
      # bracket-style annotations
      annotation_df <- data.frame(
        xmin = 1,
        xmax = seq_along(plot_cond)[-1],
        y = y_max *
          (
            0.95 +
              0.1 *
              seq_along(plot_cond[-1])
          ),
        label = contrast_df$label
      )
      p <- p +
        geom_segment(
          data = annotation_df,
          aes(
            x = xmin,
            xend = xmax,
            y = y,
            yend = y
          ),
          inherit.aes = FALSE
        ) +
        geom_text(
          data = subset(
            annotation_df,
            label != "ns"
          ),
          aes(
            x = (xmin + xmax) / 2,
            y = y * 1.04,
            label = label
          ),
          inherit.aes = FALSE,
          size = size_annot
        ) +
        geom_text(
          data = subset(
            annotation_df,
            label == "ns"
          ),
          aes(
            x = (xmin + xmax) / 2,
            y = y * 1.05,
            label = label
          ),
          inherit.aes = FALSE,
          size = size_annot * 0.6
        )
    }
  }
  return(p)
}
# run function to make barplot of plot_feat ####
# loop making plots for each feature defined in plot_feats
# run function to make barplot of plot_feat ####
plots <- purrr::map(
  plot_feats,
  function(feat) {
    # retrieve model only when statistics were run
    model <- if (run_stats) {
      models[[feat]]
    } else {
      NULL
    }
    plot_feature(
      feature = feat,
      model = model,
      all_norm_filt = all_norm_filt,
      all_mean = all_mean,
      summary_df = summary_df,
      pastel_cols = pastel_cols,
      run_stats = run_stats
    )
  }
)
# name plots by feature
names(plots) <- plot_feats
plots
# align axis in each plots, so that gird is not misaligned
aligned_plots <- cowplot::align_plots(
  plotlist = plots,
  align = "hv",
  axis = "tblr"
)
names(aligned_plots) <- plot_feats
# save data ####
write.csv(
  all_norm,
  paste("data/processed/", file_name, "_norm.csv", sep = "")
)
write.csv(
  all_norm_filt,
  paste("data/processed/", file_name, "_norm_filt.csv", sep = "")
)
write.csv(
  all,
  paste("data/processed/", file_name, ".csv", sep = "")
)
write.csv(
  summary_df,
  paste("outputs/data/", file_name, "_final_means.csv", sep = "")
)
write.csv(
  all_mean,
  paste("outputs/data/", file_name, "_batch_means.csv", sep = "")
)
# save statistical results only when statistics were run
if (run_stats) {
  write.csv(
    emm_df,
    paste0(
      "outputs/data/",
      file_name,
      "_pvals.csv"
    )
  )
}
# save plots ####
# save plots using shorthand feature name
# remove (a.u.) if in the shorthand to avoid code breaking file names
purrr::iwalk(
  aligned_plots,
  function(plot, feat) {
    file_lab <- y_lab_lookup[feat] |>
      stringr::str_remove("\\s*\\(a\\.u\\.\\)") |>
      stringr::str_replace_all("[^[:alnum:]]+", "_")
    ggsave(
      paste0(
        "outputs/figures/",
        file_name,
        "_",
        file_lab,
        "_barplot.pdf"
      ),
      plot = plot,
      width = plot_width,
      height = plot_height,
      units = "in",
      dpi = 300
    )
  }
)
rm(list = ls())

