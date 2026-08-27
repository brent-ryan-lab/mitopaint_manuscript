# Title: hits- mahalanobis distance mitopaint vis (mean per well) v1
# Step: 7.2
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 25-08-26

# load packages ####
library(data.table)
library(ggplot2)
library(ggpubr)
# set file variables ####
file_name <- "mPaintSpace2_N1_N2_N3"
md_file_name <- "mPaintSpace2_N1_N2_N3_mahal_lmme_p_PC10"
integrate_state <- "integrated"
redu_state <- "redu"
pc_use <- 10
# create function to load data ####
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
  # load md_chi_p
  md_chi_p <- as.data.frame(
    fread(
      paste("data/processed/", file_name, "_mahal_chi_p_", "PC", pc_use,".csv", sep = ""),
      header = TRUE
    )
  )
  # keep rownames as WELL_BATCH
  rownames(md_chi_p) <- md_chi_p$V1
  md_chi_p$V1 <- NULL 
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
    md_chi_p = md_chi_p,
    md_lmmme_p = md_lmmme_p
  ))
}
# load data ####
data <- load_data(file_name,
                  md_file_name,
                  classic_file_name,
                  integrate_state)
# plot all hits ####
plot_df <- data$md_lmmme_p |>
  filter(p.value < 0.05) |>
  arrange(p.value) |>
  mutate(
    Condition = factor(Condition, levels = rev(unique(Condition)))
  )
plots <- list()
plots$all_hits <- ggplot(plot_df, aes(x = -log(p.value), y = Condition)) +
  geom_point() +
  geom_vline(xintercept = -log(0.05), linetype = "dashed", colour = "red") +
  labs(x = "-log(p)",
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
# run function to plot grid of lmme_p by compound ####