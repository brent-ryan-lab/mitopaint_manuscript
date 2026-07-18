# Title: plate drift raw mitopaint visualisations (mean per well) v1
# Step: 2.2
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 17-07-2026

# load packages ####
# set file variables ####
batches_info <- list(
  N1 = list(
    file_name = "SF240215_mPaintDrift_DMSO_N1",
    batch_name = "N1"
  )
)
meta_cols <- c("Row",
                "Column",
                "Compound",	
                "Concentration")
rm_cols <- c("Timepoint",
            "Number of Analyzed Fields",
            "Time [s]",
            "Temperature",
            "Target Temperature",
            "CO2",	"Target CO2",
            "Nuclei - Number of Objects",
            "Non-border cells Selected - Number of Objects",
            "Non-border cells Selected - Nucleus Area [¬µm¬≤] - Mean per Well",
            "Non-border cells Selected - Nucleus Roundness - Mean per Well",
            "Non-border cells Selected - Cell Area [¬µm¬≤] - Mean per Well",
            "Non-border cells Selected - Cell Roundness - Mean per Well",
            "Non-border cells Selected - Intensity Cytoplasm mKeima ph7 Mean - Mean per Well",
            "Non-border cells Selected - Intensity Cytoplasm mKeima ph4/TMRM Mean - Mean per Well",
            "Non-border cells Selected - mKeima ph4/ph7 ratio - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Area [¬µm¬≤] - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Roundness - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Width [¬µm] - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Length [¬µm] - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Ratio Width to Length - Mean per Well",
            "Cell Type",	
            "Cell Count"
)
plot_var <- c("Intensity Cytoplasm Region TMRM test Mean",
              "Intensity Cytoplasm Region CellRox Deep Red test Mean")
# create a function to load data ####
load_data <- function(file_name, batch_name, rm_cols, meta_cols) {
  # load df 
  df <- as.data.frame(
    fread(
      paste(
        "data/raw/", file_name, ".csv", sep = ""), 
      skip = "Row", header = TRUE)
  )
  # remove unwanted columns
  df <- df %>%
    select(-any_of(rm_cols))
  # separate metadata 
  meta <- df %>%
    select(any_of(meta_cols))
  df <- df %>%
    select(-any_of(meta_cols))
  # populate additional metadata
  meta$Well  <- paste(meta$Column, meta$Row, sep = "_")
  meta$Batch <- batch_name
  meta$ID    <- paste(meta$Well, meta$Batch, sep = "_")
  meta$Order <- c(1:(nrow(meta)))
  meta$Condition <- paste(meta$Compound, meta$Concentration, sep = "_")
  rownames(meta) <- meta$ID
  rownames(df)   <- meta$ID
  # clean column names
  names(df) <- names(df) |>
    str_remove("^Non-border cells Selected - ") |>
    str_remove(" - Mean per Well$") |>
    str_replace("mKeima ph4/TMRM", "mt-Keima pH4") |>
    str_replace("mKeima ph7", "mt-Keima pH7") |>
    str_trim(side = "right")
  # convert all features to numeric
  df[] <- lapply(df, as.numeric)
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
              batch_name = batch_name)
  )
}
# run function to load data
batches <- map(
  batches_info,
  function(batch_info) {
    load_data(
      file_name = batch_info$file_name,
      batch_name = batch_info$batch_name,
      rm_cols = rm_cols,
      meta_cols = meta_cols
    )
  }
)
# open data to inspect
View(batches$N1$data)
# open meta to inspect
View(batches$N1$meta)
# create a function to plot plate overview ####
plot_plate <- function(feature, batch_obj) {
  
  meta <- batch_obj$meta
  df <- batch_obj$data
  
  plot_df <- meta |>
    dplyr::mutate(
      Row_num = as.integer(Row),
      Col_num = as.integer(Column),
      value = df[[feature]]
    )
  
  plate_df <- expand.grid(
    Row_num = 1:16,
    Col_num = 1:24
  ) |>
    as_tibble() |>
    dplyr::left_join(
      plot_df |>
        dplyr::select(Row_num, Col_num, value),
      by = c("Row_num", "Col_num")
    ) |>
    dplyr::mutate(
      x = Col_num,
      y = 17 - Row_num
    )
  
  ggplot(plate_df, aes(x = x, y = y)) +
    geom_point(
      aes(fill = value),
      shape = 21,
      size = 6,
      colour = "black",
      stroke = 0.4
    ) +
    scale_fill_gradientn(
      colours = c("blue", "yellow", "red"),
      na.value = "white"
    ) +
    scale_x_continuous(
      breaks = 1:24,
      labels = 1:24,
      expand = c(0.02, 0.02)
    ) +
    scale_y_continuous(
      breaks = 1:16,
      labels = 16:1,
      expand = c(0.02, 0.02)
    ) +
    coord_equal() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10)
    )
}

plots <- purrr::imap(
  setNames(plot_var, plot_var),
  ~ plot_plate(
    feature = .x,
    batch_obj = batches$N1
  )
)
