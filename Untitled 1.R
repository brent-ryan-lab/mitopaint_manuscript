# function to load mitopaint data files

load_profiling_data <- function (file_name, batch_name, rm_cols, meta_cols) {
  # load df 
  df <- as.data.frame(
    fread(
      paste(
        "data/raw/", file_name, sep = ""), 
      skip = "Row", header = TRUE)
    )
  
}