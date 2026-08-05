# mitopaint manuscript

this repo is in process, and contains the Rproj. and script files necessary to perform the data analysis and generate the visualisations for the mitopaint manuscript that is in preparation.
all contents of this repo in its current state are unpublised, and for all intents and purposes not open to use by any users outside of the Brent Ryan lab group.
note that the data (raw and processed) and outputs (data and figures) are not set to update in git commits until final publication, but scripts are reliant on these files and the structure in the directory for complete functionality.

## structure
- `scripts/` – runnable scripts
- `data/` – input data
- `data/raw/` – raw data
- `data/processed/` – processed data
- `doc/` – documents
- `outputs/` – results
- `outputs/data/` – results data files
- `outputs/figures/` – results figure files

## usage
- `1.` tidy raw mitopaint data (mean per well) v1
- `2.` plate drift correct mitopaint data (mean per well) v1
  - `2.1` plate drift correct mitopaint visualisations (mean per well) v1
  - `2.2` plate drift raw mitopaint visualisations (mean per well) v1
- `3.` robust zscore norm mitopaint data (mean per well) v1
  - `3.1` robust zscore norm mitopaint visualisations (mean per well) v1
- `4.` batch integration mitopaint data (mean per well) v1
  - `4.1` batch integration mitopaint vis (mean per well) v1
- `5.` remove redundant mitopaint data (mean per well) v1
- `5.1` remove redundant mitopaint visualisations (mean per well) v1
- `6.` dimensionality reduction mitopaint data (mean per well) v1
  - `6.1` pca- dim red mitopaint vis (mean per well) v2
    - `6.1.2` pca scree- dim red mitopaint vis (mean per well) v1
    - `6.1.3` pca biplot- dim red mitopaint vis (mean per well) v1
- `9.` classic readouts mitopaint data (mean per well) v1
  - `9.1` classic readouts pca vis (mean per well) v1


#    6.1.4 pca feature loadings- dim red mitopaint vis (mean per well) v1
#  6.2 umap- dim red mitopaint vis (mean per well) v1
# 7. mahalanobis distance mitopaint data (mean per well) v1
#  7.1 mahalanobis distance mitopaint visualisations (mean per well) v1
# 8. similarity heatmap mitopaint visualisations (mean per well) v1
#   8.1 profile heatmap mitopaint visualisations (mean per well) v1
#   9.2 classic readouts pca pearson corr heatmap vis (mean per well) v1
  
`tidy raw mitopaint data (mean per well) v1`
- for mPaintDR2_N2, mPaintDR2_N3, mPaintDR2_N4 use the following variables:
<batches_info <- list(
  N1 = list(
    file_name = "SF240627_mPaintDR2_N2",
    batch_name = "N1"
  ),
  N2 = list(
    file_name = "SF240704_mPaintDR2_N3",
    batch_name = "N2"
  ),
  N3 = list(
    file_name = "SF240711_mPaintDR2_N4",
    batch_name = "N3"
  )
)
rm_cols = c("Timepoint",
            "Number of Analyzed Fields",
            "Time [s]",
            "Temperature",
            "Target Temperature",
            "CO2",	"Target CO2",
            "Nuclei - Number of Objects",
            "Non-border cells Selected - Number of Objects",
            "Non-border cells Selected - Nucleus Area [µm²] - Mean per Well",
            "Non-border cells Selected - Nucleus Roundness - Mean per Well",
            "Non-border cells Selected - Cell Area [µm²] - Mean per Well",
            "Non-border cells Selected - Cell Roundness - Mean per Well",
            "Non-border cells Selected - Intensity Cytoplasm mKeima ph7 Mean - Mean per Well",
            "Non-border cells Selected - Intensity Cytoplasm mKeima ph4/TMRM Mean - Mean per Well",
            "Non-border cells Selected - mKeima ph4/ph7 ratio - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Area [µm²] - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Roundness - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Width [µm] - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Length [µm] - Mean per Well",
            "Non-border cells Selected - mkeima ph7 mitochondria Ratio Width to Length - Mean per Well",
            "Cell Type",	
            "Cell Count"
)
nuc_count = "Non-border cells Selected - Number of Objects"
meta_cols = c("Row",
              "Column",
              "Compound",	
              "Concentration")>

`plate drift correct mitopaint data (mean per well) v1`
- for mPaintDR2_N2, mPaintDR2_N3, mPaintDR2_N4 use the following variables
<batches_info <- list(
  N1 = list(
    file_name = "SF240627_mPaintDR2_N2",
    dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")
  ),
  N2 = list(
    file_name = "SF240704_mPaintDR2_N3",
    dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")
  ),
  N3 = list(
    file_name = "SF240711_mPaintDR2_N4",
    dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")
  )
)
p_sig <- 0.05>

`plate drift correct mitopaint visualisations (mean per well) v1`            
- for mPaintDR2_N2, mPaintDR2_N3, mPaintDR2_N4 use the following variables
<batches_info <- list(
  N1 = list(file_name = "SF240627_mPaintDR2_N2",
            dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                           "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                           "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                           "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                           "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")),
  N2 = list(file_name = "SF240704_mPaintDR2_N3",
            dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                           "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                           "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                           "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                           "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")),
  N3 = list(file_name = "SF240711_mPaintDR2_N4",
            dmso_wells =  c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                            "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                            "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                            "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                            "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9"))
)
file_name <- "mPaint_DR2_N2_3_4"
pos_control <- "CCCP_30"
pastel_cols <- lighten(c("#440154FF", "#238A8DFF", "#FDE725FF"), amount = 0.3)>

`plate drift raw mitopaint visualisations (mean per well) v1`
- for SF240215_mPaintDrift_DMSO_N1 use the following variables
<batches_info <- list(
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
legend_titles <- c(
  "TMRM Intensity",
  "CellROX Intensity"
)>

`robust zscore norm mitopaint data v1`
- for mPaintDR2_N2, mPaintDR2_N3, mPaintDR2_N4 use the following variables
<batches_info <- list(
  N1 = list(
    file_name = "SF240627_mPaintDR2_N2",
    dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")
  ),
  N2 = list(
    file_name = "SF240704_mPaintDR2_N3",
    dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")
  ),
  N3 = list(
    file_name = "SF240711_mPaintDR2_N4",
    dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")
  )
)>

`robust zscore norm mitopaint visualisations (mean per well) v1`
- for mPaintDR2_N2 use the following variables 
<batches_info <- list(
  N1 = list(
    file_name = "SF240627_mPaintDR2_N2",
    dmso_wells = c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")
  )
)>
                
`batch integration mitopaint data (mean per well) v1`
- for mPaintDR2_N2, mPaintDR2_N3, and mPaintDR2_N4 use the following variables
<batches_info <- list(
  N1 = list(
    file_name = "SF240627_mPaintDR2_N2"
  ),
  N2 = list(
    file_name = "SF240704_mPaintDR2_N3"
  ),
  N3 = list(
    file_name = "SF240711_mPaintDR2_N4"
  )
)
file_name <- "mPaintDR2_N2_N3_N4"
k_weight <- 50>

`batch integration mitopaint vis (mean per well) v1`
- for mPaintDR2_N2_N3_N4 use the following variables
<file_name <- "mPaintDR2_N2_N3_N4"
redu_state <- "redu"
integrate_state <- c("integrated", "unintegrated")
pastel_cols <- lighten(c("#440154FF", "#238A8DFF", "#FDE725FF"), amount = 0.3)
n_neighbors <- 30
n_epochs <- 500>

`remove redundant mitopaint data (mean per well) v1`
- for mPaintDR2_N2_N3_N4 use the following variables
<file_name <- "mPaintDR2_N2_N3_N4"
integrate_state <- "integrated"
dataset_name <- "mPaintDR2_N2_N3_N4"
cor_thresh <- 0.95
var_tol <- 1e-12
excl_feats <- c("Nucleus", "Nuclei", "mTagBFP2")>

`remove redundant mitopaint visualisations (mean per well) v1`
- for mPaintDR2_N2_N3_N4 use the following variables
<file_name <- "mPaintDR2_N2_N3_N4"
integrate_state <- "integrated"
dataset_name <- "mPaintDR2_N2_N3_N4"
cor_thresh <- 0.95
cor_thresh_range_scree <- seq(from = 0.90, to = 1, by = 0.01)
cor_thresh_range_pca <- seq(from = 0.65, to = 1, by = 0.05)
var_tol <- 1e-12
excl_feats <- c("Nucleus", "Nuclei", "mTagBFP2")
pastel_cols <- lighten(c("#440154FF", "#238A8DFF", "#FDE725FF"), amount = 0.3)>

`dimensionality reduction mitopaint data (mean per well) v1`
- for mPaintDR2_N2_N3_N4 use the following variables
<file_name <- "mPaintDR2_N2_N3_N4"
redu_state <- "redu"
integrate_state <- "integrated"
dims_use <- 1:50
k_param <- 15
res <- 1
perplexity <- 20
max_iter <- 4000>

`pca- dim red mitopaint vis (mean per well) v2`
- for mPaintDR2_N2_N3_N4 use the following variables
<file_name <- "mPaintDR2_N2_N3_N4"
redu_state <- "redu"
integrate_state <- "integrated"
pastel_cols <- lighten(c("#440154FF", "#238A8DFF", "#FDE725FF"), amount = 0.3)>

`pca scree- dim red mitopaint vis (mean per well) v1`
- for mPaintDR2_N2_N3_N4 use the following variables
- file variables
<file_name <- "mPaintDR2_N2_N3_N4"
redu_state <- "redu"
integrate_state <- "integrated">
- plot variables
<x_lab <- "Number of PCs"
y_lab <- "Cumulative % Variance Explained"
size_axis <- 8
size_point <- 1
plot_width <- 2.5
plot_height <- 2.5
text_nudge <- 6>

`pca biplot- dim red mitopaint vis (mean per well) v1`
- for mPaintDR2_N2_N3_N4 use the following variables
- file variables
<file_name <- "mPaintDR2_N2_N3_N4"
redu_state <- "redu"
integrate_state <- "integrated">
- plot variables
<chan_cols <- c("#FFB000", "#DC267F", "#23CC86", "#FE6100")
x_lab <- paste0(
  "PC_1 (",
  round(var[1, "Percent_Variance"], 2),
  "%)"
)
y_lab <- paste0(
  "PC_2 (",
  round(var[2, "Percent_Variance"], 2),
  "%)"
)
point_size <- 2
arrow_width <- 1.2
arrowhead_size <- 0.5
label_text_size <- 2.5
label_line_width <- 0.4
size_title <- 13
size_axis <- 10
plot_width <- 7
plot_height <- 7>

`classic readouts mitopaint data (mean per well) v1`
- for mPaintDR2_N2_N3_N4 use the following variables
- file variables
<batches_info <- list(
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
nuc_count = "Non-border cells Selected - Number of Objects"
meta_cols = c("Row",
              "Column",
              "Compound",	
              "Concentration")
file_name <- "mPaint_DR2_Classic_N2_3_4">
- plot variables
<plot_cond <- c("DMSO_0", "CCCP_30", "ROT_10")
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
plot_height <- 4>

`classic readouts pca vis (mean per well) v1`
- for mPaintDR2_N2_N3_N4 use the following variables
<file_name_paint <- "mPaintDR2_N2_N3_N4"
file_name_classic <- "mPaint_DR2_Classic_N2_3_4"
redu_state <- "redu"
integrate_state <- "integrated">