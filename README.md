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
  - `2.1` plate drift correct mitopaint visualisations (N=3) (mean per well) v1
- `3.` robust zscore norm mitopaint data (mean per well) v1
- `4.` batch integration mitopaint data (N=3) (mean per well) v1
- `5.` remove redundant mitopaint data (N=3) (mean per well) v1
- `6.` dimensionality reduction mitopaint data (N=3) (mean per well) v1
  - `6.1` PCA: dim red mitopaint vis (N=3) (mean per well) mPaintDR2_N2_N3_N4 v1

#  3.1 robust zscore norm mitopaint visualisations (mean per well) v1
#  5.1 remove redundant mitopaint visualisations (N=3) (mean per well) v1
#  6.1 pca: dim red mitopaint vis (N=3) (mean per well) v1
#    6.1.1 pca batch integration: dim red mitopaint vis (N=3) (mean per well) v1
#    6.1.2 pca scree: dim red mitopaint vis (N=3) (mean per well) v1
#    6.1.3 pca biplot: dim red mitopaint vis (N=3) (mean per well) v1
#    6.1.4 pca feature loadings: dim red mitopaint vis (N=3) (mean per well) v1
#  6.2 umap: dim red mitopaint vis (N=3) (mean per well) v1
#    6.2.1 umap batch integration: dim red mitopaint vis (N=3) (mean per well) v1
#  6.3 tsne: dim red mitopaint vis (N=3) (mean per well) v1
#    6.3.1 tsne batch integration: dim red mitopaint vis (N=3) (mean per well) v1
# 7. mahalanobis distance mitopaint data (N=3) (mean per well) v1
#  7.1 mahalanobis distance mitopaint visualisaitons (N=3) (mean per well) v1
#  8.1 similarity heatmap mitopaint visualisations (N=3) (mean per well) v1
#  8.2 profile heatmap mitopaint visualisations (N=3) (mean per well) v1
  
`tidy raw mitopaint data (mean per well) v1`
- for mPaintDR2_N2 use the following variables:
<file_name <- "SF240627_mPaintDR2_N2"
batch_name <- "N2"
rm_cols <- c("Timepoint",
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
nuc_count <- "Non-border cells Selected - Number of Objects"
meta_cols <- c("Row",
               "Column",
               "Compound",	
               "Concentration"
)>

`tidy raw mitopaint data (mean per well) v1`
- for mPaintDR2_N3 use the following variables
<file_name <- "SF240704_mPaintDR2_N3"
batch_name <- "N3"
rm_cols <- c("Timepoint",
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
nuc_count <- "Non-border cells Selected - Number of Objects"
meta_cols <- c("Row",
               "Column",
               "Compound",	
               "Concentration"
)>

`tidy raw mitopaint data (mean per well) v1`
- for mPaintDR2_N4 use the following variables
<file_name <- "SF240711_mPaintDR2_N4"
batch_name <- "N4"
rm_cols <- c("Timepoint",
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
nuc_count <- "Non-border cells Selected - Number of Objects"
meta_cols <- c("Row",
               "Column",
               "Compound",	
               "Concentration"
)>

`plate drift correct mitopaint data (mean per well) v1`
- for mPaintDR2_N2 use the following variables
<file_name <- "SF240627_mPaintDR2_N2"
p_sig <- 0.05
dmso_wells <- c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")>
                
`plate drift correct mitopaint data (mean per well) v1`
- for mPaintDR2_N3 use the following variables
<file_name <- "SF240704_mPaintDR2_N3"
p_sig <- 0.05
dmso_wells <- c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")>

`plate drift correct mitopaint data (mean per well) v1`
- for mPaintDR2_N4 use the following variables
<file_name <- "SF240711_mPaintDR2_N4"
p_sig <- 0.05
dmso_wells <- c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")>

`plate drift correct mitopaint visualisations (N=3) (mean per well) v1`            
- for mPaintDR2_N2, mPaintDR2_N3, mPaintDR2_N4 use the following variables
<file_name_N1 <- "SF240627_mPaintDR2_N2"
file_name_N2 <- "SF240704_mPaintDR2_N3"
file_name_N3 <- "SF240711_mPaintDR2_N4"
file_name <- "mPaint_DR2_N1.2.3"
dmso_wells_N1 <- c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")
dmso_wells_N2 <- c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")
dmso_wells_N3 <- c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                   "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                   "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                   "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                   "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")
pos_control <- "CCCP_30"
pastel_cols <- lighten(c("#440154FF", "#238A8DFF", "#FDE725FF"), amount = 0.3)>

`robust zscore norm mitopaint data v1`
- for mPaintDR2_N2 use the following variables
<file_name <- "SF240627_mPaintDR2_N2"
dmso_wells <- c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")>

`robust zscore norm mitopaint data v1`
- for mPaintDR2_N3 use the following variables
<file_name <- "SF240704_mPaintDR2_N3"
dmso_wells <- c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")>

`robust zscore norm mitopaint data v1`
- for mPaintDR2_N4 use the following variables
<file_name <- "SF240711_mPaintDR2_N4"
dmso_wells <- c("2_9","2_2","3_2","3_9","4_9","4_2","5_2","5_9","6_9","6_2",
                "7_2","7_9","8_9","8_2","9_2","9_9","10_9","10_2","11_2","11_9",
                "12_9","12_2","13_2","13_9","14_9","14_2","15_2","15_9",
                "16_9","16_2","17_2","17_9","18_9","18_2","19_2","19_9",
                "20_9","20_2","21_2","21_9","22_9","22_2","23_2","23_9")>
                
`batch integration mitopaint data (N=3) (mean per well) v1`
- for mPaintDR2_N2, mPaintDR2_N3, and mPaintDR2_N4 use the following variables
<file_name_N1 <- "SF240627_mPaintDR2_N2"
file_name_N2 <- "SF240704_mPaintDR2_N3"
file_name_N3 <- "SF240711_mPaintDR2_N4"
file_name <- "mPaintDR2_N2_N3_N4"
k_weight <- 50>

`remove redundant mitopaint data (N=3) (mean per well) v1`
- for mPaintDR2_N2, mPaintDR2_N3, and mPaintDR2_N4 use the following variables
<file_name <- "mPaintDR2_N2_N3_N4"
integrate_state <- "integrated"
dataset_name <- "mPaintDR2_N2_N3_N4"
cor_thresh <- 0.95
var_tol <- 1e-12
excl_feats <- c("Nucleus", "Nuclei", "mTagBFP2")>

`dimensionality reduction mitopaint data (N=3) (mean per well) v1`
- for mPaintDR2_N2, mPaintDR2_N3, and mPaintDR2_N4 use the following variables
<file_name <- "mPaintDR2_N2_N3_N4"
redu_state <- "redu"
integrate_state <- "integrated"
dims_use <- 1:50
k_param <- 15
res <- 1
perplexity <- 20
max_iter <- 4000>

`PCA: dim red mitopaint vis (N=3) (mean per well) mPaintDR2_N2_N3_N4 v1`
- use the following file variables
<file_name <- "mPaintDR2_N2_N3_N4"
redu_state <- "redu"
integrate_state <- "integrated">
- use the following plot variables 
<x_lab <- paste0(
  "PC_1 (",
  round(var[1, "Percent_Variance"], 2),
  "%)"
)
y_lab <- paste0(
  "PC_2 (",
  round(var[2, "Percent_Variance"], 2),
  "%)"
)
size_title <- 10
size_axis <- 8
size_point <- 2
plot_width <- 4
plot_height <- 3.2
pastel_cols <- lighten(c("#440154FF", "#238A8DFF", "#FDE725FF"), amount = 0.3)>