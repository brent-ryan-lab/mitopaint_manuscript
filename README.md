# mitopaint manuscript

this repo is in process, and contains the Rproj. and script files necessary to perform the data analysis and generate the visualisations in the mitopaint manuscript. 
note that the raw and processed data is not set to update in git commits until final publication, but scripts are reliant on these files for complete functionality.

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
tidy raw mitopaint data (mean per well) v1
- for mPaintDR2_N2 use the following variables
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

tidy raw mitopaint data (mean per well) v1
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

tidy raw mitopaint data (mean per well) v1
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