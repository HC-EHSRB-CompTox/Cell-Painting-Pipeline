# Computational Toxicology Research Group 96-Well Format Cell Painting Analysis Pipeline

## Smaller Scale, Same Impact: Replicating High-Throughput Phenotypic Profiling in a Medium-Throughput Lab for Use in Chemical Risk Assessment
Authors: Eunnara Cho, Stephen D. Baird, Kristin M. Eccles

Cell painting visualizes morphological changes induced by toxicity using fluorescent dyes to stain cellular structures. Coupled with high-content imaging and analysis software, cell painting allows high-throughput phenotypic profiling (HTPP) to quantify phenotypic changes and estimate points of departure for toxicity assessments. Regulatory agencies have applied HTPP in 384-well plates for chemical hazard screening. In this study, established protocols for 384-well plates were adapted for use in 96-well plates to increase accessibility for lower-throughput laboratories. U-2 OS human osteosarcoma cells in 96-well plates were exposed to 12 phenotypic reference compounds for 24 hours before fixation and staining with fluorescent dyes (golgi apparatus, endoplasmic reticulum, nucleic acids, cytoskeleton, mitochondria). Four independent chemical exposures generated four biological replicates. Stained cells were imaged on an Opera Phenix, a high-content imaging system, and the Columbus analysis software extracted numerical values for 1300 morphological features. Features were normalized to control cells, followed by Principal Component Analysis and Mahalanobis distance calculation for each treatment concentration. Mahalanobis distances were modeled to calculate benchmark concentrations (BMC) for each chemical. Most BMCs differed by less than one order of magnitude across experiments, demonstrating intra-laboratory consistency. Compared to published BMCs, 10 compounds had comparable BMCs in both plate formats. Additionally, we observed a significant inverse relationship between seeding density and Mahalanobis distances, suggesting that experimental factors like cell density may influence BMCs. Overall, we demonstrate that cell painting is adaptable across formats and laboratories, supporting efforts to develop and validate it as a complementary new approach methodology to existing toxicity tests.

## Description of the pipeline
The pipeline consists of 4 modules: `1_HC_feature_normalization`, `2_Reformat_feature_names`, `3_Global PCA and Mahalnobis distances`, and `4_Category PCA and Mahalanobis distances`.

### Required packages
caret<br/>
cowplot<br/>
factoextra<br/>
ggrepel<br/>
gridExtra<br/>
MASS<br/>
patchwork<br/>
purrr<br/>
readxl<br/>
RColorBrewer<br/>
stats<br/>
stringr<br/>
tcplfit2<br/>
tidyverse<br/>
viridis

```r
install.packages(c("caret", "cowplot", "factoextra", "ggrepel", "gridExtra", "MASS", "patchwork", "purrr", "readxl", "RColorBrewer", "stats", "stringr", "tcplfit2", "tidyverse", "viridis"))

libraries <- c("caret", "cowplot", "factoextra", "ggrepel", "gridExtra", "MASS", "patchwork", "purrr", "readxl", "RColorBrewer", "stats", "stringr", "tcplfit2", "tidyverse", "viridis")

lapply(libraries, library, character.only = TRUE)

```

### Input Data Format and File Organization
The outputs from the Columbus image analysis software  are used as inputs in the `1_HC_feature_normalization` module. The main folder contains Columbus outputs for one experiment in individual sub-folders for each plate and the plate maps (one .xlsx file per plate) are saved directly in the main folder. Each sub-folder for plates contains 96 .csv files (one .csv file per well).<br/>
Each Columbus output file contains measurements for 1384 image features for all cells in the well that passed the initial filtering criteria: nuclear area (>20, <900) and cell area (>100, <6700).

### Pipeline Modules

![module 1 and 2](https://github.com/user-attachments/assets/8093f4cc-696c-40d0-a46b-2feb3cee45f1)
<br/>
<br/>
![module 3 and 4](https://github.com/user-attachments/assets/a0002a7a-37a6-4b13-beb9-fe0ee3338195)
<br/>

```r
row <- list(conc = as.numeric(test_chem$concentration),
              resp = test_chem$mahal_dist,
              bmed = vehicle_ctrl$Median,
              cutoff = vehicle_ctrl$nMad, #MAD of the median Mahalanobis distance of the solvent control wells
              onesd = vehicle_ctrl$nMad,
              name = paste0(test_chem$chem[1]),
              assay = "Mahalanobis distance (cutoff = 1)")
  
  concRespCore(row, conthits = TRUE, aicc = FALSE, force.fit = FALSE,  bidirectional = TRUE,
               fitmodels=c("cnst", "hill", "poly1", "poly2", "pow", "exp2", "exp3","exp4", "exp5"))
```
![module 3 and 4B](https://github.com/user-attachments/assets/992f180f-dfad-4c3d-bdef-c1f6b834c168)

