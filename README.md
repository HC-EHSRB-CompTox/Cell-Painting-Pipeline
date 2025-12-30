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
```

### Input Data Format and File Organization
The outputs from the Columbus image analysis software  are used as inputs in the `1_HC_feature_normalization` module. The main folder contains Columbus outputs for one experiment in individual sub-folders for each plate and the plate maps (one .xlsx file per plate) are saved directly in the main folder. Each sub-folder for plates contains 96 .csv files (one .csv file per well).<br/>
Each Columbus output file contains measurements for 1384 image features for all cells in the well that passed the initial filtering criteria: nuclear area (>20, <900) and cell area (>100, <6700).

![module 1 and 2](https://github.com/user-attachments/assets/bdc42adb-04ee-42b0-8bb0-a6ed5dfbd641)


### Module 3: Global Principal Component Analysis (PCA) and Mahalanobis Distance Calculation
`3_Global PCA and Mahalnobis distances`

### Module 4: Categorical PCA and Mahalanobis Distance Calculation
`4_Category PCA and Mahalanobis distances`

