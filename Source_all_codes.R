# Files (.csv) from each plate are saved in separate sub-folders within the main folder. Plate maps (.xlsx) are saved in the main folder
# Data from Columbus have been filtered by nuclear area (>20, <900) and cell area (>100, <6700)

#Identify the folder where the Columbus output and plate maps are stored
folder_path <- "D:/Columbus_output_platemaps/2026-02-26 EC ESRAB iPSC[4504]"

#Create a new folder to store processed data files and set as the working directory
new_dir <- paste0("D:/Analysis_output/", Sys.Date(),"_", basename(folder_path), "_Cell painting results")
dir.create(new_dir)
setwd(new_dir)

#Create experiment name (for labeling files)
exp_name <- paste("iPSC_Hep_fem")

#Identify solvent control
ctrl_group <- "DMSO"

#########################################################################

#### 1. Feature Normalization (Join 4 plates and normalize together)
source("1_HC_feature_normalization.R")

#### 2. Reformatting of feature names and grouping by category and creating cell count plots
source("2_Reformat_feature_names.R")

#Remove list_results to free up memory
rm(list_results)

#### 3. Global Mahalanobis distance calculation and BMC analysis

source("3_Global PCA and Mahalnobis distances.R")

#### 4. Categorical Mahalanobis distance calculation and BMC analysis

source("4_Category PCA and Mahalanobis distances.R")

