# Files (.csv) from each plate are saved in separate sub-folders within the main folder. Plate maps (.xlsx) are saved in the main folder
# Data from Columbus have been filtered by nuclear area (>20, <900) and cell area (>100, <6700)

#Create a new folder to store processed data files
setwd("C:/Users/Admin/Documents/")

folder_path <- "C:/Users/Admin/Documents/Columbus_output_platemaps/Ref_chem_raw_data/Ref_chem_definitive_rep1"

new_dir <- paste0("C:/Users/Admin/Documents/Analysis_output/", Sys.Date(),"_", basename(folder_path), "_Cell painting results")
dir.create(new_dir)
setwd(new_dir)

#Create experiment name (for labeling files)
exp_name <- paste("Definitive_rep1")

#Identify control cells (e.g. "DMSO", "untreated")
ctrl_group <- "DMSO"

#########################################################################

#### 1. Feature Normalization (Joing 4 plates and normalize together)
source("C:/Users/Admin/Documents/CompTox-Cell-Painting-Pipeline/1_HC_feature_normalization.R")

test_chem_well <- as.data.frame(list_results["test_chem_well"])
colnames(test_chem_well) <- str_remove(colnames(test_chem_well), "test_chem_well.")
colnames(test_chem_well)[1] <- "WellID_PlateNo"
write.csv(test_chem_well, paste0("HepaRG_Exp1_well level_", ctrl_group ,".csv"), row.names = F)

#### 2. Reformatting of feature names and grouping by category
source("C:/Users/Admin/Documents/CompTox-Cell-Painting-Pipeline/2_Reformat_feature_names.R")

#Remove list_results to free up memory
rm(list_results)

#### 3. Global Mahalanobis distance calculation and BMC analysis
well_data <- test_chem_well

#Create new folder for global results
results_folder <- "Global_analysis_results"
dir.create(results_folder)
results_dir <- paste0(getwd(),"/", results_folder, "/")

source("C:/Users/Admin/Documents/CompTox-Cell-Painting-Pipeline/3_Global PCA and Mahalnobis distances.R")

#### 4. Categorical Mahalanobis distance calculation and BMC analysis
#Create new folder for categorical results
results_folder <- "Categorical analysis_results"
dir.create(results_folder)
results_dir <- paste0(getwd(),"/", results_folder, "/")

source("C:/Users/Admin/Documents/CompTox-Cell-Painting-Pipeline/4_Category PCA and Mahalanobis distances.R")

