

# Files (.csv) from each plate are saved in separate sub-folders within the main folder. Plate maps (.xlsx) are saved in the main folder
# Data from Columbus have been filtered by nuclear area (>20, <900) and cell area (>100, <6700)
folder_path <- "C:/Users/Admin/Documents/EC_Oct24"

#Create a new folder to store processed data files
setwd("C:/Users/Admin/Documents/")
new_dir <- paste0(Sys.Date(),"_", basename(folder_path), "_Cell painting results")
dir.create(new_dir)
setwd(new_dir)


source("C:/Users/Admin/Documents/CompTox-Cell-Painting-Pipeline/1_HC_feature_normalization.R")
#load("Normalized_CP_output_Ref_chem_definitive_rep3.RData")

source("C:/Users/Admin/Documents/CompTox-Cell-Painting-Pipeline/2_Reformat_feature_names.R")  

plates <- names(list_all)
plates <- str_replace(plates, "_", "")
  
#Remove list_all to free up memory
rm(list_all)

n <-1

while(n<=length(plates)){
  print(paste0(plates[n]))
  
  #Create new folder for global results
  results_folder <- paste0(plates[n],"_Global analysis_results")
  dir.create(results_folder)
  results_dir <- paste0(getwd(),"/", results_folder, "/")
  
  rows <- grepl(paste0(plates[n]), rownames(test_chem_well))
  well_data <- test_chem_well[rows, ]
  source("C:/Users/Admin/Documents/CompTox-Cell-Painting-Pipeline/3_Global PCA and Mahalnobis distances.R")
  
  n <- n+1
  }


n <-1

while(n<=length(plates)){
  print(paste0(plates[n]))
  
  #Create new folder for categorical results
  results_folder <- paste0(plates[n],"_Categorical analysis_results")
  dir.create(results_folder)
  results_dir <- paste0(getwd(),"/", results_folder, "/")
  
  rows <- grepl(paste0(plates[n]), rownames(test_chem_well))
  well_data <- test_chem_well[rows, ]
  
  cols <- grepl(paste0(plates[n]), colnames(test_chem_well_cat))
  feats <- as.data.frame(test_chem_well_cat$features)
  colnames(feats) <- "features"
  cat_data <- cbind("features"=feats, test_chem_well_cat[,cols])
  
  source("C:/Users/Admin/Documents/CompTox-Cell-Painting-Pipeline/4_Category PCA and Mahalanobis distances.R")
  
  print(BMC_plot)
  
  n <- n+1
  }
