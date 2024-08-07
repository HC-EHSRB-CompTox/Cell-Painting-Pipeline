setwd("C:/Users/Admin/Documents/CompTox-Cell-Painting-Pipeline")

# Files (.csv) from each plate are saved in separate sub-folders within the main folder. Plate maps (.xlsx) are saved in the main folder
# Data from Columbus have been filtered by nuclear area (>20, <900) and cell area (>100, <6700)
folder_path <- "C:/Users/Admin/Downloads/June 2024"

source("1_HC_feature_normalization.R")
source("2_Reformat_feature_names.R")  
  
plates <- names(list_all)
plates <- str_replace(plates, "_", "")

n <-1

while(n<=length(plates)){
  print(paste0(plates[n]))
  
  rows <- grepl(paste0(plates[n]), rownames(test_chem_well))
  well_data <- test_chem_well[rows, ]
  source("3_Global PCA and Mahalnobis distances.R")
  
  n <- n+1
  }


n <-1

while(n<=length(plates)){
  print(paste0(plates[n]))
  
  rows <- grepl(paste0(plates[n]), rownames(test_chem_well))
  well_data <- test_chem_well[rows, ]
  
  cols <- grepl(paste0(plates[n]), colnames(test_chem_well_cat))
  feats <- as.data.frame(test_chem_well_cat$features)
  colnames(feats) <- "features"
  cat_data <- cbind("features"=feats, test_chem_well_cat[,cols])
  
  source("4_Category PCA and Mahalanobis distances.R")
  
  n <- n+1
  }
