# Files (.csv) from each plate are saved in separate sub-folders within the main folder. Plate maps (.xlsx) are saved in the main folder
# Data from Columbus have been filtered by nuclear area (>20, <900) and cell area (>100, <6700)

#Create a new folder to store processed data files
setwd("C:/Users/Admin/Documents/")

folder_path <- "C:/Users/Admin/Documents/2024-11-05 EC PFAS"
exp_name <- paste(basename(folder_path))

new_dir <- paste0(Sys.Date(),"_", basename(folder_path), "_Cell painting results")
dir.create(new_dir)
setwd(new_dir)

#Feature Normalization
source("C:/Users/Admin/Documents/CompTox-Cell-Painting-Pipeline/1_HC_feature_normalization.R")

#Reformatting of feature names and grouping by category
source("C:/Users/Admin/Documents/CompTox-Cell-Painting-Pipeline/2_Reformat_feature_names.R")  

plates <- names(list_all)
plates <- str_replace(plates, "_", "")

#Remove list_all to free up memory
rm(list_all)

#Create new folder for global results
results_folder <- paste0("Global_analysis_results")
dir.create(results_folder)
results_dir <- paste0(getwd(),"/", results_folder, "/")

####Global Mahalanobis distance calculation and BMC analysis

n <- 1

while(n<=length(plates)){
  print(paste0(plates[n]))
  
  well_data <- test_chem_well[grepl(paste0("Plate",n), rownames(test_chem_well)), ]
  #for Nyf data
  #well_data <- test_chem_well[grepl(paste0("plate",n,"-"), rownames(test_chem_well)), ]
  ##
  
  rownames(well_data) <- rownames(test_chem_well)[grepl(paste0("Plate",n), rownames(test_chem_well))]
  
  source("C:/Users/Admin/Documents/CompTox-Cell-Painting-Pipeline/3_Global PCA and Mahalnobis distances.R")
  
  n <- n+1
  }

mahal_dist_all <- do.call(rbind, mahal_dist_all)
write.csv(mahal_dist_all, paste0(results_dir, paste0(exp_name), "_Global Mahalanobis distances.csv"))

tcpl_results_all <- do.call(rbind, tcpl_results_all)
write_csv(tcpl_results_all, file = paste0(results_dir, paste0(exp_name), "_Global fitting Mahalanobis - tcplResult.csv"))


####Categorical Mahalanobis distance calculation and BMC analysis

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
  
  #print(BMC_plot)
  
  n <- n+1
  }
