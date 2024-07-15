################################################################################
# Pipeline to normalize cell painting data
# Written By: Eunnara Cho
# Date: April 2024 
################################################################################

# load libraries
library(readxl)
library(dplyr)
library(stats)
library(tidyverse)
library(stringr)
library(purrr)
require(caret)
library(data.table)
library(tictoc)

# Files (.csv) from each plate are saved in separate sub-folders within the main folder. Plate maps (.xlsx) are saved in the main folder
# Data from Columbus have been filtered by nuclear area (>20, <900) and cell area (>100, <6700)

#Path to the main folder
#folder_path <- "C:/Users/Admin/Downloads/June 2024"

#setwd("C:/Users/Admin/Documents")

normalizeCP <- function(folder_path){
  tic()
  
#Create a new folder to store processed data files
#dir.create(paste0(Sys.Date(),"_", basename(folder_path), "_well_treatment_level_results"))
#setwd(paste0(Sys.Date(),"_", basename(folder_path), "_well_treatment_level_results"))
  
############################Load data and plate maps############################
print("Loading plate maps and feature files")

#Upload plate layouts for multiple plates and compile
plate_maps <- list.files(folder_path, pattern = "map", full.names = TRUE)

df_pm <- lapply(plate_maps, function (filename){
  df <- read_excel(filename, sheet = "plate_stacked")
  df$plate_no <- str_extract(filename, "(?<=plate_)\\d+")
  df
  }) %>%
  do.call(rbind, .)

if(any(is.na(df_pm$plate_no))){
  df_pm$plate_no <- 1
}

#Add a column for plate number and combine results for multiple plates
##Get plate number from folder name
compile_csv <- function(folder_path, plate_no) {
    rbindlist(lapply(list.files(folder_path, full.names = TRUE, pattern = "\\.csv$"), fread)) %>%
    mutate(plate_no = str_extract(plate_no, "(?<=plate_)\\d+")) %>%
    relocate(plate_no, .after = "WellName")
  }

sub_folders <- list.dirs(folder_path, recursive = FALSE)

##Compile multiple plates
print("Compiling feature files and modifying column names")

start_time <- Sys.time()
df_f <- rbindlist(lapply(sub_folders, function(sub_folder) compile_csv(sub_folder, basename(sub_folder))))
end_time <- Sys.time()
print(end_time - start_time)

df_f <- as.data.frame(df_f)

if(any(is.na(df_f$plate_no))){
  df_f$plate_no <- 1
}

#Modify columns and column names
##Delete unnecessary text in column names
colnames(df_f) <- sub(".* - ", "", colnames(df_f))

## Replace all spaces in colnames with "_"
colnames(df_f) <-  gsub(" ", "_", colnames(df_f))

##Delete unnecessary columns
ft_keep <- c("Well", "plate_no", "Nucle", "AGP", "ER", "DNA", "RNA", "Mito", "Cyto")
col_keep <- grep(paste(ft_keep, collapse = "|"), names(df_f), value = TRUE)

df_f <- df_f[, col_keep]

print(paste0(ncol(df_f)," features"))

################################################################################
########################Join feature data and plate maps########################
#Join together by plate number and well number then delete unnecessary columns
print("Adding plate information to the feature data")
df_f <- left_join(df_pm, df_f, by = c("Well" = "WellName", "plate_no"), keep = FALSE)

#Eliminate wells with fewer than 50 valid objects
df_f <- df_f %>%
  group_by(plate_no, Well) %>%
  add_tally(name = "obj_nom") %>%
  filter(obj_nom > 50) %>%
  ungroup()

#Remove columns with all NA
df_f <- df_f[, colSums(is.na(df_f)) != nrow(df_f)]

#Number of wells available for analysis
#nlevels(factor(df$))

################################################################################
######################Feature Selection and Normalization#######################

#Normalize all cell-level data to DMSO SD, median, and MAD
print("Calculating standard deviation, median, and median absolute deviation (MAD) across DMSO cells for each plate")

plate_list <- unique(df_pm$plate_no)

#Creates of list for each plate that contains normalized cell-, well-, and treatment-level data
list_all <- lapply(plate_list, function(x){

  x <- as.numeric(x) 

  test_chem <- subset(df_f, plate_no == x)
  
  test_chem_cc <- test_chem[,c(1:4)] %>% 
    group_by(plate_no, Well, Chemical, Concentration) %>%
    tally(name = "cell_count") %>%
    ungroup()

##calculate DMSO(solvent) median and median absolute deviation (MAD) using all control cells across the plates
  DMSO <- subset(df_f,  Chemical == "DMSO" & plate_no == x)

#Average cell count (cc) per well
  DMSO_avg_cc <- DMSO[, c(1,4)] %>%
    group_by(plate_no, Well) %>%
    add_tally(name = "cell_count") %>%
    ungroup()

  DMSO_avg_cc <- DMSO_avg_cc[!duplicated(DMSO_avg_cc), ] %>% 
    dplyr::select(-c(Well)) %>% group_by(plate_no) %>% summarize_all(mean)

#Remove sample info columns (will be added back in)
  DMSO <- as.data.frame(DMSO[,-c(1:4)])

##DMSO SD
  DMSO_SD <- as.data.frame(t(apply(DMSO, 2, function(x)sd(x, na.rm=TRUE))))

###Identify features with SD > 0
  DMSO_SD <- DMSO_SD[, DMSO_SD[1,]!=0]

##DMSO median
  DMSO_med <- as.data.frame(t(apply(DMSO, 2, function(x)median(x, na.rm=TRUE))))
  DMSO_med <- DMSO_med[rep(1, nrow(test_chem)), ]

##DMSO MAD
  DMSO_mad <- as.data.frame(t(apply(DMSO, 2, function(x)mad(x, constant = 1, na.rm=TRUE)))*1.4826)

##Identify features with a MAD of 0
  DMSO_mad_zero <- DMSO_mad[, DMSO_mad[1, ]==0]

  if(ncol(DMSO_mad_zero)>0){
    DMSO_mad_zero <- DMSO_mad_zero[rep(1, nrow(test_chem)), ]
    col_mad_zero <- colnames(DMSO_mad_zero)
    }

##Features with non-zero MAD
  DMSO_mad <- DMSO_mad[, DMSO_mad[1,]!=0]
  DMSO_mad <- DMSO_mad[rep(1, nrow(test_chem)), ]
  col_mad <- colnames(DMSO_mad)

##Normalize features in test samples
print(paste0("Normalizing Plate ", x ," cell-level data to DMSO median and MAD"))

  test_chem[, col_mad] <- (test_chem[,col_mad] - DMSO_med[,col_mad])/DMSO_mad[,col_mad]

  if(ncol(DMSO_mad_zero)>0){
    test_chem[, col_mad_zero] <- (test_chem[, col_mad_zero] - DMSO_med[, col_mad_zero])
    }else{
    rm(DMSO_mad_zero)
    }

  rm(DMSO, DMSO_med, DMSO_mad, DMSO_mad_zero)

#hist(DMSO[,1], main = "Histogram", xlab = "Values", ylab = "Frequency")

#1. Eliminate Features without variation

#a <- df_f[ ,-nearZeroVar(df_f)]

# to see near zero variance columns 
#NZV<- nearZeroVar(df, saveMetrics = TRUE)

# 2. Eliminate features with low reproducibility among biological replicates (Pearson correlation < 0.5)

################################################################################
####################Data Aggregation by Well and Treatment######################

#Well-level data
print("Aggregating cell-level data to well-level")

#Calculate the median of all treated wells
test_chem_well <- test_chem %>%
  group_by(plate_no, Well, Chemical, Concentration) %>% 
  summarize_all(median, na.rm=TRUE) %>%
  ungroup() %>%
  dplyr::select(-c(plate_no))

#Calculate the mean of all chemicals and concentrations
print("Calculating treatment-level mean values") 

test_chem_treat <- test_chem_well %>%
  group_by(Chemical, Concentration) %>%
  summarize_all(mean, na.rm=TRUE) %>%
  ungroup()

#z-standardization: Scale well- and treatment-level data to SD of DMSO
DMSO_SD_w <- DMSO_SD[rep(1, nrow(test_chem_well)), ]
col_SD_w <- colnames(DMSO_SD_w)
test_chem_well[,col_SD_w] <- test_chem_well[,col_SD_w]/DMSO_SD_w

test_chem_well <- test_chem_well %>%
  dplyr::select(-contains("Object")) %>%
  dplyr::select(-contains("obj"))
#  filter(Chemical != "Sorbitol")

DMSO_SD_t <- DMSO_SD[rep(1, nrow(test_chem_treat)), ]
col_SD_t <- colnames(DMSO_SD_t)
test_chem_treat[,col_SD_t] <- test_chem_treat[,col_SD_t]/DMSO_SD_t

#Concatenate chemical name, plate number, and well ID in one column
#test_chem_well$Chemical <- paste0(test_chem_well$Chemical, "_",test_chem_well$plate_no, "_",test_chem_well$Well)
#test_chem_well <- dplyr::select(test_chem_well, -c(Concentration, plate_no, Well))

print(paste0("Creating lists for Plate ", x))

list_plate <- list(
  test_chem,
  DMSO_SD,
  test_chem_cc,
  test_chem_well,
  test_chem_treat
  )

names(list_plate) <- c("test_chem", "DMSO_SD", "test_chem_cc", "test_chem_well", "test_chem_treat")

list_plate

}) %>%
  do.call(list, .)

names(list_all) <- lapply(as.numeric(plate_list), function(x) paste0("Plate_",plate_list[x]))

toc()

return(list_all)

##Well- and treatment-level data to export
#print("Saving well-level and treatment-level data as a .RData file")

#write.csv(test_chem_well, paste0(Sys.Date(), "_",basename(folder_path), "_Well-level data_1.csv"), row.names = FALSE)

#write.csv(test_chem_treat, paste0(Sys.Date(), "_",basename(folder_path), "_Treatment-level data.csv"), row.names = FALSE)

#save(list_all, file = paste0("Normalized_CP_output_",basename(folder_path), ".RData"))

}

list_all <- normalizeCP(folder_path)