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
  
############################Load data and plate maps############################
print("Loading plate maps and feature files")

#Upload plate layouts for multiple plates and compile
plate_maps <- list.files(folder_path, pattern = "map", full.names = TRUE)

df_pm <- lapply(plate_maps, function (filename){
  df <- read_excel(filename, sheet = "plate_stacked")
  df$plate_no <- str_extract(filename, "(?i)(?<=plate_)\\d+")
  df
  }) %>%
  do.call(rbind, .)

if(any(is.na(df_pm$plate_no))){
  df_pm$plate_no <- 1
}

df_pm$plate_no <- as.character(df_pm$plate_no)

###Multiple vehicle controls
ctrl_group <- unique(df_pm$Control)

#Add a column for plate number and combine results for multiple plates
##Get plate number from folder name
compile_csv <- function(folder_path, plate_no) {
    rbindlist(lapply(list.files(folder_path, full.names = TRUE, pattern = "[A-H](1[0-2]|[1-9])", recursive = TRUE), fread), use.names = TRUE) %>%
    mutate(plate_no = str_extract(plate_no,"(?<=-)(\\d+)(?=\\[)|(?i)(?<=plate[ _])\\d+")) %>%  # "(?i)(?<=plate[ _])\\d+")) %>%
    relocate(plate_no, .after = "WellName")
  }

##########################
#Randomly sample 1000 cells
compile_csv <- function(folder_path, plate_no, sample_n = 1000) {
  combined_df <- rbindlist(
    lapply(
      list.files(folder_path, full.names = TRUE, pattern = "[A-H](1[0-2]|[1-9])", recursive = TRUE),
      function(file) {
        df <- fread(file)
        set.seed(42)  # For reproducibility
        if (nrow(df) > 1000) {
          df <- df[sample(.N, 1000)]  # .N is data.table's nrow
        }
        return(df)
      }
    ),
    use.names = TRUE
  ) %>%
    mutate(plate_no = str_extract(plate_no, "(?<=-)(\\d+)(?=\\[)|(?i)(?<=plate[ _])\\d+")) %>%
    relocate(plate_no, .after = "WellName")
  
  return(combined_df)
}
################################

plate_folders <- list.dirs(folder_path, recursive = FALSE, full.names = TRUE)

##Compile multiple plates##
print("Compiling feature files and modifying column names")

start_time <- Sys.time()
df_f <- rbindlist(lapply(plate_folders, function(plate_folders) compile_csv(plate_folders, basename(plate_folders))))
df_f$plate_no <- str_extract(df_f$plate_no,"\\d+")
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
ft_keep <- c("Well", "plate_no", "Nucle", "AGP", "ER", "DNA", "RNA", "Mito", "Cyto", "position", "Position")
col_keep <- grep(paste(ft_keep, collapse = "|"), names(df_f), value = TRUE)

df_f <- df_f[, col_keep]

df_f$plate_no <- as.character(df_f$plate_no)

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
df_f <- df_f[, colSums(is.na(df_f)) < nrow(df_f)]

#Join well ID and plate number
df_f$Well <- paste(df_f$Well, df_f$plate_no, sep = "_")

df_f <- df_f %>%
  dplyr::select(-c(plate_no))

################################################################################
######################Feature Selection and Normalization#######################

test_chem <- df_f

##calculate ctrl(solvent) median and median absolute deviation (MAD) using all control cells across the plates
ctrl_cells <- subset(df_f,  Chemical %in% ctrl_group)

#Average cell count (cc) per well
ctrl_avg_cc <- ctrl_cells[, c(1,3)] %>%
  group_by(Well) %>%
  add_tally(name = "cell_count") %>%
  ungroup()

ctrl_avg_cc <- ctrl_avg_cc[!duplicated(ctrl_avg_cc), ] %>% 
  dplyr::select(-c(Well)) %>% summarize_all(mean)

#Remove sample info columns (will be added back in)
ctrl_info <- as.data.frame(ctrl_cells[,c(1:3)])
ctrl_cells <- as.data.frame(ctrl_cells[,-c(1:3)])

##ctrl median
ctrl_med <- as.data.frame(t(apply(ctrl_cells, 2, function(x)median(x, na.rm=TRUE))))
ctrl_med <- ctrl_med[rep(1, nrow(test_chem)), ]

##ctrl MAD
ctrl_mad <- as.data.frame(t(apply(ctrl_cells, 2, function(x)mad(x, constant = 1, na.rm=TRUE)))*1.4826)

##Identify features with a MAD of 0
ctrl_mad_zero <- ctrl_mad[, ctrl_mad[1, ]==0 & !is.na(ctrl_mad[1, ])]

if(ncol(ctrl_mad_zero)>0){
  ctrl_mad_zero <- ctrl_mad_zero[rep(1, nrow(test_chem)), ]
  col_mad_zero <- colnames(ctrl_mad_zero)
}

##Features with non-zero MAD
ctrl_mad <- ctrl_mad[, ctrl_mad[1,]!=0 & !is.na(ctrl_mad[1, ])]
ctrl_mad <- ctrl_mad[rep(1, nrow(test_chem)), ]
col_mad <- colnames(ctrl_mad)

#Add sample info back to ctrl_cells
ctrl_cells <- cbind(ctrl_info, ctrl_cells)

#Normalize all cell-level data to ctrl SD, median, and MAD
print("Calculating standard deviation, median, and median absolute deviation (MAD) of all ctrl cells across the plates")

test_chem_cc <- test_chem[,c(1:3)] %>% 
  group_by(Well, Chemical, Concentration) %>%
  tally(name = "cell_count") %>%
  ungroup()

##Normalize features in test samples
print(paste0("Normalizing cell-level data to ctrl median and MAD"))

test_chem[, col_mad] <- (test_chem[,col_mad] - ctrl_med[,col_mad])/ctrl_mad[,col_mad]

if(ncol(ctrl_mad_zero)>0){
  test_chem[, col_mad_zero] <- (test_chem[, col_mad_zero] - ctrl_med[, col_mad_zero])
  }else{
  rm(ctrl_mad_zero)
  }

rm(ctrl_med, ctrl_mad, ctrl_mad_zero)

################################################################################
###########################Data Aggregation by Well#############################

print("Aggregating cell-level data to well-level")

#Calculate the median of all treated wells
test_chem_well <- test_chem %>%
  group_by(Well, Chemical, Concentration) %>% 
  summarize_all(median, na.rm=TRUE) %>%
  ungroup()

#z-standardization: Scale well- and treatment-level data to SD of ctrl
print("Scaling well-level data to SD of vehicle control")

ctrl_SD_w <- test_chem_well[test_chem_well$Chemical == ctrl_group, ] %>%
  dplyr::select(-c(Well, Chemical, Concentration)) %>%
  summarise_all(mean, na.rm = TRUE)

ctrl_SD_w <- ctrl_SD_w[rep(1, nrow(test_chem_well)), ] 

col_SD_w <- colnames(ctrl_SD_w[])

test_chem_well[,col_SD_w] <- test_chem_well[,col_SD_w]/ctrl_SD_w

test_chem_well <- test_chem_well %>%
  dplyr::select(-contains("Object")) %>%
  dplyr::select(-contains("obj"))

###Save normalized well-level values
print(paste0("Creating lists for combined normalized values"))

list_results <- list(
  test_chem,
  test_chem_cc,
  test_chem_well
  )

names(list_results) <- c("test_chem", "test_chem_cc", "test_chem_well")

toc()

return(list_results)

##Well- and treatment-level data to export
print("Saving well-level data as a .csv file and list_results as a .RData file")

write.csv(test_chem_well, paste0(Sys.Date(), "_",basename(folder_path), "_Well-level data.csv"), row.names = FALSE)

}

list_results <- normalizeCP(folder_path)

save(list_results, file = paste0("Normalized_CP_output_",basename(folder_path), ".RData")) 

#################################################################################################

plate_cell_count <- data.frame(list_results$test_chem_cc)

plate_cell_count <- plate_cell_count %>%
  mutate(Well = str_trim(Well)) %>%
  separate(Well, into = c("Row", "Column"), sep = "(?<=^.)") %>%
  separate(Column, into = c("Column", "Plate"), sep = "_", fill = "right", extra = "merge")

plate_cell_count$Column <- factor(plate_cell_count$Column, levels = as.character(1:12))

plate_cell_count$Row <- factor(plate_cell_count$Row, levels = LETTERS[8:1])

cellcount_hm <- ggplot(plate_cell_count, aes(x = Column, y = Row, fill = cell_count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(title = element_text(size = 13),
        axis.title.x = element_text(size = 16, margin = margin(t=15)),
        axis.title.y = element_text(size = 16, margin = margin(r=15)),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black", margin=margin(r=5)),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13, colour = "black")) +
  labs(title = "Cell count across 9 fields of view (96-well plate)", x = "Columns", y = "Rows", fill = "Cell Count") +
  facet_grid(rows=vars(Plate))


