################################################################################
# Pipeline to analyse features from cell painting
# Written By: Kristin Eccles
# Date: January 25th, 2024 
################################################################################

# load libraries
library(readxl)
library(dplyr)
library(stats)
require(caret)

# load data
# plate layout
df_pm <- read_excel("cell_painting_p1_p2.xlsx", sheet = "plate_stacked")
# cell painting features
df_f <- read.csv("Cellpainting_features_p2.csv")

################################################################################
# Join
#join together by well number 
df <- left_join(df_pm, df_f, by = c("Well" = "WellName"), keep = FALSE)
#evaluate character column 
# remove wells that had issues from visual inspection D7, E7
df <- subset(df, ! Well == "D7" | ! Well== "E7")
################################################################################
#### Feature Selection ####
#1. Eliminate Features without variation
df_clean <- df[,-nearZeroVar(df)]
# to see near zero variance columns 
#NZV<- nearZeroVar(df, saveMetrics = TRUE)

# 2. Eliminate features with low reproducibility among biological replicates
#test_chem <- subset(df_clean,  ! Chemical == "DMSO")


# calculate DMSO(solvent) median and median absolute deviate (MAD)

DMSO <- subset(df_clean,  Chemical == "DMSO")

# Exclude wells with low cell count
DMSO_count <- median(DMSO$Nuclei...Number.of.Objects)

test_chem$norm_cell_count <- (test_chem$Nuclei...Number.of.Objects/DMSO_count)*100
#test_chem <- subset(test_chem, norm_cell_count < 50)

#calculate column medians for DMSO
DMSO_median <- apply(DMSO[,c(6:ncol(DMSO))],2,median(x, na.rm = TRUE))
#calculate column MADs
DMSO_MAD <- apply(DMSO[,c(6:ncol(DMSO))],2,mad)













