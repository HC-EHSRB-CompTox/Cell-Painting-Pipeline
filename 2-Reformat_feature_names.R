#Categorize features
#Formats feature names from Columbus to region_channel_module_category number 

#data <- read.csv("C:/Users/EUCHO/OneDrive - HC-SC PHAC-ASPC/Documents/CompTox-Cell-Painting-Pipeline/2024-04-10_Feb28exp_Well-level results2.csv")

#Load "List_plate_1_2.RData"

plates <- names(list_all)

test_chem_well <- lapply(plates, function(x){
  well_dat <- as.data.frame(list_all[[x]]["test_chem_well"]) %>%
    mutate(plate_no = gsub("_", "",x)) %>%
    relocate(plate_no, .after = "test_chem_well.Well")
  well_dat}) %>%
  do.call(rbind,.)

rm(list_all)

#List of features

variances <- apply(test_chem_well[-c(1:4)], 2, var)
test_chem_well <- test_chem_well[, variances != 0]

features <- as.data.frame(colnames(test_chem_well[-c(1:4)]))
colnames(features) <- c("features")

module <- c("Axial","Compactness","Radial","Symmetry", "Gabor", "Haralick", "SER", "Intensity", "Profile","Ratio", "Length", "Width", "Roundness", "Area")
region <- c("Cell", "Nucleus", "Cytoplasm", "Cyto", "Membrane", "Ring")
channel <- c("AGP", "Mito", "DNA", "RNA", "ER")

features$region <- str_match(features$features, paste(region, collapse = "|"))
features$channel <- str_match(features$features, paste(channel, collapse = "|"))
features$module <- str_match(features$features, paste(module, collapse = "|"))

#Combine Gabor, Haralick, and SER in the Texture module
features['module'][features['module']=="Gabor"|features['module']=="Haralick"|features['module']=="SER"] <- "Texture"
features['module'][features['module']=="Ratio"|features['module']=="Length"|features['module']=="Width"|features['module']=="Roundness"|features['module']=="Area"] <- "Morph"

#Add category number
features <- features %>%
  group_by(region, channel, module) %>%
  mutate(cat_group = cur_group_id()) %>%
  ungroup()

features <- features %>%
  group_by(region, channel, module, cat_group) %>%
  mutate(ft_no = row_number()) %>%
  ungroup()

ft_colnames <- features %>%
  dplyr::select(region, channel, module, cat_group, ft_no)

#Concatenate region, channel, and module
ft_colnames2 <- as.data.frame(paste0(ft_colnames$region,"_",ft_colnames$channel,"_",ft_colnames$module))
colnames(ft_colnames2) <- c("features")

#Concatenate region, channel, module, and category group number
ft_colnames <- as.data.frame(paste0(ft_colnames$region,"_",ft_colnames$channel,"_",ft_colnames$module,"_",ft_colnames$cat_group, "_", ft_colnames$ft_no))

#Add chemical and concentration columns back in and combine with feature names
#sample_info <- as.data.frame(colnames(test_chem_well[c(1:2)]))

sample_info <- as.data.frame(c("Well", "Plate", "Chemical", "Concentration"))

colnames(ft_colnames) <- "features"
colnames(sample_info) <- "features"

fts <- rbind(sample_info, ft_colnames)

#Replace the original features names with the reformatted feature names
colnames(test_chem_well) <- t(fts) # run code up to here for BMDExpress formatting

#Concatenate chemical name, plate number, and well ID in one column
test_chem_well$Chemical <- paste0(test_chem_well$Chemical, "_",test_chem_well$Concentration, "_", test_chem_well$Well, "_", test_chem_well$Plate)
test_chem_well <- dplyr::select(test_chem_well, -c(Concentration, Well, Plate))

##Sets sample info as row names
rows <- test_chem_well$Chemical
test_chem_well <- test_chem_well[,-1]
rownames(test_chem_well) <- rows

#Features in rows and samples in columns
test_chem_well_cat <- t(test_chem_well)
#colnames(test_chem_well_t) <- test_chem_well_t[1,]
test_chem_well_cat <- as.data.frame(test_chem_well_cat)

test_chem_well_cat <- cbind(ft_colnames2, test_chem_well_cat)
rownames(test_chem_well_cat) <- NULL

