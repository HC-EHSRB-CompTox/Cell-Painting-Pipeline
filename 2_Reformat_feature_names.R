library(viridis)
library(RColorBrewer)
library(ggplot2)
library(cowplot)

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

#Plot cell count for each chemical

cell_count <- lapply(plates, function(x){
  cells <- as.data.frame(list_all[[x]]["test_chem_cc"])
  cells}) %>%
  do.call(rbind,.)

colnames(cell_count) <- gsub("test_chem_cc.", "", colnames(cell_count))

cell_count[cell_count$Chemical == "DMSO",]$Concentration <- 1

plate_no <- unique(cell_count$plate_no)

plate_plot <- lapply(plate_no, function(x){
  
  plate_data <- cell_count[cell_count$plate_no == x, ]
  
  # Create an ordered factor for Chemical with "DMSO" first
  plate_data$Chemical <- factor(plate_data$Chemical, 
                                levels = c("DMSO", setdiff(unique(plate_data$Chemical), "DMSO")), 
                                ordered = TRUE)
  
  y_min <- min(plate_data[plate_data$Chemical == "DMSO",]$cell_count)
  y_max <- max(plate_data[plate_data$Chemical == "DMSO",]$cell_count)
  
  plate_plot <- ggplot(plate_data, aes(x=Concentration, y=cell_count, colour=Chemical)) + 
    geom_point() + 
    geom_ribbon(aes(ymin = y_min, ymax = y_max), fill = "grey70", alpha = 0.5) +        
    geom_hline(yintercept = median(plate_data[plate_data$Chemical == "DMSO",]$cell_count), linetype = "dashed", color = "red") +
    scale_y_continuous() +
    ylim(min(plate_data$cell_count)-50, max(plate_data$cell_count)+50) +
    scale_x_log10() +
    theme(legend.position="none") +
    scale_color_brewer(palette = "Paired") +
    geom_text(aes(label = Well), size = 2, vjust = -1, hjust = 0.5) +
    facet_wrap(vars(Chemical), scales = "free") +
    ggtitle(paste0("Plate ", x))
  
    plate_plot
  }) 

nrow <- length(plate_no)

combined_plot <- do.call(plot_grid, c(plate_plot, nrow = nrow))

ggsave("Cell count by chemical_all plates.jpeg",
       combined_plot,
       width = 20, height = 40, units = "cm"
  ) 

DMSO_plot <- ggplot(cell_count[cell_count$Chemical == "DMSO",], aes(x= as.character(plate_no), y=cell_count, colour=plate_no)) + 
  geom_point() + 
  scale_y_continuous() +
  theme(legend.position="none") +
  xlab("Plate Number") +
  ylab("Cell count") +
  geom_text(aes(label = Well), size = 2, vjust = -1, hjust = 0.5) +
  ggtitle("Vehicle Control (DMSO)")

ggsave("Cell count_DMSO_all plates.jpeg",
       DMSO_plot,
       width = 25, height = 25, units = "cm"
  ) 

#List of features

###For processing Nyffeler's data
if(FALSE){
# test_chem_well <- df %>%
#   select(-c("N_objects", "N_fields"))

  test_chem_well <- test_chem_well %>%
    group_by(PlateID) %>%
    mutate(Plate = cur_group_id()) %>%
    ungroup() %>%
    select(Plate, everything(), -c("PlateID"))
  
  test_chem_well$Plate <- paste0("plate",test_chem_well$Plate)
  
  plates <- unique(test_chem_well$Plate)
  
#  test_chem_well <- test_chem_well %>%
#    group_by(Plate, Well) %>%
#    mutate(Well_2 = row_number()) %>%
#    ungroup()
  
#  Concentration <- rep(0, nrow(test_chem_well))
  
#  Well <- paste0(test_chem_well$Well, "_", test_chem_well$Well_2)
  
#  test_chem_well <- test_chem_well %>%
#    select(-c(Well_2))
  
#  test_chem_well <- cbind(Concentration = Concentration,Well, test_chem_well)
  
  test_chem_well$Chemical <- gsub("Dimethyl sulfoxide", "DMSO", test_chem_well$Chemical)
  test_chem_well$Concentration[is.na(test_chem_well$Concentration)] <- 0
  
}
####

variances <- apply(test_chem_well, 2, var)
variances[1:4] <- 1

test_chem_well <- test_chem_well[, variances != 0 & !is.na(variances)]

features <- as.data.frame(colnames(test_chem_well[-c(1:4)]))
colnames(features) <- c("features")

module <- c("Axial","Compactness","Radial","Symmetry", "Gabor", "Haralick", "SER", "Intensity","intensity","Profile","Ratio", "Length", "Width", "Roundness", "Area", "position")
region <- c("Nuclei", "nuclei","Cell", "Nucleus", "Cytoplasm", "Cyto", "Membrane", "Ring")
channel <- c("AGP", "Mito", "DNA", "RNA", "ER","Shape", "position", "Position")

features$region <- str_match(features$features, paste(region, collapse = "|"))
features$channel <- str_match(features$features, paste(channel, collapse = "|"))
features$module <- str_match(features$features, paste(module, collapse = "|"))

#Combine Gabor, Haralick, and SER in the Texture module
features['module'][features['module']=="Gabor"|features['module']=="Haralick"|features['module']=="SER"] <- "Texture"
features['module'][features['module']=="Ratio"|features['module']=="Length"|features['module']=="Width"|features['module']=="Roundness"|features['module']=="Area"] <- "Morph"
features['region'][features['region']=="nuclei"] <- "Nucleus"

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
sample_info <- as.data.frame(colnames(test_chem_well[c(1:4)]))

sample_info <- as.data.frame(c("Well", "Plate", "Chemical", "Concentration"))

colnames(ft_colnames) <- "features"
colnames(sample_info) <- "features"

fts <- rbind(sample_info, ft_colnames)

#Replace the original features names with the reformatted feature names
colnames(test_chem_well) <- t(fts) # run code up to here for BMDExpress formatting

#Concatenate chemical name, plate number, and well ID in one column

test_chem_well$Chemical <- paste0(test_chem_well$Chemical, "_",test_chem_well$Concentration, "_", test_chem_well$Well, "_", test_chem_well$Plate)
test_chem_well <- test_chem_well[, !names(test_chem_well) %in% c("Plate", "Well","Concentration")]

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

