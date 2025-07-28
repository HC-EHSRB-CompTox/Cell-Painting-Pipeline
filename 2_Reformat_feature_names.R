library(viridis)
library(RColorBrewer)
library(ggplot2)
library(cowplot)

#Categorize features
#Formats feature names from Columbus to region_channel_module_category number 
test_chem_well <- as.data.frame(list_results["test_chem_well"])

test_chem_well$plate_no <- str_extract(test_chem_well$test_chem_well.Well, "(?<=_)(\\d+)")

test_chem_well <- test_chem_well %>%
  relocate(plate_no, .after = "test_chem_well.Well")

#Plot cell count for each chemical
cell_count <- as.data.frame(list_results["test_chem_cc"])

colnames(cell_count) <- gsub("test_chem_cc.", "", colnames(cell_count))

cell_count$plate_no <- str_extract(cell_count$Well, "(?<=_)(\\d+)")

cell_count[cell_count$Chemical == ctrl_group,]$Concentration <- 1

cell_count_avg <- cell_count %>%
  group_by(Chemical, Concentration) %>%
  summarise(Average_cc = mean(cell_count),
            SD = sd(cell_count))

cell_count_avg$relative_cc <- cell_count_avg$Average_cc/cell_count_avg[cell_count_avg$Chemical == ctrl_group, ]$Average_cc

cell_count_avg$SD_rel_cc <- cell_count_avg$relative_cc*sqrt((cell_count_avg$SD/cell_count_avg$Average_cc)^2 + (cell_count_avg[cell_count_avg$Chemical == ctrl_group, ]$SD/cell_count_avg[cell_count_avg$Chemical == ctrl_group, ]$Average_cc)^2)

#cut-off for reduction in cell count
cytotoxic_conc <- cell_count_avg %>%
  filter(relative_cc < 0.2)

relative_cell_count <- ggplot(cell_count_avg[cell_count_avg$Chemical != ctrl_group,], aes(x = Concentration, y = relative_cc, colour = Chemical)) +
    geom_point() +
    geom_errorbar(aes(ymin = relative_cc - SD_rel_cc, ymax = relative_cc + SD_rel_cc), width = 0.05) +
    theme_bw() + 
    theme(legend.position="none") +
    scale_x_log10() +
    ylim(0, 1.6) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    facet_wrap(vars(Chemical), scales = "free")

plate_no <- unique(cell_count$plate_no)

plate_plot <- lapply(plate_no, function(x){
  
  plate_data <- cell_count[cell_count$plate_no == x, ]
  
  # Create an ordered factor for Chemical with ctrl_group first
  plate_data$Chemical <- factor(plate_data$Chemical, 
                                levels = c(ctrl_group, setdiff(unique(plate_data$Chemical), ctrl_group)), 
                                ordered = TRUE)
  
  y_min <- min(plate_data[plate_data$Chemical == ctrl_group,]$cell_count)
  y_max <- max(plate_data[plate_data$Chemical == ctrl_group,]$cell_count)
  
  plate_plot <- ggplot(plate_data, aes(x=Concentration, y=cell_count, colour=Chemical)) + 
    geom_point() + 
   # geom_errorbar(aes(ymin = relative_cc - SD_rel_cc, ymax = relative_cc + SD_rel_cc), width = 0.05) +
    geom_ribbon(aes(ymin = y_min, ymax = y_max), fill = "grey70", alpha = 0.5) +        
    geom_hline(yintercept = median(plate_data[plate_data$Chemical == ctrl_group,]$cell_count), linetype = "dashed", color = "red") +
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
       width = 20, height = 60, units = "cm"
  ) 

ggsave("Relative cell count.jpeg",
       relative_cell_count,
       width = 50, height = 30, units = "cm"
) 

ctrl_plot <- ggplot(cell_count[cell_count$Chemical == ctrl_group,], aes(x= as.character(plate_no), y=cell_count, colour=plate_no)) + 
  geom_point() + 
  scale_y_continuous() +
  theme(legend.position="none") +
  xlab("Plate Number") +
  ylab("Cell count") +
  geom_text(aes(label = Well), size = 2, vjust = -1, hjust = 0.5) +
  ggtitle(paste0("Vehicle Control (",ctrl_group, ")"))

ggsave("Cell count_ctrl_all plates.jpeg",
       ctrl_plot,
       width = 25, height = 25, units = "cm"
  ) 


#Exclude wells with >50% reduction in relative cell count
test_chem_well <- test_chem_well %>%
  filter(!paste(test_chem_well.Chemical, test_chem_well.Concentration) %in% paste(cytotoxic_conc$Chemical, cytotoxic_conc$Concentration))

###########################################################################

variances <- apply(test_chem_well, 2, var)
variances[1:4] <- 1

test_chem_well <- test_chem_well[, variances != 0 & !is.na(variances)]

features <- as.data.frame(colnames(test_chem_well[-c(1:4)]))
colnames(features) <- c("features")

module <- c("Axial","Compactness","Radial","Symmetry", "Gabor", "Haralick", "SER", "Intensity","intensity","Profile","Ratio", "Length", "Width", "Roundness", "Area", "position", "Position")
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

test_chem_well$Chemical <- paste0(test_chem_well$Chemical, "_",test_chem_well$Concentration, "_", test_chem_well$Well)
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

