#Principal Component Analysis of well-level data to reduce dimensionality
#Global Mahalanobis distance calculations

library(stats)
library(factoextra)
library(tcplfit2)
library(patchwork)
library(cowplot)
library(MASS)
library(gridExtra)

#test_chem_well.RData as input

#Eliminate feats with no variance across samples
variances <- apply(well_data, 2, var)
no_var <- which(variances == 0)

if(length(no_var) != 0 & is.numeric(no_var)){
  well_data <- well_data[, -no_var]
}

#Principal component analysis (PCA)
pca <- prcomp(well_data, center = TRUE, scale = TRUE)
#pca_x <- pca$x

#Scree plot
scree_plot <- fviz_eig(pca, addlabels = T, ylim = c(0,75), main = paste0("Scree plot"))  
scree_plot

#Calculate the number of PCs that describe >95% variance
cumulative_prop <- cumsum(pca$sdev^2)/sum(pca$sdev^2)

PC_90 <- length(which(cumulative_prop<0.90))+1
PC_95 <- length(which(cumulative_prop<0.95))+1
PC_99 <- length(which(cumulative_prop<0.99))+1

#Rotation Matrix
a <- ncol(pca$rotation)

pca_x <- as.data.frame(pca$x)
pca_x <- as.data.frame(cbind(chem = rownames(well_data), pca_x)) 

ctrl_pc <- pca_x %>%
  filter(grepl(ctrl_group, chem))

ctrl_pc <- ctrl_pc[,-c(1)]

#DMSO_mean <- colMeans(DMSO_pc)
ctrl_mean <- apply(ctrl_pc, 2, mean)

dat <- as.matrix(well_data) %*% pca$rotation[,1:PC_95]

#Covariance Matrix
Cov <- cov(dat)

#Checkpoint
det(Cov)
isSymmetric(Cov)

solved_cov <- ginv(Cov)

#Mahalanobis distance determination

mahal_dist <- mahalanobis(dat, ctrl_mean, solved_cov, inverted = T) 

mahal_dist <- as.data.frame(mahal_dist)

mahal_dist <- cbind(chem = rownames(well_data), mahal_dist)

sample <- as.data.frame(mahal_dist$chem)
colnames(sample) <- "sample"

sample <- tibble(sample) %>%
  separate(sample, c("chem", "concentration", "Well", "Plate"), sep = "_", remove = T)

mahal_dist <- cbind(sample, mahal_dist = mahal_dist$mahal_dist)

#mahal_dist <- mahal_dist %>%
#  dplyr::select(-c("Plate"))

mahal_dist$concentration <- as.numeric(mahal_dist$concentration)

mahal_dist <- mahal_dist %>%
  group_by(chem) %>%
  mutate(conc_count = length(unique(concentration)))

mahal_ctrl <- mahal_dist %>%
  filter(chem == ctrl_group)

mahal_dist <- mahal_dist %>%
  filter(conc_count > 3) 

mahal_dist <- rbind(mahal_dist, mahal_ctrl)

mahal_dist$exp <- exp_name

####################
##Plot Mahalanobis distances

chem_list <- unique(mahal_dist$chem)
plate <- unique(mahal_dist$Plate)

mahal_dist[mahal_dist$chem == ctrl_group,]$concentration <- 1

lapply(plate, function(plate){
  
  ggsave(paste0(results_dir, paste0(exp_name), "_", paste0("Plate",plate), "_Global Mahalanobis.jpeg"),
       
     ggplot(mahal_dist[mahal_dist$Plate == plate,], aes(x=concentration, y=mahal_dist, colour=chem)) + 
      geom_point() + 
      geom_hline(yintercept = median(mahal_ctrl$mahal_dist), linetype = "dashed", color = "maroon4") +
      geom_ribbon(aes(ymin = median(mahal_ctrl$mahal_dist) - mad(mahal_ctrl$mahal_dist), ymax = median(mahal_ctrl$mahal_dist) + mad(mahal_ctrl$mahal_dist)), fill = "thistle3", colour = NA, alpha = 0.5) +
      scale_y_continuous() +
      scale_x_log10() +
      #scale_color_manual(values = viridis(12)) +
      xlab("Concentration (\u03BCM)") +
      ylab("Mahalanobis Distance") +
      theme_classic() +
      theme(legend.position="none",
            axis.title.x = element_text(size = 17, margin = margin(t=15)),
            axis.title.y = element_text(size = 17, margin = margin(r=15)),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12)) +
      facet_wrap(vars(chem), scales = "free") +
      theme( strip.text.x = element_text(size = 12, face = "bold")),
      
      width = 50, height = 35, units = "cm") 
  
  })

ggsave(paste0(results_dir, paste0(exp_name), "_Global Mahalanobis.jpeg"),
       
       ggplot(mahal_dist, aes(x=concentration, y=mahal_dist, colour=chem)) + 
         geom_point() + 
         geom_hline(yintercept = median(mahal_ctrl$mahal_dist), linetype = "dashed", color = "maroon4") +
         geom_ribbon(aes(ymin = median(mahal_ctrl$mahal_dist) - mad(mahal_ctrl$mahal_dist), ymax = median(mahal_ctrl$mahal_dist) + mad(mahal_ctrl$mahal_dist)), fill = "thistle3", colour = NA, alpha = 0.5) +
         scale_y_continuous() +
         scale_x_log10() +
         #scale_color_manual(values = viridis(12)) +
         xlab("Concentration (\u03BCM)") +
         ylab("Mahalanobis Distance") +
         theme_classic() +
         theme(legend.position="none",
               axis.title.x = element_text(size = 17, margin = margin(t=15)),
               axis.title.y = element_text(size = 17, margin = margin(r=15)),
               axis.text.x = element_text(size = 12),
               axis.text.y = element_text(size = 12)) +
         facet_wrap(vars(chem), scales = "free") +
         theme( strip.text.x = element_text(size = 12, face = "bold")),
       
       width = 50, height = 35, units = "cm") 

#####Tcplfit2 to derive benchmark concentrations#####

test_chem_res <- mahal_dist %>%
  filter(chem !=ctrl_group)

chem_list <- unique(test_chem_res$chem)

DMSO_dist <- mahal_dist %>%
  filter(chem == ctrl_group)

vehicle_ctrl <- mahal_dist %>%
  filter(chem == ctrl_group) %>%
  summarise(Median = median(mahal_dist, na.rm=T), 
            nMad = mad(mahal_dist, constant=1.4826, na.rm=T))

#Function to model concentration-response using concRespCore of tcplfit2
conc_res_modeling <- function(test_chem, vehicle_ctrl){ 
  
  row <- list(conc = as.numeric(test_chem$concentration),
              resp = test_chem$mahal_dist,
              bmed = vehicle_ctrl$Median,
              cutoff = vehicle_ctrl$nMad,
              onesd = vehicle_ctrl$nMad,
              name = paste0(test_chem$chem[1]),
              assay = "Mahalanobis distance (cutoff = 1)")
  
  concRespCore(row, conthits = TRUE, aicc = FALSE, force.fit = FALSE,  bidirectional = TRUE,
               fitmodels=c("cnst", "hill", "poly1", "poly2", "pow", "exp2", "exp3","exp4", "exp5"))
  
}

#Model concentration-response for each chemical and plot individually  
tcpl_results <- lapply(chem_list, function(x){
  chem_data <- test_chem_res %>%
    filter(chem == x)
  tcpl_chem <- conc_res_modeling(chem_data, vehicle_ctrl)
  
  if(tcpl_chem$hitcall > 0.9){
    jpeg(file = paste0(results_dir, "_", x, "_tcplfit.jpg"))
    concRespPlot(tcpl_chem, ymin=min(chem_data$mahal_dist)-5, ymax=max(chem_data$mahal_dist)+10, draw.error.arrows = FALSE)
    dev.off()
  }
  
  tcpl_chem
}) %>%
      do.call(rbind, .)

ggsave(paste0(results_dir, "_Global_BMC.jpeg"), 
       
       ggplot(tcpl_results, aes(x=bmd, y=name, colour=name)) +
         geom_point(size=2) +
         geom_errorbar(aes(xmin = bmdl, xmax = bmdu), width = 0.2) +
         scale_x_log10() +
         theme(legend.position="none") +
         xlab("Median Best BMC (1SD above vehicle control) \u03BCM") +
         ylab("Chemicals"),
       
       width = 40, height = 20, units = "cm"
)

tcpl_results$exp <- exp_name

#Save tcpl results
write_csv(tcpl_results, file = paste0(results_dir, "_Global fitting Mahalanobis - tcplResult.csv"))

########################################################################################
mahal_dist <- mahal_dist %>%
  separate(Well, into = c("Row", "Column"), sep = "(?<=\\D)(?=\\d)")

mahal_dist$Column <- factor(mahal_dist$Column, levels = as.character(1:12))

mahal_dist$Row <- factor(mahal_dist$Row, levels = LETTERS[8:1])

mahaldist_hm <- ggplot(mahal_dist, aes(x = Column, y = Row, fill = mahal_dist)) +
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
  labs(title = "Mahalanobis distance of each well from the plate mean", x = "Columns", y = "Rows", fill = "Mahalanobis\nDistance") +
  facet_grid(rows = vars(Plate))

mahaldist_hm

#plots <- grid.arrange(cellcount_hm, mahaldist_hm, ncol = 2)

#ggsave(plots, filename = "Cell count and Mahalnobis distance heatmaps.jpeg",
#       height = 10, width =20)
