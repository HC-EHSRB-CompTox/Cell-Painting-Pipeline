#Principal Component Analysis of well-level data to reduce dimensionality
#Global Mahalanobis 

library(stats)
library(factoextra)
library(tcplfit2)

#test_chem_well.RData as input

#Eliminate feats with no variance across samples
variances <- apply(test_chem_well, 2, var)
no_var <- which(variances == 0)

test_chem_well <- test_chem_well[, -no_var]

#Principal component analysis (PCA)
pca <- prcomp(test_chem_well, center = TRUE, scale = TRUE)
#pca_x <- pca$x

#Scree plot
scree_plot <- fviz_eig(pca, addlabels = T, ylim = c(0,75), main = paste0("Scree plot"))  
scree_plot

#Calculate the number of PCs that describe >95% variance
cumulative_prop <- cumsum(pca$sdev^2)/sum(pca$sdev^2)

PC_90 <- length(which(cumulative_prop<0.90))+1
PC_95 <- length(which(cumulative_prop<0.95))+1
PC_99 <- length(which(cumulative_prop<0.99))+1

#Plot PCs (Nyffeler 2021 code)
a <- length(cumulative_prop)
  
plot(x=1:a, y=cumulative_prop, col="gray50", pch=19, cex=0.5, type="p",
     ylim=c(0,1), xlab="# of components", ylab="Proportion of variance retained", main="Principal components of HTPP U-2OS data")
  #horizontal part
  segments(x0=30, y0=0.90, x1 = PC_90, col="blue", lty='dashed')
  segments(x0=30, y0=0.95, x1 = PC_95, col="blue", lty='solid', lwd=2)
  segments(x0=30, y0=0.99, x1 = PC_99, col="blue", lty='dotted')
  #vertical part
  segments(x0=PC_90, y0=0.1, y1 = 0.90, col="blue", lty='dashed')
  segments(x0=PC_95, y0=0.1, y1 = 0.95, col="blue", lty='solid', lwd=2)
  segments(x0=PC_99, y0=0.1, y1 = 0.99, col="blue", lty='dotted')
  text(x=c(PC_90, PC_95, PC_99), y=0.05, labels=c(PC_90, PC_95, PC_99), srt=90)
  text(x=0, y=c(0.9, 0.95, 0.99),  labels=paste0(c(90,95,99), "%"), cex=0.7)

#Rotation Matrix
a <- ncol(pca$rotation)

pca_x <- as.data.frame(pca$x)
pca_x <- as.data.frame(cbind(chem = rownames(test_chem_well), pca_x)) 

DMSO_pc <- pca_x %>%
  filter(grepl("DMSO", chem))

DMSO_pc <- DMSO_pc[,-c(1)]

DMSO_mean <- colMeans(DMSO_pc)

dat <- as.matrix(test_chem_well) %*% pca$rotation[,1:PC_95]

#Covariance Matrix
Cov <- cov(dat)

#Checkpoint
det(Cov)
isSymmetric(Cov)

#Mahalanobis distance determination

mahal_dist <- mahalanobis(dat, DMSO_mean, Cov, inverted = F) 

mahal_dist <- as.data.frame(mahal_dist)

mahal_dist <- cbind(chem = rownames(test_chem_well), mahal_dist)

sample <- as.data.frame(mahal_dist$chem)
colnames(sample) <- "sample"

sample <- tibble(sample) %>%
  separate(sample, c("chem", "concentration", "Well", "Plate"), sep = "_", remove = T)

mahal_dist <- cbind(sample, mahal_dist = mahal_dist$mahal_dist)

mahal_dist <- mahal_dist %>%
  select(-c("Well", "Plate"))

mahal_dist$concentration <- as.numeric(mahal_dist$concentration)

#Plot Mahalanobis distances 
ggplot(mahal_dist, aes(x=concentration, y=mahal_dist, colour=chem)) +
  geom_point() +
  scale_y_continuous() +
  facet_wrap(vars(chem), scales = "free")

#####Tcplfit2 to derive benchmark concentrations#####

test_chem_res <- mahal_dist %>%
  filter(chem !="DMSO")

chem_list <- unique(test_chem_res$chem)

vehicle_ctrl <- mahal_dist %>%
  filter(chem == "DMSO") %>%
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
            fitmodels=c("cnst", "hill", "poly1", "poly2", "pow", "exp2", "exp3","exp4", "exp5"),
            bmr_scale = 1.349)
  
  }

#Model concentration-response for each chemical and plot individually  
tcpl_results <- lapply(chem_list, function(x){
    chem_data <- test_chem_res %>%
      filter(grepl(x, chem))
    tcpl_chem <- conc_res_modeling(chem_data, vehicle_ctrl)
    concRespPlot(tcpl_chem, ymin= min(chem_data$mahal_dist)-5, ymax=max(chem_data$mahal_dist)+10,  draw.error.arrows = FALSE)
    
    tcpl_chem
    }) %>%
      do.call(rbind, .)

#Plot all BMC, BMCU, and BMCL together
BMC_plot <- ggplot(tcpl_results, aes(x=bmd, y=name, colour=name)) +
    geom_point(size=2) +
    geom_errorbar(aes(xmin = bmdl, xmax = bmdu), width = 0.2) +
    scale_x_log10() +
    xlab("Median Best BMC (1SD above vehicle control) \u03BCM") +
    ylab("Chemicals")
  
BMC_plot

#Save tcpl results
save(tcpl_results, file = "Global fitting Mahalanobis - tcplResult_cutoff_1.RData")


