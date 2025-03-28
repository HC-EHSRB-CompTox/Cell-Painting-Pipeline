#Category level Principal component analysis, Mahalanobis distance calculations, and Concentration-response analysis

library(stats)
library(factoextra)
library(tcplfit2)
library(ggrepel)
library(cowplot)

#Filter out categories with fewer than 2 features
cat_lev_data <- cat_data %>%
  group_by(features) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  filter(count > 1)

cat_lev_data <- cat_lev_data[,-c(ncol(cat_lev_data))]

#Create a list of categories with more than 2 features (ft_cats)
samples <- colnames(cat_lev_data)[2:ncol(cat_lev_data)]
ft_cats <- unique(cat_lev_data$features)

#Eliminate categories with no variance across samples

cat_lev_data <- lapply(ft_cats, function(x){
  cat_dat <- well_data[, grep(x, colnames(well_data))]
  variances <- apply(cat_dat, 2, var)
  cat_dat <- cat_dat[, variances != 0] 
  cat_dat
  }) %>%
  do.call(cbind,.)

#PCA by category
pca_cat <- lapply(ft_cats, function(x){
  cat_dat <- cat_lev_data[, grep(x, colnames(cat_lev_data))]
  pca <- prcomp(cat_dat, center = TRUE, scale = TRUE)
  return (pca)
}) %>%
  setNames(ft_cats)

if(FALSE){
n <- 1

while(n<length(ft_cats)){
  if(!is.null(pca_cat[[n]])){
    scree_plot <- fviz_eig(pca_cat[[n]], addlabels = T, ylim = c(0,100), main = paste0(names(pca_cat[n])))
  }
  #ggsave(scree_plot, filename = paste0(names(pca_cat[n]), "_screeplot.jpg"))
  n <- n+1
}
}

#Category Mahalanobis distances
  
cat_mahal_dist <- lapply(ft_cats, function(x){

  pca <- pca_cat[[x]]
  
  cumulative_prop <- cumsum(pca$sdev^2)/sum(pca$sdev^2)
  PC_95 <- length(which(cumulative_prop<0.95))+1
  
  pca_x <- as.data.frame(pca$x)
  pca_x <- as.data.frame(cbind(chem = samples, pca_x)) 
  
  DMSO_pc <- pca_x %>% filter(grepl("DMSO", chem))
  DMSO_pc <- DMSO_pc[,-c(1)]
  DMSO_mean <- colMeans(DMSO_pc)
  
  cat_dat <- well_data[, grep(x, colnames(well_data))]
  zero_col <- sapply(cat_dat, function(col) all(col == 0))
  cat_dat <- cat_dat[,!zero_col]
  
  dat <- as.matrix(cat_dat) %*% pca$rotation[,1:PC_95]
  
  #Covariance Matrix
  Cov <- cov(dat)
  
  #Checkpoint
  det(Cov)
  isSymmetric(Cov)
  
  #Mahalanobis distance determination
  
  mahal_dist <- mahalanobis(dat, DMSO_mean, Cov, inverted = F) 
  
  mahal_dist <- as.data.frame(mahal_dist)
  
  mahal_dist <- cbind(chem = rownames(well_data), mahal_dist)
  
  sample_info <- mahal_dist$chem
  sample_info <- tibble(sample_info) %>%
    separate(sample_info, c("chem", "concentration", "well"), sep = "_", remove = T)
  
  mahal_dist <- cbind(sample_info, mahal_dist = mahal_dist$mahal_dist)
  return(mahal_dist)
  }) %>%
  setNames(ft_cats)


# Separate categorical Mahalanobis distances by chemical

chemicals <- unique(cat_mahal_dist[[1]]["chem"])
chemicals <- chemicals$chem

mahal_dist_chem <- lapply(chemicals, function(x){
  chemical <- x
  
  chem_dat <- lapply(ft_cats, function(y){
    chem_dat <- filter(cat_mahal_dist[[y]], chem == chemical)
    return(chem_dat)}) %>%
    setNames(ft_cats)
  
  chem_dat <- do.call(rbind, chem_dat)
  
  return(chem_dat)}) %>%
  setNames(chemicals) %>%
  do.call(rbind, .)

mahal_dist_chem$feat <- rownames(mahal_dist_chem)
mahal_dist_chem$feat <- str_extract(mahal_dist_chem$feat, "(?<=\\.).*?(?=\\.)")
rownames(mahal_dist_chem) <- NULL

#Plot Mahalanobis distances 

lapply (chemicals, function(x){

  plot <- ggplot(mahal_dist_chem[mahal_dist_chem$chem == x,], aes(x=concentration, y=mahal_dist, colour=feat)) +
    geom_point() +
    scale_y_continuous() +
    xlab("Concentration \u03BCM") +
    ylab("Mahalanobis Distance") +
    facet_wrap(vars(feat), scales = "free") + 
    theme(legend.position="none")

  ggsave(paste0(results_dir, plates[n],"_", paste0(x, collapse = "_"), "_Categorical mahalnobis.jpeg"), 
         plot,
         width = 40, height = 40, units = "cm")
  })


#####Tcplfit2 to derive benchmark concentrations#####

#Function to model concentration-response using concRespCore function in tcplfit2
conc_res_modeling <- function(test_chem, vehicle_ctrl){ 
  
  row <- list(conc = as.numeric(test_chem$concentration),
              resp = test_chem$mahal_dist,
              bmed = vehicle_ctrl$Median,
              cutoff = vehicle_ctrl$nMad,
              onesd = vehicle_ctrl$nMad,
              name = paste0(test_chem$chem[1], "_",test_chem$feat[1]),
              assay = "Mahalanobis distance (cutoff = 1)")
  
  concRespCore(row, conthits = TRUE, aicc = TRUE, force.fit = FALSE,  bidirectional = TRUE,
               fitmodels=c("cnst", "hill", "poly1", "poly2", "pow", "exp2", "exp3","exp4", "exp5"))
}


#Concentration-response modeling of each category

test_chem_res <- mahal_dist_chem %>%
  filter(chem !="DMSO")

chem_list <- unique(test_chem_res$chem)
feats <- unique(test_chem_res$feat)

vehicle_ctrl <- mahal_dist_chem %>%
  filter(chem == "DMSO") %>%
  summarise(Median = median(mahal_dist, na.rm=T), 
            nMad = mad(mahal_dist, constant=1.4826, na.rm=T))

#tcplfit2
tcpl_results_cat <- lapply(chemicals, function(x){
  
  chem_data <- test_chem_res[test_chem_res$chem == x,]
  
  tcpl_results <- lapply(feats, function(y){
    chem_data <- chem_data %>%
      filter(grepl(y, feat))
    tcpl_chem <- conc_res_modeling(chem_data, vehicle_ctrl)
    
    #Individual feature concentration-response curves
    #concRespPlot(tcpl_chem, ymin= min(chem_data$mahal_dist)-5, ymax=max(chem_data$mahal_dist)+10,  draw.error.arrows = FALSE)
  
    tcpl_chem
    }) %>%
    do.call(rbind,.)

  return(tcpl_results)
  }) %>%
  do.call(rbind,.)

#split name column into chemical name and feature name
split_col <- str_split_fixed(tcpl_results_cat$name, "_", 2)
tcpl_results_cat <- cbind("chem"=split_col[,1], "feat"=split_col[,2], tcpl_results_cat[,-1])

#BMC filtering criteria
 #Eliminate BMCs without ACC and BMCU/BMCL
a <- tcpl_results_cat[!is.na(tcpl_results_cat$bmdu) & !is.na(tcpl_results_cat$acc) & !is.na(tcpl_results_cat$bmdl),]

a <- a %>% filter(chem != "Sorbitol" & chem != "Staurosporine") %>% filter(hitcall>0.99)

#Model error parameter cutoff
b <- a %>%
  filter(er<0.99)

chem <- paste(unique(b$chem), collapse ="_")

#Plot all BMC, BMCU, and BMCL together
BMC_plot <- ggplot(b, aes(x=bmd, y=feat, colour=feat)) +
  geom_point(size=2) +
  geom_errorbar(aes(xmin = bmdl, xmax = bmdu), width = 0.2) +
  scale_x_log10() +
  xlab("Best BMC (1SD above vehicle control) \u03BCM") +
  ylab("Category") +
  theme(legend.position="none",
        axis.text.y = element_text(size = 6)) +
  facet_wrap(vars(chem))

ggsave(paste0(results_dir, plates[n],"_", paste0(chem_list, collapse = "_"), "_Categorical BMC.jpeg"),
       BMC_plot,
       width = 40, height = 20, units = "cm")

#Save tcpl results
write_csv(tcpl_results_cat, file = paste0(results_dir, plates[n],"_",chem,"_Category Mahalanobis - tcplResult.csv"))


