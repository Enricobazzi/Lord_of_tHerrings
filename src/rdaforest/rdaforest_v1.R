#### LOAD LIBRARIES ####

# install.packages("gradientForest", repos="http://R-Forge.R-project.org")
# install.packages("~/opt/RDA-forest-main/RDAforest_2.6.8.tar.gz")
library(RDAforest)
library(sdmpredictors)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(geosphere)
library(adegenet)
library(ggrepel)
library(terra)
library(ggrepel)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### LISTS AND FUNCTIONS ####

# magma palette for map
map.pal <- rev(viridis(100, option = "magma"))
# list of samples to include
samples <- c("AAL1", "AAL2", "AAL3", "AF29", "AF30", "AF31", "AF8", "AK1", "AK2", "AK3", "AM27", "AM29", "AM33", "AM8", "BF16", "BF18", "BF19", "BF21", "BM14", "BM15", "BM16", "BM19", "Fehmarn3", "Fehmarn44", "Fehmarn6", "Gavle100", "Gavle54", "Gavle98", "MHER001", "MHER003", "MHER006", "MHER007", "MHER008", "MHER009", "MHER010", "MHER011", "MHER012", "MHER013", "MHER014", "MHER015", "MHER016", "MHER017", "MHER018", "MHER019", "MHER020", "MHER021", "MHER022", "MHER023", "MHER024", "MHER025", "MHER026", "MHER027", "MHER028", "MHER029", "MHER030", "MHER034", "MHER035", "MHER036", "MHER037", "MHER038", "MHER039", "MHER044", "MHER045", "MHER046", "MHER052", "MHER053", "MHER054", "MHER055", "MHER056", "MHER057", "MHER058", "MHER059", "MHER061", "MHER062", "MHER063", "MHER065", "MHER066", "NorthSea13", "NorthSea19", "NorthSea34", "NSSH33", "NSSH34", "NSSH36", "Z12", "Z14", "Z4")

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### LOAD DATA ####

# load data table with coordinates and samples (manually selected for now)
full_table <- read.table("data/samples_table.csv", sep = ",", header = TRUE, na.strings = "UNKNOWN")
data_table <- full_table[full_table$sample_id %in% samples, ]

# geno: genetic data
gt_data <- read.PLINK("data/rdaforest/raws/modern_samples.step3.biall_acanaf.02missmax_0fmiss_005maf.raw")
gt_data <- as.matrix(gt_data)
gt_data <- gt_data[rownames(gt_data) %in% samples, ]
geno <- gt_data
# filter candidate SNPs
cand_snps <- read.table("data/rdaforest/cand_snps.csv", sep = ",", header = TRUE)
#geno <- gt_data[, cand_snps$snp]

#PCA <- rda(geno, scale=T)
#PCs <- data.frame(scores(PCA, choices=1:4, display="sites", scaling=0))
#PCs$region <- data_table$region
#PCs$sample <- data_table$sample_id
#PCs$spawn <- data_table$spawn
#
#ggplot()+
#  geom_point(data = PCs, aes(x = PC1, y = PC2, fill = spawn), shape = 21)

# latlon: spatial coordinates
latlon <- data.frame(
  x = as.numeric(data_table$x),
  y = as.numeric(data_table$y)
)
rownames(latlon) <- samples

# env: environmental values for samples
env <- read.table("data/rdaforest/etable.csv", sep=",")[, c(2,3)]
# eventually subset but for now only useful samples are in etable.csv
rownames(env) <- samples

# envc: dataframe from raster layers of environmental variables
t_layer <- load_layers("MS_biogeo13_sst_mean_5m")
s_layer <- load_layers("MS_biogeo08_sss_mean_5m")
envc <- as.data.frame(crop(t_layer,extent(-30, 30, 40, 75)), xy = TRUE)
envc <- cbind(envc, as.data.frame(crop(s_layer,extent(-30, 30, 40, 75))))
colnames(envc) <- c("x", "y", "sst_mean", "sss_mean")
# reduce resolution - speeds things up
envc <- as.data.frame(terra::aggregate(terra::rast(envc), fact = 2, fun = mean), xy = T)


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### RUN GRADIENT FOREST ####

# genetic distances:
cordist <- 1 - cor(t(geno))

# ordination:
ord <- capscale(cordist ~ 1)
# plot(ord$CA$eig/sum(ord$CA$eig),xlab="PC",ylab="proportion of variance explained") # -> first 4 pcs should be good
pcs.d <- scores(ord, scaling = 1, display = "sites", choices = c(1:4)) 

# converting lat, lon to great circle distances
latlon.gcd <- gcd.dist(latlon)[[1]]
distGCD <- gcd.dist(latlon)[[2]]

# make GF - use more PCs than needed
gf <- makeGF(ord, env, pcs2keep = c(1:20))
# gf$result
# rescaling to proportion of total variance (based on eigenvalues in ord1)
eigen.var <- (ord$CA$eig / sum(ord$CA$eig))[names(gf$result)]
# total variance explained by model
# sum(eigen.var*gf$result)

# setting the number of PCs to keep
tokeep <- 4
# computing properly scaled importances:
imps <- data.frame(importance_RDAforest(gf, ord))
names(imps)="R2"
imps$var=row.names(imps)
imps$var=factor(imps$var, levels = imps$var[order(imps$R2)])
# plot_gf_turnovers(gf ,imps$var)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### FORMING PREDICTIONS ####

# run model
oj <- ordinationJackknife(
  Y = cordist,
  X = env[, c("sss_mean", "sst_mean")],
  newX = envc,
  nreps = 30,
  top.pcs = tokeep,
  extra = 0.1
)
# spots on the map that are within modeled parameter range:
goods <- oj$goodrows
# predictor data restricted to only those spots:
ras2 <- envc[which(goods), ]
xy2 <- envc[goods, c("x","y")]
names(xy2) <- c("lon","lat")
rfpreds <- oj$predictions.direct
turnovers <- oj$predictions.turnover
bests <- names(oj$median.importance)

# plots
rast_rfpreds <- rasterFromXYZ(cbind(xy2, rfpreds))
rast_rfpreds[[1]] <- (rast_rfpreds[[1]]-rast_rfpreds[[1]]@data@min) / (rast_rfpreds[[1]]@data@max-rast_rfpreds[[1]]@data@min)*255
rast_rfpreds[[2]] <- (rast_rfpreds[[2]]-rast_rfpreds[[2]]@data@min) / (rast_rfpreds[[2]]@data@max-rast_rfpreds[[2]]@data@min)*255
rast_rfpreds[[3]] <- (rast_rfpreds[[3]]-rast_rfpreds[[3]]@data@min) / (rast_rfpreds[[3]]@data@max-rast_rfpreds[[3]]@data@min)*255
plotRGB(rast_rfpreds, r = 1, g = 2, b = 3, bgalpha = 0)
points(latlon$x, latlon$y)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### HABITAT SUITABILITY TEST ####

for (i in 1:length(rownames(geno))){
  print(rownames(geno)[i])
  # points(latlon[rownames(geno)[i],1],latlon[rownames(geno)[i],2],pch=0, col = "red")
  
  # data to build a model
  geno0=geno[-i,]
  # our lone wolf
  geno.i=geno[i,]
  # distance matrix for model building
  cordist0=1-cor(t(geno0))
  # adding the lone wolf to the matrix, to predict its gPCs later
  cordist.i=1-cor(t(rbind(geno0,geno.i)))
  
  # Building the landscape model without the lone wolf (will take a while...)
  oj0 <- ordinationJackknife(
    Y = cordist0,
    X = env[-i, c("sss_mean", "sst_mean")],
    newX = envc,
    nreps = 10,
    top.pcs = tokeep,
    extra = 0.1
  )
  
  # step one: getting gPCs for the other (model-building) wolves
  ord0 = capscale(cordist0 ~ 1)
  sc0 = scores(ord0, scaling = 1, display = "sites", choices=c(1:5))
  # predicting gPCs for the same old wolves plus the new "lone wolf"
  scp = predict(ord0, cordist.i, type='sp', scaling="sites")[,1:5]
  # making sure the predictions for the new wolf align (mostly in terms of gPC sign) with the model
  pro = procrustes(sc0,scp[-nrow(cordist.i),],scale=FALSE)
  # at last, getting the predicted gPCs for the lone wolf
  wi = as.vector(as.numeric(pro$rotation %*% scp[nrow(scp),]))[1:tokeep]
  sc = adapt_scale(oj0$predictions.direct)[[2]]
  # calculating environmental mismatch
  home = env_mismatch(X = wi,Y = oj0, sy = envc[,1:2], sc = sc)
  # plot habitat mismatch for our wolf:
  # plot(terra::rast(home),col=map.pal)
  # this is where the wolf actually came from:
  # points(latlon[rownames(geno)[i],1],latlon[rownames(geno)[i],2],pch=0, col = "red")
  
  # CALCULATE DISTANCE BETWEEN POINT AND MIN VAL OF MISMATCH
  mism <- terra::rast(home)
  which(values(mism) == min(home$env.mismatch))
  p1 <- data.frame(x=latlon[rownames(geno)[i],1], y=latlon[rownames(geno)[i],2])
  p2 <- data.frame(xyFromCell(mism, which(values(mism) == min(home$env.mismatch))))
  #distances <- sqrt((p2$x - p1[1])^2 + (p2$y - p1[2])^2)
  #closest_point <- p2[which.min(distances), ]
  #d <- round(geosphere::distHaversine(p1, closest_point)/1000, 1)
  ## Midpoint for labeling
  #xm <- mean(as.numeric(p1[1], closest_point[1]))
  #ym <- mean(as.numeric(p1[2], closest_point[2]))
  
  # plot!
  pdf(paste0("plots/rdaforest/env_mismatch/", rownames(geno)[i], ".mismatch.pdf"))
  plot(mism, col=map.pal)
  points(p1, pch = 0, col = "blue")
  #points(rbind(p1, p2), pch = 0, col = "red")
  #lines(rbind(p1, p2), col = "blue", lwd = 2)
  text(p1[1], p1[2], labels = rownames(geno)[i], pos = 1, col = "blue", cex = 1.1)
  dev.off()
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### HABITAT SUITABILITY WITH IND DATA ####

# read individual genetic data
sample <- "ND157"
hist.samples <- c("HER001", "HER002", "HER003", "HER004", "HER005", "HER006", "HER007", "HER008", "HER009", "HER012", "HER013", "HER015", "HER016", "HER018", "HER035", "HER036", "HER037", "HER038", "HER039", "HER041", "HER042", "HER043", "HER045", "HER046", "HER047", "HER051", "HER056", "HER059", "HER065", "HER067", "HER069", "HER102", "HER104", "HER106", "HER109", "HER110", "HER111", "HER112", "HER113", "HER114", "HER115", "HER116", "HER117", "HER118", "HER121", "HER123", "HER124", "HER127", "HER128", "HER130", "HER132", "HER133", "HER134", "HER135", "HER136", "ND001", "ND009", "ND010", "ND011", "ND021", "ND022", "ND023", "ND024", "ND026", "ND028", "ND029", "ND030", "ND031", "ND032", "ND034", "ND036", "ND038", "ND039", "ND040", "ND157")

for (sample in hist.samples){
  
  ind.data <- read.PLINK(paste0("data/rdaforest/raws/", sample, ".02missmax_0fmiss_005maf.snps.raw"))
  ind.data <- as.matrix(ind.data)
  colnames(ind.data) <- sub("_[A-Za-z0-9]+$", "", colnames(ind.data))
  # geno.i <- ind.data[, intersect(colnames(ind.data), sub("_[A-Za-z0-9]+$", "", cand_snps$snp)), drop = FALSE]
  # geno.i <- ind.data[, intersect(colnames(ind.data), sub("_[A-Za-z0-9]+$", "", candtau$SNP)), drop = FALSE]
  geno.i <- ind.data
  
  if (ncol(geno.i) < 100){
    print("no snps! go next")
    next
  }
  
  geno0 <- geno
  colnames(geno0) <- sub("_[A-Za-z0-9]+$", "", colnames(geno0))
  geno0 <- geno0[, colnames(geno.i)]
  
  # distance matrix for model building
  cordist0=1-cor(t(geno0))
  if (any(is.na(cordist0)) == TRUE){
    print("no variation in at least one sample! go next")
    next
  }
  # adding the lone wolf to the matrix, to predict its gPCs later
  cordist.i=1-cor(t(rbind(geno0,geno.i)))
  if (any(is.na(cordist.i)) == TRUE){
    print("no variation in at least one sample! go next")
    next
  }
  # Building the landscape model without the lone wolf (will take a while...)
  oj0 <- ordinationJackknife(
    Y = cordist0,
    X = env[, c("sss_mean", "sst_mean")],
    newX = envc,
    nreps = 10,
    top.pcs = tokeep,
    extra = 0.1
  )
  
  # step one: getting gPCs for the other (model-building) wolves
  ord0 = capscale(cordist0 ~ 1)
  sc0 = scores(ord0, scaling = 1, display = "sites", choices=c(1:5))
  # predicting gPCs for the same old wolves plus the new "lone wolf"
  scp = predict(ord0, cordist.i, type='sp', scaling="sites")[,1:5]
  # making sure the predictions for the new wolf align (mostly in terms of gPC sign) with the model
  pro = procrustes(sc0, scp[-nrow(cordist.i),], scale=FALSE)
  # at last, getting the predicted gPCs for the lone wolf
  wi = as.vector(as.numeric(pro$rotation %*% scp[nrow(scp),]))[1:tokeep]
  sc = adapt_scale(oj0$predictions.direct)[[2]]
  # calculating environmental mismatch
  home = env_mismatch(X = wi,Y = oj0, sy = envc[,1:2], sc = sc)
  # plot habitat mismatch for our wolf:
  pdf(paste0("plots/rdaforest/env_mismatch/", sample, ".mismatch.pdf"))
  plot(terra::rast(home),col=map.pal)
  # this is where the wolf actually came from:
  points(as.numeric(full_table[full_table$sample_id == sample, "x"]),
         as.numeric(full_table[full_table$sample_id == sample, "y"]),
         pch = 0, col = "blue")
  text(as.numeric(full_table[full_table$sample_id == sample, "x"]),
       as.numeric(full_table[full_table$sample_id == sample, "y"]),
       labels = sample, pos = 1, col = "blue", cex = 1.1)
  dev.off()
  ##
  sample_df <- as.data.frame(scp)
  sample_df$population <- c(data_table$region, "HIST")
  col <- rep(NA, length(rownames(sample_df)))
  col[grep("BALTIC", sample_df$population)] <- "#E31A1C"
  col[grep("NE-ATL", sample_df$population)] <- "#3B528BFF"
  col[grep("TRANS", sample_df$population)] <- "#A035AF"
  col[grep("HIST", sample_df$population)] <- "green3"
  sample_df$color <- col
  # plot!
  sample_labels <- rownames(sample_df)
  ax_expl_constrain = round(x = (ord0$CA$eig / sum(ord0$CA$eig)) * 100,
                            digits = 2)
  x = 1
  y = 2
  rda_plot <- ggplot() +
    geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.3) +
    geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.3) +
    geom_point(aes(x = sample_df[,x], y = sample_df[,y]),
               shape = 21, cex = 3, fill = sample_df$color, alpha = 0.8) +
    geom_text_repel(aes(x = sample_df[, x],
                        y = sample_df[, y],
                        label = sample_labels),
                    family = "Verdana",
                    size = 2,
                    box.padding = 0.3,
                    point.padding = 0.2,
                    max.overlaps = Inf) +
    theme_bw(base_size = 14, base_family = "Verdana") +
    #theme(panel.background = element_blank(), panel.grid = element_blank(),
    #      plot.background = element_blank()) +
    theme(panel.grid = element_blank()) +
    xlab(paste0(colnames(sample_df)[x], " - ", ax_expl_constrain[x], "%")) +
    ylab(paste0(colnames(sample_df)[y], " - ", ax_expl_constrain[y], "%"))
  ggsave(filename = paste0("plots/rdaforest/rda/", sample, ".rdaforest.pca.png"),
         plot = rda_plot, height =  10, width = 10)
}

# sample <- "ND032"
# ind.data <- read.PLINK(paste0("data/rdaforest/raws/", sample, ".02missmax_0fmiss_005maf.snps.raw"))
# ind.data <- as.matrix(ind.data)
# colnames(ind.data) <- sub("_[A-Za-z0-9]+$", "", colnames(ind.data))
# # geno.i <- ind.data[, intersect(colnames(ind.data), sub("_[A-Za-z0-9]+$", "", cand_snps$snp)), drop = FALSE]
# geno.i <- ind.data
# geno0 <- geno
# colnames(geno0) <- sub("_[A-Za-z0-9]+$", "", colnames(geno0))
# geno0 <- geno0[, colnames(geno.i)]
# cordist.i=1-cor(t(rbind(geno0,geno.i)))
# hc <- hclust(as.dist(cordist.i), method = "average")
# plot(hclust(as.dist(cordist.i), method = "average"))


####
pdf("plots/rdaforest/modern.hc.allsnps.v2.pdf", width = 15, height = 5)
hc <- hclust(as.dist(1-cor(t(geno))), method = "average")
plot(hc)
dev.off()

pdf("plots/rdaforest/modern.hc.candsnps.pdf", width = 15, height = 5)
hc <- hclust(as.dist(1-cor(t(geno[, cand_snps$snp]))), method = "complete")
plot(hc)
dev.off()

# geno0 <- geno[, cand_snps$snp]
# geno0 <- geno[, candtau$SNP]
geno0 <- geno
colnames(geno0) <- sub("_[A-Za-z0-9]+$", "", colnames(geno0))

hist.samples <- c("HER065", "HER109", "HER127", "HER132", "HER133", "HER135", "ND022", "ND157")
for (sample in hist.samples){
  ind.data <- read.PLINK(paste0("data/rdaforest/raws/", sample, ".02missmax_0fmiss_005maf.snps.raw"), quiet = T)
  ind.data <- as.matrix(ind.data)
  colnames(ind.data) <- sub("_[A-Za-z0-9]+$", "", colnames(ind.data))
  ind.data <- ind.data[, (colnames(ind.data) %in% colnames(geno0)), drop = FALSE]
  if (ncol(ind.data) < 500){
    print(paste(sample, " has too few snps going next"))
    next
  }
  geno0 <- geno0[, (colnames(geno0) %in% colnames(ind.data))]
  geno0 <- rbind(geno0, ind.data)
  print(paste(sample, " :", ncol(geno0)))
}

hc <- hclust(as.dist(1-cor(t(geno0))), method = "complete")
plot(hc)

# "simple" pca
PCA <- rda(geno0[rownames(geno0[!data_table$region == "BALTIC", ]), ], scale=T)
sample_df <- data.frame(scores(PCA, choices=c(1:4), display="sites", scaling="none"))
sample_df$population <- column_to_rownames(full_table, var = "sample_id")[rownames(sample_df), "region"]
# colors
col <- rep(NA, length(rownames(sample_df)))
col[grep("BALTIC", sample_df$population)] <- "#E31A1C"
col[grep("NE-ATL", sample_df$population)] <- "#3B528BFF"
col[grep("TRANS", sample_df$population)] <- "#A035AF"
sample_df$color <- col
# labels
sample_labels <- rownames(sample_df)
# variance explained by rda axes
rda_ax_expl_constrain = round(x = (PCA$CA$eig / sum(PCA$CA$eig)) * 100,
                              digits = 2)
# plot!
x = 1
y = 3
pca_plot <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.3) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.3) +
  geom_point(aes(x = sample_df[,x], y = sample_df[,y]),
             shape = 21, cex = 3, fill = sample_df$color, alpha = 0.8) +
  geom_text_repel(aes(x = sample_df[, x],
                      y = sample_df[, y],
                      label = sample_labels),
                  family = "Verdana",
                  size = 2,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  max.overlaps = Inf) +
  theme_bw(base_size = 14, base_family = "Verdana") +
  theme(panel.grid = element_blank()) +
  xlab(paste0(colnames(sample_df)[x], " - ", rda_ax_expl_constrain[x], "%")) +
  ylab(paste0(colnames(sample_df)[y], " - ", rda_ax_expl_constrain[y], "%"))
pca_plot

