#### LOAD LIBRARIES ####

# install.packages("gradientForest", repos="http://R-Forge.R-project.org")
# install.packages("~/opt/RDA-forest-main/RDAforest_2.6.8.tar.gz")
library(RDAforest)
library(sdmpredictors)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(terra)
library(geosphere)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### LISTS AND FUNCTIONS ####

# magma palette for map
map.pal <- rev(viridis(100, option = "magma"))

# unwanted populations: pacific and NW atlantic
unwanted <- c("HWS1_Japan_SeaOfJapan",
              "HWS2_PechoraSea_BarentsSea",
              "HWS3_WhiteSea_WhiteSea",
              "HWS4_KandalakshaBay_WhiteSea",
              "HWS5_KandalakshaBay_WhiteSea",
              "HWS6_Balsfjord_Atlantic",
              "PB8_Pacific_Pacific_Spring",
              "HGS9_Greenland_Atlantic_Spring",
              "DalFB_Atlantic_Spring",
              "DalInB_Atlantic_Spring",
              "DalNsS_Atlantic_Spring",
              "DalBoB_Atlantic_Autumn",
              "DalGeB_Atlantic_Autumn",
              "DalNsF_Atlantic_Autumn"
              )

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### LOAD DATA ####


# geno: genetic data - use 001 subset for now
freq <- read.table("data/selection_scans_poolseq/subsets_freqs/60.Neff.001.freq")[, -(1:2)]
colnames(freq) <- strsplit(readLines("data/published_data/60.Neff.freq", n = 1), "\t")[[1]][-c(1:2)]
freq <- freq[ , (!names(freq) %in% unwanted)]
#freq <- freq[cand$snp, ] # candidate snps
#freq <- freq[cand[cand$predictor == "sss_mean",]$snp, ] # S candidate snps
#freq <- freq[cand[cand$predictor == "sst_mean",]$snp, ] # T candidate snps
freq <- t(na.omit(freq))
geno <- freq

# latlon: spatial coordinates
latlon <- read.table("data/selection_scans_poolseq/poolseqdata_info.csv",
                     sep = ",", header = TRUE) %>%
  mutate(location = Sample.name, x = Longitude, y = Latitude) %>%
  dplyr::select(location, x, y) %>%
  distinct() %>%
  column_to_rownames(var = "location") %>%
  drop_na()
latlon <- latlon[rownames(freq), ]

# env: environmental values for samples
tem <- read.table("data/selection_scans_poolseq/sst_mean.colnames.txt", header = TRUE)
sal <- read.table("data/selection_scans_poolseq/sss_mean.colnames.txt", header = TRUE)
env <- pivot_longer(tem, colnames(tem), values_to = "sst_mean") %>%
  column_to_rownames("name")
env$sss_mean <- as.numeric(sal[1, ])
env <- env[rownames(freq), ]

# envc: dataframe from raster layers of environmental variables
t_layer <- load_layers("MS_biogeo13_sst_mean_5m")
s_layer <- load_layers("MS_biogeo08_sss_mean_5m")
envc <- as.data.frame(crop(t_layer,extent(-30, 30, 40, 75)), xy = TRUE)
envc <- cbind(envc, as.data.frame(crop(s_layer,extent(-30, 30, 40, 75))))
colnames(envc) <- c("x", "y", "sst_mean", "sss_mean")
# reduce resolution - speeds things up
envc <- as.data.frame(aggregate(terra::rast(envc), fact = 2, fun = mean), xy = T)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### RUN GRADIENT FOREST ####


# genetic distances:
cordist <- 1 - cor(t(geno))

# ordination:
ord <- capscale(cordist ~ 1)
# plot(ord$CA$eig/sum(ord$CA$eig),xlab="PC",ylab="proportion of variance explained") # -> first 2 pcs should be good
pcs.d <- scores(ord, scaling = 1, display = "sites", choices = c(1:2)) 

# converting lat, lon to great circle distances
latlon.gcd <- gcd.dist(latlon)[[1]]
distGCD <- gcd.dist(latlon)[[2]]

# make GF - use more PCs than needed
gf <- makeGF(ord, env, pcs2keep = c(1:10))
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
  nreps = 10,
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
  distances <- sqrt((p2$x - p1[1])^2 + (p2$y - p1[2])^2)
  closest_point <- p2[which.min(distances), ]
  d <- round(geosphere::distHaversine(p1, closest_point)/1000, 1)
  # Midpoint for labeling
  xm <- mean(as.numeric(p1[1], closest_point[1]))
  ym <- mean(as.numeric(p1[2], closest_point[2]))
  
  # plot!
  pdf(paste0("plots/selection_scans_poolseq/env_mismatch_ST/", rownames(geno)[i], ".mismatch.pdf"))
  plot(mism, col=map.pal)
  points(rbind(p1, p2), pch = 0, col = "red")
  lines(rbind(p1, p2), col = "blue", lwd = 2)
  text(xm, ym, labels = paste0(d, " km"), pos = 3, col = "blue", cex = 1.1)
  dev.off()
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### HABITAT SUITABILITY WITH IND DATA ####

ind.data <- read.table("Z12.subset_001.freq", col.names = c("snp", "freq")) %>%
  column_to_rownames("snp")
ind.data <- ind.data[colnames(freq), ]

geno.i <- na.omit(ind.data)
geno0 <- freq[, !is.na(ind.data)]
tokeep <- 4

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
pro = procrustes(sc0, scp[-nrow(cordist.i),], scale=FALSE)
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
distances <- sqrt((p2$x - p1[1])^2 + (p2$y - p1[2])^2)
closest_point <- p2[which.min(distances), ]
d <- round(geosphere::distHaversine(p1, closest_point)/1000, 1)
# Midpoint for labeling
xm <- mean(as.numeric(p1[1], closest_point[1]))
ym <- mean(as.numeric(p1[2], closest_point[2]))

plot(mism, col=map.pal)
points(rbind(p1, p2), pch = 0, col = "red")
lines(rbind(p1, p2), col = "blue", lwd = 2)
text(xm, ym, labels = paste0(d, " km"), pos = 3, col = "blue", cex = 1.1)
