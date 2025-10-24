# install.packages("gradientForest", repos="http://R-Forge.R-project.org")
# install.packages("~/opt/RDA-forest-main/RDAforest_2.6.8.tar.gz")
library(RDAforest)
library(sdmpredictors)
library(tidyverse)
library(RColorBrewer)
library(viridis)
map.pal <- rev(viridis(100, option = "magma"))
# locations
locations <- c('A_Kalix_Baltic_Spring', 'B_Vaxholm_Baltic_Spring',
               'DalBoB_Atlantic_Autumn', 'DalFB_Atlantic_Spring',
               'DalGeB_Atlantic_Autumn', 'DalInB_Atlantic_Spring',
               'DalNsF_Atlantic_Autumn', 'DalNsS_Atlantic_Spring',
               'G_Gamleby_Baltic_Spring', 'HGS10_Downs_EnglishChannel_Winter',
               'HGS11_RingkobingFjord_NorthSea_Spring',
               'HGS12_BornholmBasin_Baltic_Autumn', 'HGS15_NSSH_Atlantic_Spring',
               'HGS16_Orkney_NorthSea_Autumn', 'HGS17_IsleOfMan_IrishSea_Autumn',
               'HGS18_CelticSea_Atlantic_AutumnWinter',
               'HGS19_TeelinBay_Atlantic_Winter', 'HGS1_Riga_Baltic_Spring',
               'HGS20_CapeWrath_Atlantic_Spring', 'HGS21_Hebrides_Atlantic_Mixed',
               'HGS22_CapeWrath_Atlantic_Autumn', 'HGS23_Clyde_Atlantic_Spring',
               'HGS24_Landvik_Atlantic_Spring', 'HGS25_Lindas_Atlantic_Spring',
               'HGS26_Lusterfjorden_Atlantic_Spring', 'HGS27_Gloppen_Atlantic_Spring',
               'HGS2_Riga_Baltic_Spring', 'HGS3_Riga_Baltic_Autumn',
               'HGS4_Riga_Baltic_Autumn', 'HGS5_Schlei_Baltic_Autumn',
               'HGS6_Schlei_Baltic_Spring', 'HGS71_Rugen_Baltic_Spring',
               'HGS72_Rugen_Baltic_Spring', 'HGS8_KattegatNorth_Atlantic_Spring',
               'HGS9_Greenland_Atlantic_Spring', 'H_Fehmarn_Baltic_Autumn',
               'J_Traslovslage_Baltic_Spring', 'LandvikS17_Norway_Baltic_Spring',
               'N_NorthSea_Atlantic_Autumn', 'O_Hamburgsund_Atlantic_Spring',
               'PB10_Skagerrak_Atlantic_Spring', 'PB11_Kalmar_Baltic_Spring',
               'PB12_Karlskrona_Baltic_Spring', 'PB1_HastKar_Baltic_Spring',
               'PB2_Iceland_Atlantic_Spring', 'PB4_Hudiksvall_Baltic_Spring',
               'PB5_Galve_Baltic_Spring', 'PB6_Galve_Baltic_Summer',
               'PB7_Galve_Baltic_Autumn', 'PB9_Kattegat_Atlantic_Spring',
               'PN3_CentralBaltic_Baltic_Spring', 'Q_Norway_Atlantic_Atlantic_Spring',
               'TysklandS18_Germany_Baltic')

# load genetic data - use 001 subset for now
freq <- read.table("data/selection_scans_poolseq/subsets_freqs/60.Neff.001.freq")[, -(1:2)]
colnames(freq) <- strsplit(readLines("data/published_data/60.Neff.freq", n = 1), "\t")[[1]][-c(1:2)]
freq <- freq[ , (names(freq) %in% locations)]
freq <- t(na.omit(freq))

# ecotypes:
ecotype <- data.frame(
  location = rownames(freq),
  ecotype = sapply(strsplit(rownames(freq), "_"), function(x) {
    paste(tail(x, 2), collapse = "_")})
)

# latlon
data_table <- read.table("data/selection_scans_poolseq/poolseqdata_info.csv",
                         sep = ",", header = TRUE)
data_table <- data_table[(data_table$Sample.name %in% locations), ]
latlon <- data_table %>%
  mutate(location = Sample.name, x = Longitude, y = Latitude) %>%
  dplyr::select(location, x, y) %>%
  distinct() %>% column_to_rownames(var = "location")
latlon <- latlon[rownames(freq), ]

# env
tem <- read.table("data/selection_scans_poolseq/sst_mean.txt", col.names = locations)
sal <- read.table("data/selection_scans_poolseq/sss_mean.txt", col.names = locations)
env <- pivot_longer(tem,
                       locations,
                       values_to = "sst_mean") %>% column_to_rownames("name")
env$sss_mean <- as.numeric(sal[1,])

# envc
t_layer <- load_layers("MS_biogeo13_sst_mean_5m")
s_layer <- load_layers("MS_biogeo08_sss_mean_5m")
envc <- as.data.frame(crop(t_layer,extent(-80, 30, 40, 75)), xy = TRUE)
envc <- cbind(envc, as.data.frame(crop(s_layer,extent(-80, 30, 40, 75))))
colnames(envc) <- c("x", "y", "sst_mean", "sss_mean")

#### DATA EXPLORATION
geno <- freq
# distances:
cordist=1-cor(t(geno))
# ordination:
ord=capscale(cordist~1)
pcs.d=scores(ord,scaling=1,display="sites",choices=c(1:2))
# plot proportion of variance explained
plot(ord$CA$eig/sum(ord$CA$eig),xlab="PC",ylab="proportion of variance explained")
# plot mds1-mds2
so=data.frame(scores(ord,scaling=1,display="sites"))
ggplot(so,aes(MDS1,MDS2,color=ecotype$ecotype))+
  geom_point()+coord_equal()+theme_bw()+scale_color_brewer(palette = "Paired")
### EXPLORE IBD
# converting lat, lon to great circle distances
GCD=gcd.dist(latlon)
latlon.gcd=GCD[[1]]
distGCD=GCD[[2]]
plot(as.dist(cordist)~distGCD,pch=16,cex=0.6,col=rgb(0,0,0,alpha=0.2))
protest(capscale(distGCD~1),capscale(cordist~1))

#### CLEANING PREDICTORS
cor(env)

#### EXPLORATORY RDA-FOREST ANALYSIS
### RUN GF
gf=makeGF(ord,env,pcs2keep=c(1:5))
gf$result
# rescaling to proportion of total variance (based on eigenvalues in ord1)
eigen.var=(ord$CA$eig/sum(ord$CA$eig))[names(gf$result)]
# total variance explained by model
sum(eigen.var*gf$result)
### IMPORTANCE OF PREDICTORS
# setting the number of PCs to keep
tokeep=3
# computing properly scaled importances:
imps=data.frame(importance_RDAforest(gf,ord))
# some data frame housekeeping...
names(imps)="R2"
imps$var=row.names(imps)
# reordering predictors by their importances:
imps$var=factor(imps$var,levels=imps$var[order(imps$R2)])
# plotting
ggplot(imps,aes(var,R2))+geom_bar(stat="identity")+coord_flip()+theme_bw()
### TURNOVER
gf3=makeGF(ord,env,pcs2keep=c(1:2))
plot_gf_turnovers(gf3,imps$var[1:2])
### PREDICTOR SELECTION
mm=mtrySelJack(Y=cordist,X=env,nreps=10,prop.positive.cutoff=0.5, top.pcs=tokeep,
               mintry = 1, maxtry = 2)
mm$goodvars <- c("sss_mean", "sst_mean") # cheat both variables in

#### Assessing confidence in importances; forming predictions
oj=ordinationJackknife(
  Y=cordist,
  X=env[,mm$goodvars],
  newX=envc,
  nreps=5,
  top.pcs=tokeep,
  extra=0.1
  )

for(i in 1:length(mm$goodvars)){
  plot_turnover(oj,envc[oj$goodrows,],names(oj$median.importance)[i])
}

#### Plotting maps of predicted adaptive neighborhoods
# spots on the map that are within modeled parameter range:
goods=oj$goodrows
# predictor data restricted to only those spots:
ras2=envc[which(goods),]
xy2=envc[goods,c("x","y")]
names(xy2)=c("lon","lat")
rfpreds=oj$predictions.direct
turnovers=oj$predictions.turnover
bests=names(oj$median.importance)

# plotting
plot_adaptation(rfpreds,ras2[,bests],xy2,main="turnovers",
                # options affecting PCA plot:
                rangeExp=1.5,
                scal=10,
                jitscale=0.05,
                # options affecting map and PCA colors:
                color.scheme="011",
                lighten=0.8,
                # options affecting clustering:
                cluster.guide = turnovers,
                nclust=5
)

aa <- cbind(xy2, turnovers)
aa$third <- aa$sst_mean - aa$sss_mean
bb <- cbind(xy2, rfpreds)

cl <- clara(turnovers, k=10)
cl <- cl$clustering
aa <- cbind(xy2, cl)

newrast <- rasterFromXYZ(aa)
plot(newrast, col = brewer.pal(10, "Set1"))

newrast <- rasterFromXYZ(bb)
# display.brewer.all()

newrast[[1]] <- (newrast[[1]]-newrast[[1]]@data@min) / (newrast[[1]]@data@max-newrast[[1]]@data@min)*255
newrast[[2]] <- (newrast[[2]]-newrast[[2]]@data@min) / (newrast[[2]]@data@max-newrast[[2]]@data@min)*255
newrast[[3]] <- (newrast[[3]]-newrast[[3]]@data@min) / (newrast[[3]]@data@max-newrast[[3]]@data@min)*255
plotRGB(newrast, r=1, g=2, b=3, bgalpha=0)
points(latlon$x, latlon$y)


###########
###########
###########
###########

i=40
locations[i]
# data to build a model
geno0=geno[-i,]
# our lone wolf
geno.i=geno[i,]
# distance matrix for model building
cordist0=1-cor(t(geno0))
# adding the lone wolf to the matrix, to predict its gPCs later
cordist.i=1-cor(t(rbind(geno0,geno.i)))

# Building the landscape model without the lone wolf (will take a while...)
envc_lowres <- as.data.frame(aggregate(terra::rast(envc), fact = 2, fun = mean), xy = T)
oj0=ordinationJackknife(
  Y=cordist0,
  X=env[-i,mm$goodvars],
  newX=envc_lowres,
  covariates=latlon.gcd[-i,],
  nreps=8,
  top.pcs=tokeep,
  extra=0.1
  )

# step one: getting gPCs for the other (model-building) wolves
ord0=capscale(cordist0~1)
sc0=scores(ord0, scaling=1, display="sites",choices=c(1:40))
# predicting gPCs for the same old wolves plus the new "lone wolf"
scp=predict(ord0,cordist.i,type='sp',scaling="sites")[,1:40]
# making sure the predictions for the new wolf align (mostly in terms of gPC sign) with the model
pro=procrustes(sc0,scp[-nrow(cordist.i),],scale=FALSE)
# at last, getting the predicted gPCs for the lone wolf
wi=as.vector(as.numeric(pro$rotation %*% scp[nrow(scp),]))[1:tokeep]
sc=adapt_scale(oj0$predictions.direct, quantiles = c(0.95))
# calculating environmental mismatch
home=env_mismatch(X=wi,Y=oj0,sy=envc_lowres[,1:2],sc=sc)
# plot habitat mismatch for our wolf:
plot(terra::rast(home),col=map.pal)
# this is where the wolf actually came from:
points(latlon[locations[i],1],latlon[locations[i],2],pch=0, col = "red")

# PLOT WITH MAX AND MEAN VALUES (SIMPLIFY)
mism <- terra::rast(home)
#t_both[t_both < 10] <- 10
mism[mism > min(home$env.mismatch)*1.5] <- max(home$env.mismatch)
mism[mism > min(home$env.mismatch)*1.1 & mism < min(home$env.mismatch)*1.5] <- mean(home$env.mismatch)
plot(mism,col=map.pal)
points(latlon[locations[i],1],latlon[locations[i],2],pch=0, col = "red")

# CALCULATE DISTANCE BETWEEN POINT AND MIN VAL OF MISMATCH
mism <- terra::rast(home)
which(values(mism) == min(home$env.mismatch))
minpoint <- xyFromCell(mism, which(values(mism) == min(home$env.mismatch)))
d <- sqrt((minpoint[1] - latlon[locations[i],1])^2 + (minpoint[2] - latlon[locations[i],2])^2)
d
