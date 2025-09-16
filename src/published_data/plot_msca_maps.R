library(sdmpredictors)
library(raster)
library(terra)

setwd("~/Documents/Lord_of_tHerrings/")
my_colors <- colorRampPalette(c("#5E85B8", "#EDF0C0", "#C13127"))
greys <- colorRampPalette(c("grey0", "grey50", "grey99"))

list_datasets(marine = TRUE)
a <- list_layers(marine = TRUE, datasets = "MARSPEC")
# 480   Bio-ORACLE           BO_sstmax                   Sea surface temperature (maximum)
# 481   Bio-ORACLE          BO_sstmean                      Sea surface temperature (mean)
# 482   Bio-ORACLE           BO_sstmin                   Sea surface temperature (minimum)
# 483   Bio-ORACLE         BO_sstrange                     Sea surface temperature (range)

t <- load_layers("MS_biogeo15_sst_max_5m")
t <- projectRaster(t, crs = CRS("+init=EPSG:3035"), method = "bilinear")
s <- load_layers("MS_biogeo08_sss_mean_5m")
s <- projectRaster(s, crs = CRS("+init=EPSG:3035"), method = "bilinear")
# plot(crop(t,extent(-15, 25, 65, 85)), col = my_colors(30), axes = T, box = F, colNA = NA)

cod_ext <- extent(2750000, 5500000, 2500000, 6500000)
t_cod <- crop(t, cod_ext)
plot(t_cod, col = my_colors(30), axes = F, box = F, colNA = NA)
s_cod <- crop(s, cod_ext)
sardine_ext <- extent(1200000, 6500000, 500000, 3500000)
t_sardine <- crop(t, sardine_ext)
plot(t_sardine, col = my_colors(30), axes = F, box = F, colNA = NA)
s_sardine <- crop(s, sardine_ext)

png("plots/t_cod.png", height = 10, width = 10, units = "in", res = 300)
plot(t_cod, col = my_colors(30), axes = F, box = F, colNA = NA)
dev.off()

png("plots/t_sardine.png", height = 10, width = 10, units = "in", res = 300)
plot(t_sardine, col = my_colors(30), axes = F, box = F, colNA = NA)
dev.off()

both_ext <- extent(800000, 8500000, 500000, 6000000)
t_both <- crop(t, both_ext)
hist(getValues(t_both))
t_both[t_both < 10] <- 10
t_both[t_both > 30] <- 30

plot(t_both, col = my_colors(30), axes = F, box = F, colNA = NA)

png("plots/both_ext.png", height = 8, width = 22, units = "in", res = 300)
plot(t_both, col = my_colors(30), axes = F, box = F, colNA = NA)
dev.off()


####
t <- load_layers("MS_biogeo13_sst_mean_5m")
her_points <- read.table("data/published_data_samples_coords_parsed.csv", sep = ',', header = T)
# her_points <- her_points[!grepl("MHER", her_points$id), ]
coordinates(her_points) <- ~ y + x
crs(her_points) <- crs(t)
t <- projectRaster(t, crs = CRS("+init=EPSG:3035"), method = "bilinear")
her_points <- spTransform(her_points, CRS("+init=EPSG:3035"))

plot(crop(t,extent(-2500000, 11000000, 2500000, 12000000)), col = my_colors(30), axes = T, box = F, colNA = NA)
plot(her_points, pch = 21, cex = 0.5, add = T)

t_her <- crop(t,extent(-140, 150, 30, 80))
plot(t_her, col = my_colors(30), axes = T, box = F, colNA = NA)
plot(crop(t,extent(-80, 80, 30, 90)), col = my_colors(30), axes = T, box = F, colNA = NA)
plot(her_points, pch = 21, fill = "green", add = T)

library(geobuffer)
library(tidyverse)

e <- load_layers(c(
  "MS_biogeo08_sss_mean_5m",
  "MS_biogeo09_sss_min_5m",
  "MS_biogeo10_sss_max_5m",
  "MS_biogeo11_sss_range_5m",
  "MS_biogeo12_sss_variance_5m",
  "MS_biogeo13_sst_mean_5m",
  "MS_biogeo14_sst_min_5m",
  "MS_biogeo15_sst_max_5m",
  "MS_biogeo16_sst_range_5m",
  "MS_biogeo17_sst_variance_5m"
))
samples <- her_points$id
var_df <- data.frame()
for (sample in samples){
  coord_sample <- her_points %>% filter(her_points$id == sample)
  # in a dataframe with columns x and y
  coords <- data.frame(x=as.numeric(coord_sample$y),
                       y=as.numeric(coord_sample$x))
  # get buffer around coordinates
  buff <- geobuffer_pts(xy = coords, dist_m = 50000)
  # create df for sample's row in general table
  row <- data.frame(sample=sample)
  # extract each variable's mean values in buffer
  for (n in 1:nlayers(e)){
    l <- e[[n]]
    l_means <- raster::extract(l, buff, na.rm = T, df = T, fun = mean)
    row <- cbind(row, l_means[2])
  }
  var_df <- rbind(var_df, row)
}
write.table(var_df, file = "data/var_df.csv", sep = ",", row.names = F)

####
library(corrplot)
library(tidyverse)

var_df <- read.table("data/var_df.csv", sep = ",", header = T)
vmat <- as.matrix(var_df[,c(2:10)])
rownames(vmat) <- var_df$sample
colnames(vmat) <- colnames(var_df)[c(2:10)]
cor_matrix = cor(x=vmat, y=NULL)
corrplot(cor_matrix, type = "upper", 
         method = "square", 
         addCoef.col = "black", 
         tl.col = "black")

vmat_uncorr <- as.matrix(var_df[,c(
  "MS_biogeo08_sss_mean_5m",
  "MS_biogeo11_sss_range_5m",
  "MS_biogeo13_sst_mean_5m"
)])
cor_matrix = cor(x=vmat_uncorr, y=NULL)
corrplot(cor_matrix, type = "upper", 
         method = "square", 
         addCoef.col = "black", 
         tl.col = "black")

plot(var_df$MS_biogeo13_sst_mean_5m, var_df$MS_biogeo08_sss_mean_5m, 
     cex = (var_df$MS_biogeo11_sss_range_5m - min(var_df$MS_biogeo11_sss_range_5m)) / diff(range(var_df$MS_biogeo11_sss_range_5m)) * 2 + 1,
     pch = 19, xlab = "MS_biogeo13_sst_mean_5m", ylab = "MS_biogeo08_sss_mean_5m", main = "MS_biogeo11_sss_range_5m encoded in size")

ggplot(var_df, aes(x=sample, y=MS_biogeo08_sss_mean_5m)) + 
  geom_bar(stat = "identity")





