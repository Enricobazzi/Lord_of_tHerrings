library(RDAforest)
library(tidyverse)
library(sdmpredictors)
library(terra)
library(biscale)
library(cowplot)
library(viridis)

map.pal <- rev(viridis(100, option = "magma"))

get_samples_from_file <- function(file_path) {
  samples <- read.table(file_path)[, 1] |> as.character()
  return(samples)
}

read_matrix_from_file <- function(file_path) {
  mat <- as.matrix(read.table(file_path))
  return(mat)
}

get_samples_matrix <- function(matrix_file, samples, all_samples) {
  mat <- read_matrix_from_file(matrix_file)
  if (nrow(mat) == length(samples) && ncol(mat) == length(samples)) {
    rownames(mat) <- samples
    colnames(mat) <- samples
    return(mat)
  } else {
    sample_indices <- which(all_samples %in% samples)
    filtered_mat <- mat[sample_indices, sample_indices]
    rownames(filtered_mat) <- all_samples[sample_indices]
    colnames(filtered_mat) <- all_samples[sample_indices]
    return(filtered_mat)
  }
}

get_env_from_file <- function(file_path) {
  etable <- read.table(file_path, sep = ",", header = TRUE, na.strings = c("NA"))
  rownames(etable) <- etable$sample
  # remove sample x y columns
  etable <- etable[, !(colnames(etable) %in% c("sample", "x", "y"))]
  return(etable)
}

get_latlon_from_file <- function(file_path) {
  etable <- read.table(file_path, sep = ",", header = TRUE, na.strings = c("NA"))
  latlon <- data.frame(
    lon = etable$x,
    lat = etable$y
  )
  rownames(latlon) <- etable$sample
  return(latlon)
}

build_envc <- function(
  layer_names = c("MS_biogeo13_sst_mean_5m", "MS_biogeo08_sss_mean_5m"),
  extent_coords = c(-30, 30, 40, 75),
  layer_labels = NULL) {
  # extent_coords: c(xmin, xmax, ymin, ymax)
  # Load first layer and crop to get base dataframe with x,y
  first_layer <- load_layers(layer_names[1])
  envc <- as.data.frame(crop(first_layer, extent(extent_coords)), xy = TRUE)
  # If no labels provided, use layer names
  if (is.null(layer_labels)) {
    layer_labels <- layer_names
  }
  # Load and append remaining layers if any
  if (length(layer_names) > 1) {
    for (i in 2:length(layer_names)) {
      layer <- load_layers(layer_names[i])
      layer_df <- as.data.frame(crop(layer, extent(extent_coords)))
      envc <- cbind(envc, layer_df)
    }
  }
  # Set column names: x, y, then layer labels
  colnames(envc) <- c("x", "y", layer_labels)
  envc <- as.data.frame(terra::aggregate(terra::rast(envc), fact = 2, fun = mean), xy = TRUE)
  return(envc)
}

get_dist_matrix <- function(matrix_file, samples, all_samples) {
  cov_matrix <- get_samples_matrix(matrix_file, samples, all_samples)
  dist_matrix <- as.dist(1 - cov2cor(cov_matrix))
  dist_matrix <- as.matrix(dist_matrix)
  return(dist_matrix)
}

get_sample_latlon <- function(sample) {
  latlon <- read.table(
    "data/samples_table.csv",
    sep = ",",
    header = TRUE,
    na.strings = c("NA")
  )[, c("sample_id", "x", "y")]
  rownames(latlon) <- latlon$sample_id
  latlon <- latlon[, c("x", "y")]
  colnames(latlon) <- c("lon", "lat")
  return(latlon[sample, ])
}

# dummy arguments
all_samples_name <- "wp1_all"
model_samples_name <- "wp1_modern"
sites_name <- "sf7_sites"

# load samples - full set (complete matrix) and model set (used for modelling)
all_samples_file <- paste0("data/angsd_matrix/bamlists/", all_samples_name, ".sample_list.txt")
model_samples_file <- paste0("data/angsd_matrix/bamlists/", model_samples_name, ".sample_list.txt")
all_samples <- get_samples_from_file(all_samples_file)
model_samples <- get_samples_from_file(model_samples_file)
all_samples <- all_samples[all_samples != "HER135"] # Exclude sample HER135 - should be fixed in the future
model_samples <- model_samples[model_samples != "HER135"] # Exclude sample HER135 - should be fixed in the future

# load etable and latlon from etable file
etable_file <- paste0("data/rdaforest/etable.", model_samples_name, ".csv")
env <- get_env_from_file(etable_file)
latlon <- get_latlon_from_file(etable_file)
env <- env[model_samples != "HER135", ] # Exclude sample HER135 - should be fixed in the future
latlon <- latlon[model_samples != "HER135", ] # Exclude sample HER135 - should be fixed in the future
# converting lat, lon to great circle distances
latlon.gcd <- gcd.dist(latlon)[[1]]
distGCD <- gcd.dist(latlon)[[2]]

# build envc
envc <- build_envc(layer_labels = colnames(env))

# build dist matrix from covariance matrix file
matrix_file <- paste0("data/angsd_matrix/", all_samples_name, ".", sites_name, ".pcangsd.cov")
cordist <- get_dist_matrix(matrix_file, model_samples, all_samples)

# run RDA forest
ord <- capscale(cordist ~ 1)
pcs.d <- scores(ord, scaling = 1, display = "sites", choices = c(1:10))
# make GF - use more PCs than needed
gf <- makeGF(ord, env, pcs2keep = c(1:5))
# rescaling to proportion of total variance (based on eigenvalues in ord1)
eigen.var <- (ord$CA$eig / sum(ord$CA$eig))[names(gf$result)]
# total variance explained by model
sum(eigen.var * gf$result)
# setting the number of PCs to keep
tokeep <- 5

# computing properly scaled importances:
imps <- data.frame(importance_RDAforest(gf, ord))
names(imps) <- "R2"
imps$var <- row.names(imps)
imps$var <- factor(imps$var, levels = imps$var[order(imps$R2)])
# plot_gf_turnovers(gf ,imps$var)

#### FORMING PREDICTIONS ####
# run model predictions with ordination jackknife
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
xy2 <- envc[goods, c("x", "y")]
names(xy2) <- c("lon", "lat")
rfpreds <- oj$predictions.direct
turnovers <- oj$predictions.turnover
bests <- names(oj$median.importance)

# create raster from predictions
rast_rfpreds <- rasterFromXYZ(cbind(xy2, rfpreds))
# Normalize rasters to 0-255 range
for (i in 1:3) {
  min_val <- rast_rfpreds[[i]]@data@min
  max_val <- rast_rfpreds[[i]]@data@max
  rast_rfpreds[[i]] <- (rast_rfpreds[[i]] - min_val) / (max_val - min_val) * 255
}
pdf("rdaforest_wp1_modern_sf7_sites_predictions.pdf", width = 8, height = 6)
plotRGB(rast_rfpreds, r = 1, g = 2, b = 3, bgalpha = 0)
points(latlon$lon, latlon$lat)
dev.off()

# plot turnovers
rast_turnovers <- rasterFromXYZ(cbind(xy2, turnovers))
# ... to be continued

## run rda forest using leave one out strategy
#for (i in seq_along(model_samples)){
for (i in 137:165){
  sample <- model_samples[i]
  print(sample)
  # remove the sample from the distance matrix
  cordist0 <- cordist[-which(rownames(cordist) == sample), -which(colnames(cordist) == sample)]
  # copy the full distance matrix placing the removed sample at the end
  cordist.i <- cordist
  cordist.i <- cordist.i[
    c(setdiff(rownames(cordist.i), sample), sample),
    c(setdiff(colnames(cordist.i), sample), sample)
  ]
  # build the landscape model without the lone wolf (will take a while...)
  oj0 <- ordinationJackknife(
    Y = cordist0,
    X = env[-i, c("sss_mean", "sst_mean")],
    newX = envc,
    nreps = 10,
    top.pcs = tokeep,
    extra = 0.1
  )
  # step one: getting gPCs for the other (model-building) wolves
  ord0 <- capscale(cordist0 ~ 1)
  sc0 <- scores(ord0, scaling = 1, display = "sites", choices = c(1:tokeep))
  # predicting gPCs for the same old wolves plus the new "lone wolf"
  scp <- predict(ord0, cordist.i, type = "sp", scaling = "sites")[, 1:tokeep]
  # making sure the predictions for the new wolf align (mostly in terms of gPC sign) with the model
  pro <- procrustes(sc0, scp[-nrow(cordist.i), ], scale = FALSE)
  # at last, getting the predicted gPCs for the lone wolf
  wi <- as.vector(as.numeric(pro$rotation %*% scp[nrow(scp), ]))[1:tokeep]
  sc <- adapt_scale(oj0$predictions.direct)[[2]]
  # calculating environmental mismatch
  home <- env_mismatch(X = wi, Y = oj0, sy = envc[, 1:2], sc = sc)
  mism <- rasterFromXYZ(home)
  sample_coords <- get_sample_latlon(sample)
  pdf(paste0("plots/rdaforest/env_mismatch/", model_samples_name, ".", sites_name, "/", sample, ".mismatch.pdf"))
  plot(mism, col = map.pal)
  points(sample_coords$lon, sample_coords$lat, pch = 0, col = "#3535ff", cex = 1.5, lwd = 2)
  text(sample_coords$lon, sample_coords$lat, labels = sample, pos = 1, cex = 1.1, col = "#3535ff")
  dev.off()
  # save mismatch raster
  writeRaster(mism,
  filename = paste0(
    "data/rdaforest/env_mismatch/", model_samples_name, ".", sites_name, "/", sample, ".mismatch.tif"
    ), overwrite = TRUE
  )
}

#### PREDICTIONS FOR REST OF SAMPLES ####
## run rda forest with the rest of samples (not used for model building)
rest_samples <- setdiff(all_samples, model_samples)
rest_samples <- rest_samples[rest_samples != "HER135"] # Exclude sample HER135
for (i in seq_along(rest_samples)) {
#for (i in 145:150) {
  sample <- rest_samples[i]
  print(sample)
  # build matrix with the new "lone wolf" added at the end
  cordist.i <- cordist
  sample_distances <- get_dist_matrix(matrix_file, c(model_samples, sample), all_samples)
  cordist.i <- sample_distances[
    c(model_samples, sample),
    c(model_samples, sample)
  ]
  # build the landscape model with the model samples only
  oj0 <- ordinationJackknife(
    Y = cordist,
    X = env[, c("sss_mean", "sst_mean")],
    newX = envc,
    nreps = 10,
    top.pcs = tokeep,
    extra = 0.1
  )
  # predicting gPCs for the model samples
  ord0 <- capscale(cordist ~ 1)
  sc0 <- scores(ord0, scaling = 1, display = "sites", choices = c(1:tokeep))
  # predicting gPCs for the same old wolves plus the new "lone wolf"
  scp <- predict(ord0, cordist.i, type = "sp", scaling = "sites")[, 1:tokeep]
  # making sure the predictions for the new wolf align (mostly in terms of gPC sign) with the model
  pro <- procrustes(sc0, scp[-nrow(cordist.i), ], scale = FALSE)
  # at last, getting the predicted gPCs for the lone wolf
  wi <- as.vector(as.numeric(pro$rotation %*% scp[nrow(scp), ]))[1:tokeep]
  sc <- adapt_scale(oj0$predictions.direct)[[2]]
  # calculating environmental mismatch
  home <- env_mismatch(X = wi, Y = oj0, sy = envc[, 1:2], sc = sc)
  mism <- rasterFromXYZ(home)
  sample_coords <- get_sample_latlon(sample)
  pdf(paste0("plots/rdaforest/env_mismatch/", all_samples_name, ".", sites_name, "/", sample, ".mismatch.pdf"))
  plot(mism, col = map.pal)
  points(sample_coords$lon, sample_coords$lat, pch = 0, col = "#3535ff", cex = 1.5, lwd = 2)
  text(sample_coords$lon, sample_coords$lat, labels = sample, pos = 1, cex = 1.1, col = "#3535ff")
  dev.off()
  # save mismatch raster
  writeRaster(mism,
    filename = paste0(
      "data/rdaforest/env_mismatch/", all_samples_name, ".", sites_name, "/", sample, ".mismatch.tif"
    ), overwrite = TRUE
  )
}
