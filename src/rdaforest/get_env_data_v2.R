library(sdmpredictors)
library(raster)
library(geobuffer)
library(tidyverse)

# dummy argument
samples_name <- "wp1_modern"

get_samples_from_file <- function(file_path) {
  samples <- read.table(file_path)[, 1] |> as.character()
  return(samples)
}

get_sample_coords <- function(data_table, sample_id) {
  coords <- data.frame(
    x = as.numeric(data_table[data_table$sample_id == sample_id, ]$x),
    y = as.numeric(data_table[data_table$sample_id == sample_id, ]$y)
  )
  return(coords)
}

get_sample_values <- function(sample_id, coords, t_layer, s_layer, buffer_dist = 50000) {
  buff <- geobuffer_pts(xy = as.matrix(coords), dist_m = buffer_dist)
  evalues <- data.frame(
    sample = sample_id,
    x = coords$x,
    y = coords$y,
    sst_mean = raster::extract(t_layer, buff, na.rm = TRUE, df = FALSE, fun = mean),
    sss_mean = raster::extract(s_layer, buff, na.rm = TRUE, df = FALSE, fun = mean)
  )
  return(evalues)
}

build_etable <- function(data_table, t_layer, s_layer) {
  etable <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(etable) <- c("sample", "x", "y", "sst_mean", "sss_mean")
  for (sample in data_table$sample_id){
    print(sample)
    coords <- get_sample_coords(data_table, sample)
    ecol <- get_sample_values(sample, coords, t_layer, s_layer)
    etable <- rbind(etable, ecol)
  }
  return(etable)
}

# load climate layers - temperature and salinity
t_layer <- load_layers("MS_biogeo13_sst_mean_5m")
s_layer <- load_layers("MS_biogeo08_sss_mean_5m")

# samples list file and load samples
samples_file <- paste0("data/angsd_matrix/bamlists/", samples_name, ".sample_list.txt")
samples <- get_samples_from_file(samples_file)

# load data table with coordinates and samples
data_table <- read.table("data/samples_table.csv", sep = ",", header = TRUE, na.strings = c("UNKNOWN", "NA"))
data_table <- data_table[data_table$sample_id %in% samples, ]

# build environmental data table
etable <- build_etable(data_table, t_layer, s_layer)

# save environmental data table
output_file <- paste0("data/rdaforest/etable.", samples_name, ".csv")
write.table(etable, output_file, sep = ",", row.names = FALSE, quote = FALSE)
