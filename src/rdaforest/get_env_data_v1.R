library(sdmpredictors)
library(raster)
library(geobuffer)
library(tidyverse)

# load climate layers - temperature and salinity
t_layer <- load_layers("MS_biogeo13_sst_mean_5m")
s_layer <- load_layers("MS_biogeo08_sss_mean_5m")

# load data table with coordinates and samples (manually selected for now)
data_table <- read.table("data/samples_table.csv", sep = ",", header = TRUE, na.strings = "UNKNOWN")
samples <- c("AAL1", "AAL2", "AAL3", "AF29", "AF30", "AF31", "AF8", "AK1", "AK2", "AK3", "AM27", "AM29", "AM33", "AM8", "BF16", "BF18", "BF19", "BF21", "BM14", "BM15", "BM16", "BM19", "Fehmarn3", "Fehmarn44", "Fehmarn6", "Gavle100", "Gavle54", "Gavle98", "MB1", "MB2", "MB3", "MHER001", "MHER003", "MHER006", "MHER007", "MHER008", "MHER009", "MHER010", "MHER011", "MHER012", "MHER013", "MHER014", "MHER015", "MHER016", "MHER017", "MHER018", "MHER019", "MHER020", "MHER021", "MHER022", "MHER023", "MHER024", "MHER025", "MHER026", "MHER027", "MHER028", "MHER029", "MHER030", "MHER034", "MHER035", "MHER036", "MHER037", "MHER038", "MHER039", "MHER044", "MHER045", "MHER046", "MHER052", "MHER053", "MHER054", "MHER055", "MHER056", "MHER057", "MHER058", "MHER059", "MHER061", "MHER062", "MHER063", "MHER065", "MHER066", "NorthSea13", "NorthSea19", "NorthSea34", "NSSH33", "NSSH34", "NSSH36", "Z12", "Z14", "Z4")
data_table <- data_table[data_table$sample_id %in% samples, ]

# empty table to fill with data from each sample
etable <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(etable) <- c("sample", "sst_mean", "sss_mean")
for (sample in samples){
  coords <- data.frame(
    x = as.numeric(data_table[data_table$sample_id == sample, ]$x),
    y = as.numeric(data_table[data_table$sample_id == sample, ]$y)
  )
  buff <- geobuffer_pts(xy = as.matrix(coords), dist_m = 50000)
  ecol <- data.frame(
    sample = sample,
    sst_mean = raster::extract(t_layer, buff, na.rm = TRUE, df = FALSE, fun = mean),
    sss_mean = raster::extract(s_layer, buff, na.rm = TRUE, df = FALSE, fun = mean)
  )
  etable <- rbind(etable, ecol)
}
