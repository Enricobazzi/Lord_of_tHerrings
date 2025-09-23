library(sdmpredictors)
library(raster)
library(geobuffer)
library(tidyverse)

# load climate layers - temperature and salinity
t_layer <- load_layers("MS_biogeo13_sst_mean_5m")
s_layer <- load_layers("MS_biogeo08_sss_mean_5m")

# load data table with coordinates
data_table <- read.table("data/selection_scans_poolseq/poolseqdata_info.csv",
                         sep = ",", header = TRUE)
# remove pacific herring
pacific_pops <- c("HWS1_Japan_SeaOfJapan", "HWS2_PechoraSea_BarentsSea",
                 "HWS3_WhiteSea_WhiteSea", "HWS4_KandalakshaBay_WhiteSea",
                 "HWS5_KandalakshaBay_WhiteSea", "HWS6_Balsfjord_Atlantic",
                 "PB8_Pacific_Pacific_Spring")
data_table <- data_table %>% filter(!Sample.name %in% pacific_pops)

# empty table to fill with data from each location
locations <- unique(data_table$Sample.name)
etable <- data.frame(matrix(nrow = 2, ncol = 0))
rownames(etable) <- c("sst_mean", "sss_mean")

# extract data around each point
for (location in locations){
  coords <- data_table %>% filter(Sample.name == location) %>%
    mutate(x = Longitude, y = Latitude) %>%
    select(x, y) %>% distinct()
  buff <- geobuffer_pts(xy = coords, dist_m = 50000)
  ecol <- data.frame(
    location = location,
    sst_mean = raster::extract(t_layer, buff, na.rm = T, df = F, fun = mean),
    sss_mean = raster::extract(s_layer, buff, na.rm = T, df = F, fun = mean)
    ) %>%
    pivot_longer(cols = -location, names_to = "variable", values_to = "value") %>%
    pivot_wider(names_from = location, values_from = value) %>%
    column_to_rownames("variable")
  etable <- cbind(etable, ecol)
}

# save tables in format for baypass (no column and row names)
write.table(etable[1,], file = "data/selection_scans_poolseq/sst_mean.txt",
            col.names = F, row.names = F, sep = " ")
write.table(etable[2,], file = "data/selection_scans_poolseq/sss_mean.txt",
            col.names = F, row.names = F, sep = " ")