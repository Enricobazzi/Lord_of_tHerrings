library(sdmpredictors)
library(raster)
library(geobuffer)
library(tidyverse)

# list of locations
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
etable <- data.frame(matrix(nrow = 2, ncol = 0))
rownames(etable) <- c("sst_mean", "sss_mean")

# extract data around each point
for (location in locations){
  coords <- data_table %>%
    filter(Sample.name == location) %>%
    mutate(x = Longitude, y = Latitude) %>%
    dplyr::select(x, y) %>%
    distinct()
  buff <- geobuffer_pts(xy = coords, dist_m = 50000)
  ecol <- data.frame(
    location = location,
    sst_mean = raster::extract(t_layer, buff, na.rm = TRUE, df = FALSE, fun = mean),
    sss_mean = raster::extract(s_layer, buff, na.rm = TRUE, df = FALSE, fun = mean)
    ) %>%
    pivot_longer(cols = -location, names_to = "variable", values_to = "value") %>%
    pivot_wider(names_from = location, values_from = value) %>%
    column_to_rownames("variable")
  etable <- cbind(etable, ecol)
}

# save tables in format for baypass (no column and row names)
write.table(etable[1, ], file = "data/selection_scans_poolseq/sst_mean.txt",
            col.names = FALSE, row.names = FALSE, sep = " ")
write.table(etable[2, ], file = "data/selection_scans_poolseq/sss_mean.txt",
            col.names = FALSE, row.names = FALSE, sep = " ")
