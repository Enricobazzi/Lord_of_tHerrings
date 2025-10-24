library(tidyverse)
library(sf)

# function to plot basic map
plot_map <- function(world_map_sp){
  # Convert to sf object
  world_map_sf <- st_as_sf(world_map_sp)
  # Plot
  wm <- ggplot() +
    geom_sf(data = world_map_sf, fill = "grey80", color = "grey0") +
    coord_sf(
      xlim = c(-80, 30),
      ylim = c(40, 72),
      expand = FALSE,
      default_crs = st_crs(world_map_sf)
    )
  return(wm)
}

# list of locations
data_table <- read.table("data/selection_scans_poolseq/poolseqdata_info.csv",
                         sep = ",", header = TRUE)

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

# Load shapefile as sp object
world_map_sp <- rgdal::readOGR("~/Documents/Lynxtrogression_v2/plots/ms_figures/ne_50m_land/ne_50m_land.shp")

for (n in 1:500){
  nn <- sprintf("%03d", n)
  print(nn)
  # read omega matrix
  mat_omega <- as.matrix(read.table(paste0(
    "data/selection_scans_poolseq/baypass_core/",
    nn,
    "_mat_omega.out"
    )))
  rownames(mat_omega) <- locations
  colnames(mat_omega) <- locations
  # dendrogram from matrix
  dist_mat <- as.dist(1 - abs(cov2cor(mat_omega)))  # convert covariance to distance
  hc <- hclust(dist_mat)
  pdf(paste0("plots/selection_scans_poolseq/baypass_core/", nn,"_omegamat.pdf"),
      width = 12, height = 8)
  plot(hc)
  dev.off()
  # plot points on map colored by group
  groups <- data.frame(group=cutree(hc, k=10)) %>% rownames_to_column(var="location")
  groups$group <- as_factor(groups$group)
  coords <- data.frame()
  for (location in groups$location){
    coord <- data_table %>%
      filter(Sample.name == location) %>%
      mutate(x = Longitude, y = Latitude) %>%
      dplyr::select(x, y) %>%
      distinct()
    coords <- rbind(coords, coord)
  }
  groups <- cbind(groups, coords)
  map <- plot_map(world_map_sp) +
    geom_point(data = groups, aes(x, y, fill = group), shape = 21, size = 3) + 
    scale_fill_discrete()
  ggsave(paste0("plots/selection_scans_poolseq/baypass_core/", nn,"_hcmap.pdf"),
         plot = map, width = 12, height = 8)
}

###



####

library(FactoMineR)
library(tidyverse)
locations <- strsplit(readLines("data/published_data/60.Neff.freq", n = 1), "\t")[[1]][-c(1:2)]

a <- read.table("data/selection_scans_poolseq/subsets_freqs/60.Neff.001.freq")[, -(1:2)]
colnames(a) <- locations
pacific_pops <- c("HWS1_Japan_SeaOfJapan", "HWS2_PechoraSea_BarentsSea",
                   "HWS3_WhiteSea_WhiteSea", "HWS4_KandalakshaBay_WhiteSea",
                   "HWS5_KandalakshaBay_WhiteSea", "HWS6_Balsfjord_Atlantic",
                   "PB8_Pacific_Pacific_Spring")
a <- a[ , !(names(a) %in% pacific_pops)]
ab <- t(na.omit(a))
pc_ab <- FactoMineR::PCA(ab, graph = F, scale.unit=F)
pcdf <- as.data.frame(pc_ab$ind$coord)
cols <- c()
for (loc in rownames(pcdf)){
  ifelse(grepl("Baltic", loc), cols <- c(cols, "blue"),
         ifelse(grepl("Atlantic", loc), cols <- c(cols, "red"), cols <- c(cols, "green")))
}
ggplot() +
  geom_point(data = pcdf, aes(x=Dim.1, y=Dim.2), shape=21, fill=cols)

aaa <- merge(groups, pcdf %>% rownames_to_column(var="location"), by = "location")
ggplot() +
  geom_point(data = aaa, aes(x=Dim.1, y=Dim.2, fill=group), shape=21, size = 3) + 
  scale_fill_discrete()
