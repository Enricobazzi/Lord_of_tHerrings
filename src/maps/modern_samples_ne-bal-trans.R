library(tidyverse)
library(sf)
library(ggrepel)

plot_map <- function(world_map_sp){
  # Load shapefile as sp object
  # Convert to sf object
  world_map_sf <- st_as_sf(world_map_sp)
  # Plot
  wm <- ggplot() +
    geom_sf(data = world_map_sf, fill = "grey80", color = "grey0") +
    coord_sf(
      xlim = c(-80, 30),
      ylim = c(40, 75),
      expand = FALSE,
      default_crs = st_crs(world_map_sf)
    )
  return(wm)
}

tab <- read.table("data/samples_table.csv", sep = ",", header = TRUE, na.strings = "UNKNOWN")
world_map_sp <- rgdal::readOGR("~/Documents/Lynxtrogression_v2/plots/ms_figures/ne_50m_land/ne_50m_land.shp")

modern <- tab %>%
  filter(
    species == "clupea_harengus",
    (region != "E-PACIFIC" & region != "NW-ATL") | is.na(region),
    year >= 1930
  ) %>%
  mutate(
    x = as.numeric(x),
    y = as.numeric(y)
  )

modern_summary <- modern %>%
  count(x, y, spawn)

plot_map(world_map_sp) +
  geom_point(
    data = modern_summary,
    aes(x = x, y = y, fill = spawn, size = n), 
    shape = 21) +
  geom_label_repel(
    data = modern_summary,
    aes(x = x, y = y, label = n),
    size = 3,
    box.padding = 0.3,
    label.size = 0.2,       # border thickness
    fill = "white",         # background color
    color = "black"         # text color
  ) +
  xlab("") + ylab("")

# writeLines(modern$sample_id, "data/rdaforest/modern_samples.txt")

s <- c("AAL1", "AAL2", "AAL3", "AF29", "AF30", "AF31", "AF8", "AK1", "AK2", "AK3", "AM27", "AM29", "AM33", "AM8", "BF16", "BF18", "BF19", "BF21", "BM14", "BM15", "BM16", "BM19", "Fehmarn3", "Fehmarn44", "Fehmarn6", "Gavle100", "Gavle54", "Gavle98", "MHER001", "MHER003", "MHER006", "MHER007", "MHER008", "MHER009", "MHER010", "MHER011", "MHER012", "MHER013", "MHER014", "MHER015", "MHER016", "MHER017", "MHER018", "MHER019", "MHER020", "MHER021", "MHER022", "MHER023", "MHER024", "MHER025", "MHER026", "MHER027", "MHER028", "MHER029", "MHER030", "MHER034", "MHER035", "MHER036", "MHER037", "MHER038", "MHER039", "MHER044", "MHER045", "MHER046", "MHER052", "MHER053", "MHER054", "MHER055", "MHER056", "MHER057", "MHER058", "MHER059", "MHER061", "MHER062", "MHER063", "MHER065", "MHER066", "NorthSea13", "NorthSea19", "NorthSea34", "NSSH33", "NSSH34", "NSSH36", "Z12", "Z14", "Z4")
s <- c("201550141", "201550145", "201604591", "201650794", "201650795", "201650797", "201650798", "201650799", "1552004723", "2016507910", "2016507911", "2016507912", "2016507913", "2016507915", "2016507916", "2016507918", "2016507919", "2016507920", "AAL1", "AAL2", "AAL3", "AF29", "AF30", "AF31", "AF8", "AK1", "AK2", "AK3", "AM27", "AM29", "AM33", "AM8", "BF16", "BF18", "BF19", "BF21", "BM14", "BM15", "BM16", "BM19", "Fehmarn3", "Fehmarn44", "Fehmarn6", "Gavle100", "Gavle54", "Gavle98", "MHER001", "MHER002", "MHER003", "MHER005", "MHER006", "MHER007", "MHER008", "MHER009", "MHER010", "MHER011", "MHER012", "MHER013", "MHER014", "MHER015", "MHER016", "MHER017", "MHER018", "MHER019", "MHER020", "MHER021", "MHER022", "MHER023", "MHER024", "MHER025", "MHER026", "MHER027", "MHER028", "MHER029", "MHER030", "MHER034", "MHER035", "MHER036", "MHER037", "MHER038", "MHER039", "MHER044", "MHER045", "MHER046", "MHER052", "MHER053", "MHER054", "MHER055", "MHER056", "MHER057", "MHER058", "MHER059", "MHER061", "MHER062", "MHER063", "MHER064", "MHER065", "MHER066", "NorthSea13", "NorthSea19", "NorthSea34", "NSSH33", "NSSH34", "NSSH36", "Z12", "Z14", "Z4")
modern <- modern[modern$sample_id %in% s, ]



