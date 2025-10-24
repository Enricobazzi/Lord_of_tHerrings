library(tidyverse)
library(raster)
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

# world map
world_map_sp <- rgdal::readOGR("~/Documents/Lynxtrogression_v2/plots/ms_figures/ne_50m_land/ne_50m_land.shp")

# samples table
samples_table <- read.table("data/samples_table.csv", sep = ",", header = TRUE)
wg_samples <- samples_table[!is.na(samples_table$wg.depth), ]
wg_samples$year[wg_samples$year == "UNKNOWN" & grepl("^ND", wg_samples$sample_id)] <- 1850
wg_samples$year[wg_samples$year == "UNKNOWN" & grepl("^MB", wg_samples$sample_id)] <- 2000
wg_samples$x[wg_samples$x == "UNKNOWN"] <- NA
wg_samples$y[wg_samples$y == "UNKNOWN"] <- NA
wg_samples$year <- as.numeric(wg_samples$year)
wg_samples$x <- as.numeric(wg_samples$x)
wg_samples$y <- as.numeric(wg_samples$y)

# add century category
wg_samples$century <- ifelse((wg_samples$year >= 600 & wg_samples$year) < 700, "7th",
                             ifelse((wg_samples$year >= 700 & wg_samples$year) < 800, "8th",
                             ifelse((wg_samples$year >= 800 & wg_samples$year) < 900, "9th",
                             ifelse((wg_samples$year >= 900 & wg_samples$year) < 1000, "10th",
                             ifelse((wg_samples$year >= 1000 & wg_samples$year) < 1100, "11th",
                             ifelse((wg_samples$year >= 1100 & wg_samples$year) < 1200, "12th",
                             ifelse((wg_samples$year >= 1200 & wg_samples$year) < 1300, "13th",
                             ifelse((wg_samples$year >= 1200 & wg_samples$year) < 1300, "13th",
                             ifelse((wg_samples$year >= 1200 & wg_samples$year) < 1300, "13th",
                             ifelse((wg_samples$year >= 1300 & wg_samples$year) < 1400, "14th",
                             ifelse((wg_samples$year >= 1400 & wg_samples$year) < 1500, "15th",
                             ifelse((wg_samples$year >= 1500 & wg_samples$year) < 1600, "16th", 
                             ifelse((wg_samples$year >= 1600 & wg_samples$year) < 1700, "17th",
                             ifelse((wg_samples$year >= 1700 & wg_samples$year) < 1800, "18th",
                             ifelse((wg_samples$year >= 1800 & wg_samples$year) < 1900, "19th",
                             ifelse((wg_samples$year >= 1900 & wg_samples$year) < 2000, "20th",
                             ifelse((wg_samples$year >= 2000 & wg_samples$year) < 2100, "21st",
                                    "none")))))))))))))))))
plot_map(world_map_sp) +
  geom_point(data = wg_samples[grepl("^ND", wg_samples$sample_id), ],
             aes(x, y, fill = year), shape = 21, size = 3) + 
  scale_fill_continuous()
