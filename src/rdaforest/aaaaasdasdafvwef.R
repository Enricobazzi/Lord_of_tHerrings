library(biscale)
library(sf)
library(tidyverse)
library(cowplot)

plot(turnovers$sst_mean, turnovers$sss_mean)
plot(env$sst_mean, env$sss_mean)

rast_turnovers <- rasterFromXYZ(cbind(xy2, turnovers))
world_map_sp <- rgdal::readOGR("~/Documents/Lynxtrogression_v2/plots/ms_figures/ne_50m_land/ne_50m_land.shp")

crs(rast_turnovers) <- crs(world_map_sp)
rast_turnovers <- terra::project(rast(rast_turnovers), "EPSG:6933")

world_map_sf <- st_as_sf(world_map_sp)
world_map_sf <- st_transform(world_map_sf, 6933)

# Convert your raster brick to terra object (if it's a raster::RasterBrick)
#r_terra <- rast(aggregate(rast_turnovers, fact = 2, fun = mean))  # or rast(r_agg) if you aggregated earlier
r_terra <- rast_turnovers

# Create hexagonal grid (in the same CRS as raster)
hex_sf <- st_make_grid(st_as_sfc(st_bbox(r_terra)),
                       cellsize = 45681.84,  # e.g., 50 km (adjust as needed)
                       square = FALSE) |> 
  st_sf()
# Convert sf hex grid to terra vector
plot(r_terra)
plot(hex_vect)
hex_vect <- vect(hex_sf)
# Extract mean raster values for each hex
hex_mean <- terra::extract(r_terra, hex_vect, fun = mean, na.rm = TRUE)

# Join the results back to the hex polygons
hex_sf$sst_mean <- hex_mean[, 2]  # second column is the mean value
hex_sf$sss_mean <- hex_mean[, 3]  # second column is the mean value
ndimes <- 4
data <- bi_class(na.omit(hex_sf), x = sst_mean, y = sss_mean, style = "fisher", dim = ndimes)

a <- ggplot() +
  geom_sf(data = data, aes(fill = bi_class), color = "white", size = 0.005, show.legend = FALSE) +
  bi_scale_fill(pal = "DkViolet2", dim = 4, rotate_pal = FALSE, flip_axes = FALSE) + 
  geom_sf(data = world_map_sf, fill = "grey80", color = "grey0") +
  bi_theme() +
  coord_sf(
    xlim = c(st_bbox(r_terra)[1], st_bbox(r_terra)[3]),
    ylim = c(st_bbox(r_terra)[2], st_bbox(r_terra)[4]),
    expand = FALSE,
    default_crs = st_crs(world_map_sf)
  )
a
b <- bi_legend(pal = "DkViolet2",
               xlab = "Temperature",
               ylab = "Salinity",
               size = 9, dim = ndimes,
               rotate_pal = FALSE, flip_axes = FALSE)
b
finalPlot <- ggdraw() +
  draw_plot(a, 0, 0, 1, 1) +
  draw_plot(b, 0.2, .65, 0.2, 0.2)
finalPlot
