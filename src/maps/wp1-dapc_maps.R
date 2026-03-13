library(tidyverse)
library(sf)
library(cowplot)

# Function to create line from point to rectangle in direction of centroid
make_line_to_rect <- function(pt, rect) {
  if (is.null(rect) || length(rect) == 0) {
    return(st_geometrycollection())
  }
  
  # point coordinates
  p <- st_coordinates(pt)[1, 1:2]
  
  # rectangle centroid
  c <- st_coordinates(st_centroid(rect))[1, 1:2]
  
  # direction vector
  dir_vec <- c - p
  
  # extend line far enough
  bb <- st_bbox(rect)
  scale <- max(bb$xmax - bb$xmin, bb$ymax - bb$ymin) * 5
  far_pt <- p + dir_vec / sqrt(sum(dir_vec^2)) * scale
  long_line <- st_sfc(st_linestring(rbind(p, far_pt)), crs = st_crs(pt))
  
  # intersection with rectangle boundary
  hit <- st_intersection(long_line, st_boundary(rect))
  
  if (length(hit) == 0) {
    return(st_geometrycollection())
  } else {
    # compute distances from point to all intersections
    hit_coords <- st_coordinates(hit)[, 1:2, drop = FALSE]
    dists <- sqrt((hit_coords[,1] - p[1])^2 + (hit_coords[,2] - p[2])^2)
    closest <- hit_coords[which.min(dists), ]
    st_linestring(rbind(p, closest))
  }
}

prep_pie <- function(df, min_prop = 1e-4) {
  df |>
    mutate(prop = value / sum(value)) |>
    mutate(prop = ifelse(prop < min_prop, 0, prop)) |>
    filter(prop > 0)
}

# load world map shapefile
world_map_sf <- sf::read_sf("data/maps/ne_50m_land/ne_50m_land.shp")
world_map_sf <- st_transform(world_map_sf, crs = 4326)
# set limits
xl <- c(-27, 25)
yl <- c(47, 70)
# sample metadata
sample_data_file <- "data/samples_table.csv"
sample_data <- read.table(sample_data_file, sep = ",",
                          header = TRUE, na.strings = "UNKNOWN")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## PART 1 : DAPC pie charts for WP1 samples
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# load dapc_df and DAPC object
dapc_df <- read.table("data/angsd_matrix/dapc_df_wp1_final_sf7_sites.tsv",
                      header = TRUE, sep = "\t")
dapc_df$loc <- str_replace_all(dapc_df$loc, " ", "")

for (time_period in unique(dapc_df$time)){
  time_df <- dapc_df %>%
    filter(time == time_period)
  for (location in unique(time_df$loc)){

    a <- time_df |>
      gather(post_group, posterior, X1:X5) |>
      select(loc, post_group, posterior)

    P <- a |>
      filter(loc == location) |>
      group_by(post_group) |>
      summarize(Postsum = sum(posterior)) |>
      ungroup() |>
      mutate(Postprop = Postsum/sum(Postsum)) 

    PP <- data.frame(
      group = P$post_group,
      value = P$Postprop
    )

    PP2 <- PP |>
      mutate(prop = value / sum(value)) |>   # normalize
      mutate(prop = ifelse(prop < 1e-4, 0, prop)) |>  # kill tiny slices
      filter(prop > 0)
    
    PPP <- ggplot(PP2, aes(x = 1, y = prop, fill = group)) +
      geom_col(width = 0.1, color = "grey20", linewidth = 2) +
      coord_polar(theta = "y") +
      scale_fill_manual(values = c(
        "X1" = "#1b639e",
        "X5" = "#02d97c",
        "X3" = "#b370b2",
        "X2" = "#e73f29",
        "X4" = "#a6721e",
        "Other" = "#7a7a7a"
      )) +
      theme_void() +
      theme(legend.position = "none")
    
    ggsave(
      filename = paste0("plots/angsd_matrix/pie_charts/", time_period, ".", location, ".png"),
      plot = PPP,
      height = 3, width = 3,
      device = ragg::agg_png
    )
  }
}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## PART 2 : Maps with historical, sillperioder and modern pie charts
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# dataframe with x y coordinates for each location
modern_coords <- data.frame(
  loc = c("Iceland&FaroeIslands", "NorwegianSea&Fjords", "NorthSea", "Britain&Ireland", "Skagerrak-Kattegat"),
  x = c(-11,   15,  10,  -8,  18),
  y = c(64.5,  67,  54,  55, 58)
)
sillperioder_coords <- data.frame(
  loc = c("Iceland&FaroeIslands", "NorwegianSea&Fjords", "NorthSea", "Britain&Ireland", "Skagerrak-Kattegat"),
  x = c(-11,   0,  10, -11,  18),
  y = c(64.5,  65,  52,  55, 58)
)
historical_coords <- data.frame(
  loc = c("Iceland&FaroeIslands", "NorwegianSea&Fjords", "NorthSea", "Britain&Ireland", "Skagerrak-Kattegat"),
  x = c(-11,   0,  10,  -8,  18),
  y = c(64.5,  65,  52,  55, 58)
)

# sf version of the same
modern_points_sf <- modern_coords %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(world_map_sf))
sillperioder_points_sf <- sillperioder_coords %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(world_map_sf))
historical_points_sf <- historical_coords %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(world_map_sf))

# starting map
basemap <- ggplot() +
  geom_sf(data = world_map_sf, fill = "grey85", color = "grey0") +
  xlab("") +
  ylab("") +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

# PLOT MAP!

for (time_period in unique(dapc_df$time)){
  
  # time dataframe
  time_df <- dapc_df %>%
    filter(time == time_period)
  
  # metadata of samples in time period
  time_metadata <- sample_data %>%
    filter(sample_id %in% time_df$old_id)
  time_metadata$loc <- time_df$loc[match(time_metadata$sample_id, time_df$old_id)]
  
  # spatial coordinates of sampling location
  time_locations <- time_metadata %>%
    select(x, y, loc) %>%
    distinct() %>%
    st_as_sf(coords = c("x", "y"), crs = 4326) %>%
    st_transform(st_crs(world_map_sf))
  
  # create one rectangle per each loc with x min = xmin - 1, x max = xmax + 1, y min = ymin -1, y max = ymax + 1
  loc_rectangles <- time_locations %>%
    group_by(loc) %>%
    summarise(
      geometry = st_as_sfc(
        st_bbox(geometry) + c(-1, -0.3, 1, 0.3)
      ),
      .groups = "drop"
    )
  
  # create lines from pie charts to rectangles
  if (time_period == "modern"){
    points_sf <- modern_points_sf
    loc_coords <- modern_coords
  } else if (time_period == "sillperioder"){
    points_sf <- sillperioder_points_sf
    loc_coords <- sillperioder_coords
  } else if (time_period == "historical"){
    points_sf <- historical_points_sf
    loc_coords <- historical_coords
  }
  rect_df <- loc_rectangles %>%
    mutate(rect_geom = st_geometry(.)) %>%
    st_drop_geometry()
  points_joined <- points_sf %>%
    left_join(rect_df, by = "loc")
  lines_sfc <- mapply(
    make_line_to_rect,
    points_joined$geometry,
    points_joined$rect_geom,
    SIMPLIFY = FALSE
  )
  lines_to_edge <- st_sf(
    points_joined %>% select(-geometry),
    geometry = st_sfc(lines_sfc, crs = st_crs(points_sf))
  )
  
  # create new basemap
  basemap2 <- basemap +
    geom_sf(data = time_locations, color = "black", shape = 1, size = 2, stroke = 0.6) +
    geom_sf(data = loc_rectangles, fill = NA, color = "red", size = 0.3) +
    geom_sf(data = lines_to_edge,
            color = "red", linewidth = 0.3) +
    coord_sf(
      xlim = xl, ylim = yl, expand = FALSE, default_crs = NULL
    ) +
    # add title in the middle of the plot
    ggtitle(label = paste0(str_to_title(time_period))) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12)
    )
  
  for (location in unique(time_df$loc)){
    # add pie charts
    filename <- paste0("plots/angsd_matrix/pie_charts/", time_period, ".", location, ".png")
    popx <- loc_coords$x[loc_coords$loc == location]
    popy <- loc_coords$y[loc_coords$loc == location]
    basemap2 <- basemap2 +
      draw_image(filename, x = popx, y = popy, vjust = 0.5, hjust = 0.5, scale = 8) 
  }
  ggsave(
    filename = paste0("plots/angsd_matrix/", time_period, ".dapc_map.png"),
    plot = basemap2,
    height = 5, width = 5
  )
  ggsave(
    filename = paste0("plots/angsd_matrix/", time_period, ".dapc_map.pdf"),
    plot = basemap2,
    height = 5, width = 5
  )
 }

