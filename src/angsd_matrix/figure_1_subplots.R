library(tidyverse)
library(ggforce)
library(adegenet)
library(cowplot)
library(sf)

## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## FUNCTIONS:

# get sample names for a dataset from the sample_list file:
get_samples_from_dataset <- function(dataset) {
  file_path <- paste0("data/angsd_matrix/bamlists/", dataset, ".sample_list.txt")
  samples <- read.table(file_path)[, 1] |> as.character()
  return (samples)
}

# get metadata table
get_metadata <- function(dapc_dataset) {
  samples <- get_samples_from_dataset(dapc_dataset)
  sample_data_file <- "~/Documents/Silly-periods/data/samples_table.csv"
  sample_data <- read.table(sample_data_file, sep = ",",
                            header = TRUE, na.strings = "UNKNOWN")
  sample_data <- sample_data[sample_data$sample_id %in% samples, ]
  return (sample_data)
}

# read best k from table
read_best_k <- function(dapc_dataset, sites_name) {
  best_k <- read.table(
    paste0("data/angsd_matrix/dapc/", dapc_dataset, ".", sites_name, ".k_table.txt"),
    header = T)[1,1]
  return (best_k)
}

# load dapc object
load_dapc <- function(pcangsd_dataset, dapc_dataset, sites_name, k) {
  DAPC <- readRDS(paste0("data/angsd_matrix/dapc/final_dapc.",
                         pcangsd_dataset, ".", dapc_dataset, ".", sites_name, ".k", k, ".rds"))
  return (DAPC)
}

# get the dictionary appropriate for the sites
get_groups_dict <- function(sites_name){
  # manually selected by inspecting the figure from plot_dapc_v2.R
  if (sites_name == "supplementary_file_7.v2") {
    groups_dict <- c(
      "1" = "North Atlantic",
      "2" = "North Sea",
      "3" = "Norwegian",
      "4" = "Britain & Ireland",
      "5" = "Skagerrak & Kattegat",
      "6" = "Baltic"
    )
  } else if (sites_name == "baltic_v_atlantic.v2") {
    groups_dict <- c(
      "1" = "Marine",
      "2" = "Brackish",
      "3" = "Intermediate"
    )
  } else if (sites_name == "spring_v_autumn.v2") {
    groups_dict <- c(
      "1" = "Autumn",
      "2" = "Spring",
      "3" = "Spring",
      "4" = "Spring"
    )
  }
}

# stock assignment
assign_ecotype <- function(GeneticCluster, SpawnTime) {
  
  if (
    GeneticCluster == "North Atlantic" &
    SpawnTime == "Autumn"
  ) {
    ecotype <- "North Atlantic Autumn-Spawner"
  } else if (
    GeneticCluster == "Norwegian" &
    SpawnTime == "Spring"
  ) {
    ecotype <- "Norwegian Spring-Spawner"
  } else if (
    GeneticCluster == "North Sea" &
    SpawnTime == "Autumn"
  ) {
    ecotype <- "North Sea Autumn-Spawner"
  } else if (
    GeneticCluster == "Skagerrak & Kattegat" &
    SpawnTime == "Spring"
  ) {
    ecotype <- "Skagerrak & Kattegat Spring-Spawner"
  } else if (
    GeneticCluster == "Baltic" &
    SpawnTime == "Spring"
  ) {
    ecotype <- "Baltic Spring-Spawner"
  } else if (
    GeneticCluster == "Britain & Ireland" &
    SpawnTime == "Autumn"
  ) {
    ecotype <- "Ireland & Britain Autumn-Spawner"
  } else if (
    GeneticCluster == "Baltic" &
    SpawnTime == "Autumn"
  ) {
    ecotype <- "Baltic Autumn-Spawner"
  } else if (
    GeneticCluster == "Britain & Ireland" &
    SpawnTime == "Spring"
  ) {
    ecotype <- "Ireland & Britain Spring-Spawner"
  } else if (
    GeneticCluster == "Skagerrak & Kattegat" &
    SpawnTime == "Autumn"
  ) {
    ecotype <- "Skagerrak & Kattegat Autumn-Spawner"
  } else if (
    GeneticCluster == "North Atlantic" &
    SpawnTime == "Spring"
  ) {
    ecotype <- "North Atlantic Spring-Spawner"
  } else if (
    GeneticCluster == "North Sea" &
    SpawnTime == "Spring"
  ) {
    ecotype <- "North Sea Spring-Spawner"
  } else {
    ecotype <- "Other"
  }
  return (ecotype)
}

# stock assignment in table
generate_ecotype_column <- function(eco_df) {
  mapply(
    function(x1, x2) {
      # Example rule:
      assign_ecotype(x1, x2)
    },
    eco_df$GeneticGroup, eco_df$SpawnTime,
    USE.NAMES = FALSE
  )
}

# read NS inv table
read_ns_column <- function(pcangsd_dataset, dapc_dataset){
  a <- read.table(
    paste0("data/angsd_matrix/dapc/ns_inv_table.", pcangsd_dataset, ".", dapc_dataset, ".csv"),
    sep = ",", header = T)
  return (a$Inversion)
}

# get fancy period names
get_fancy_period <- function(period) {
  period_dict <- c(
    "ah" = "Before 1700s",
    "17sp" = "1747–1805 Sillperiod",
    "18rh" = "Between Sillperiods",
    "18sp" = "1877–1906 Sillperiod",
    "mh" = "Present-day"
  )
  return (unname(period_dict[period]))
}

# get fancy region names
get_fancy_region <- function(region) {
  return(gsub("_", " ", region))
}

# get group name
get_group <- function(group, sites_name) {
  groups_dict <- get_groups_dict(sites_name)
  return (unname(groups_dict[group]))
}

# get groups from dataset and sites
grp_assigns <- function(pcangsd_dataset, dapc_dataset, sites_name) {
  best_k <- read_best_k(dapc_dataset, sites_name)
  dapc_obj <- load_dapc(pcangsd_dataset, dapc_dataset, sites_name, best_k)
  return (sapply(dapc_obj$assign, get_group, sites_name = sites_name))
}

# stock palette of colors
stock_palette <- c(
  "Norwegian Spring-Spawner" = "#02d97c",
  "North Atlantic Autumn-Spawner" = "#fed439",
  "North Atlantic Spring-Spawner" = "#fff5c8",
  "North Sea Autumn-Spawner" = "#ff7f00",
  "North Sea Spring-Spawner" = "#ffd199",
  "Ireland & Britain Spring-Spawner" = "#a6916e",
  "Ireland & Britain Autumn-Spawner" = "#a6721e",
  "Baltic Spring-Spawner" = "#9ceee7",
  "Baltic Autumn-Spawner" = "#00b8a9",
  "Skagerrak & Kattegat Spring-Spawner" = "#a55ba4",
  "Skagerrak & Kattegat Autumn-Spawner" = "#e0c7e0"
)

# doughnut plot
plot_donut <- function(df, period, region){
  
  zz <- df |> filter(Period == period, Region == region)
  
  plot_table <- zz |>
    count(Stock, name = "Count") |>
    mutate(
      fraction = Count / sum(Count),
      ymax = cumsum(fraction) * 2 * pi,
      ymin = lag(ymax, default = 0)
    )
  
  p <- ggplot(plot_table) +
    # annotate("point", x = 0, y = 0, shape = 19,
    #          size = 1, fill = "white", color = "white") +
    geom_arc_bar(aes(x0 = 0, y0 = 0,
                     r0 = 0.00,        # inner radius: smaller = bigger hole, larger = smaller hole
                     r  = 1.00,        # outer radius
                     start = ymin, end = ymax, fill = Stock),
                 color = "grey20", linewidth = 0.05) +
    # annotate("text", x = 0, y = 0,
    #          label = sum(nrow(zz)),
    #          fontface = "bold", family = "Verdana", size = 0.5) +
    scale_fill_manual(values = stock_palette) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none")
  return (p)
}


# NS allele frequencies 
plot_ns_freq <- function(df, period, region){
  
  zz <- z |> filter(Period == period, Region == region)
  # Count alleles
  alleles <- zz %>%
    summarise(
      North = 2 * sum(Inversion == "North", na.rm = TRUE) +
        1 * sum(Inversion == "Het", na.rm = TRUE),
      South = 2 * sum(Inversion == "South", na.rm = TRUE) +
        1 * sum(Inversion == "Het", na.rm = TRUE)
    ) %>%
    tidyr::pivot_longer(everything(),
                        names_to = "Allele",
                        values_to = "Count") %>%
    mutate(
      frac = Count / sum(Count),
      end = 2 * pi * cumsum(frac),
      start = lag(end, default = 0)
    )
  
  p <- ggplot(alleles) +
    geom_arc_bar(
      aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
        start = start, end = end, fill = Allele), color = "grey20", linewidth = 0.05
    ) +
    coord_fixed() +
    scale_fill_manual(values = c("North" = "blue3", "South" = "red3")) +
    theme_void() +
    theme(legend.position = "none")
  return (p)
}

# pie chart of salinity
plot_salinity_pie <- function(df, period, region){
  
  zz <- df |> filter(Period == period, Region == region)
  
  plot_table <- zz |>
    count(Salinity, name = "Count") |>
    mutate(
      fraction = Count / sum(Count),
      ymax = cumsum(fraction) * 2 * pi,
      ymin = lag(ymax, default = 0)
    )
  
  p <- ggplot(plot_table) +
    geom_arc_bar(aes(x0 = 0, y0 = 0,
                     r0 = 0.00,        # inner radius: smaller = bigger hole, larger = smaller hole
                     r  = 1.00,        # outer radius
                     start = ymin, end = ymax, fill = Salinity),
                 color = "grey20", linewidth = 0.05) +
    scale_fill_manual(values = c("Marine" = "white", "Intermediate" = "grey",
                                 "Brackish" = "black")) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none")
  return (p)
}

# plot N
plot_n <- function(df, period, region){
  zz <- df |> filter(Period == period, Region == region)
  
  p <- ggplot(data.frame(x = 0, y = 0),
              aes(x, y)) +
    geom_point(shape = 21, color = "grey20", fill = "white", stroke = 0.05, size = 3) +
    annotate("text", x = 0, y = 0,
             label = sum(nrow(zz)),
             fontface = "bold", family = "Verdana", size = 1) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none")
  
  return (p)
}

plot_map <- function(df, period){
  # subeset
  zz <- df |> filter(Period == period,
                     (Region == "Skagerrak & Kattegat" | Region == "Norway" | Region == "Britain & Ireland"))
  
  # load world map shapefile
  world_map_sf <- sf::read_sf("data/maps/ne_50m_land/ne_50m_land.shp")
  world_map_sf <- st_transform(world_map_sf, crs = 4326)
  # set limits
  xl <- c(-15, 22)
  yl <- c(50, 70)
  # point sf from coords
  points_sf <- st_as_sf(
    zz,
    coords = c("x", "y"),
    crs = 4326,
    remove = FALSE
  )
  points_sf <- points_sf[!duplicated(st_geometry(points_sf)), ]
  
  # make bbox polygon in EPSG:4326
  bbox_4326 <- st_as_sfc(
    st_bbox(
      c(xmin = xl[1], xmax = xl[2],
        ymin = yl[1], ymax = yl[2]),
      crs = st_crs(4326)
    )
  )
  
  # transform everything to EPSG:3035
  world_map_3035 <- st_transform(world_map_sf, 3035)
  points_3035 <- st_transform(points_sf, 3035)
  bbox_3035 <- st_transform(bbox_4326, 3035)
  
  # starting map
  basemap <- ggplot() +
    geom_sf(data = world_map_3035, fill = "#E3E3E3", color = "grey0", linewidth = 0.1) +
    # geom_sf(data = bbox_3035, fill = NA, color = "red", linewidth = 1) +
    geom_sf(data = points_3035, aes(shape = Region), size = 1.5, stroke = 0.3, show.legend = F) +
    coord_sf(crs = st_crs(3035), xlim = st_bbox(bbox_3035)[c("xmin", "xmax")],
             ylim = st_bbox(bbox_3035)[c("ymin", "ymax")]) +
    scale_shape_manual(values = c(
      "Britain & Ireland" = 15,
      "Skagerrak & Kattegat" = 19,
      "Norway" = 17
    )) +
    xlab("") +
    ylab("") +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "#F7F7F7"),
      plot.background = element_blank(),
      panel.border = element_rect(color = "#A8A8A8", fill = NA, linewidth = 0.5),
      plot.title = element_text(hjust = 0.5, size = 12),
      plot.margin = margin(10, 5, 0, 0),
      axis.title = element_text(size = 6, face = "bold"),
      axis.text = element_text(size = 6),
      axis.ticks = element_line(linewidth = 0.1)
    )
  return (basemap)
}

plot_stocks_legend <- function(z) {
  p <- ggplot(
    data = z |> filter(Region == "Skagerrak & Kattegat" | Region == "Norway" | Region == "Britain & Ireland"),
    aes(x, y, fill = Stock)) +
    geom_point(shape = 21, size = 5) +
    scale_fill_manual(values = stock_palette) +
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  legend <- ggdraw(get_legend(p))
  return (legend)
}

plot_ns_legend <- function(z) {
  
  zz <- z |> filter(Period == "Present-day", Region == "Skagerrak & Kattegat")
  # Count alleles
  alleles <- zz %>%
    summarise(
      North = 2 * sum(Inversion == "North", na.rm = TRUE) +
        1 * sum(Inversion == "Het", na.rm = TRUE),
      South = 2 * sum(Inversion == "South", na.rm = TRUE) +
        1 * sum(Inversion == "Het", na.rm = TRUE)
    ) %>%
    tidyr::pivot_longer(everything(),
                        names_to = "Allele",
                        values_to = "Count") %>%
    mutate(
      frac = Count / sum(Count),
      end = 2 * pi * cumsum(frac),
      start = lag(end, default = 0)
    )
  
  p <- ggplot(alleles) +
    geom_arc_bar(
      aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
          start = start, end = end, fill = Allele), color = "grey20", linewidth = 0.05
    ) +
    coord_fixed() +
    scale_fill_manual(values = c("North" = "blue3", "South" = "red3")) +
    theme_void()
  
  legend <- ggdraw(get_legend(p))
  return (legend)
}

plot_salinity_legend <- function(z) {
  
  zz <- z |> filter(Period == "Present-day", Region == "Skagerrak & Kattegat")
  
  plot_table <- zz |>
    count(Salinity, name = "Count") |>
    mutate(
      fraction = Count / sum(Count),
      ymax = cumsum(fraction) * 2 * pi,
      ymin = lag(ymax, default = 0)
    )
  
  p <- ggplot(plot_table) +
    geom_arc_bar(aes(x0 = 0, y0 = 0,
                     r0 = 0.00,        # inner radius: smaller = bigger hole, larger = smaller hole
                     r  = 1.00,        # outer radius
                     start = ymin, end = ymax, fill = Salinity),
                 color = "grey20", linewidth = 0.05) +
    scale_fill_manual(values = c("Marine" = "white", "Intermediate" = "grey",
                                 "Brackish" = "black")) +
    coord_fixed() +
    theme_void()
  
  legend <- ggdraw(get_legend(p))
  return (legend)
}

## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## GO:

# name of the dataset used to generate matrix in PCANGSD
pcangsd_dataset <- "full_herr"
# name of the dataset used for DAPC analysis
dapc_dataset <- "wp1_final_bal"

# build complete DF
metadata <- get_metadata(dapc_dataset)
# table with one row per individual
z <- data.frame(
  sample = metadata$sample_id,
  id = metadata$new.id,
  Region = get_fancy_region(metadata$region),
  Period = get_fancy_period(metadata$period),
  Year = metadata$year,
  x = as.numeric(metadata$x),
  y = as.numeric(metadata$y),
  GeneticGroup = grp_assigns(pcangsd_dataset, dapc_dataset, "supplementary_file_7.v2"),
  SpawnTime = grp_assigns(pcangsd_dataset, dapc_dataset, "spring_v_autumn.v2"),
  Salinity = grp_assigns(pcangsd_dataset, dapc_dataset, "baltic_v_atlantic.v2"),
  Inversion = read_ns_column(pcangsd_dataset, dapc_dataset)
)
z <- z |> mutate(Stock = generate_ecotype_column(z))

# loop through periods
periods <- unique(z$Period)

for (period in periods){
  # and places
  regions <- unique(z[which(z$Period == period), ]$Region)
  for (region in regions){
    ## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
    # 1. plot stocks pie:
    pie <- plot_donut(z, period, region)
    ggsave(
      filename = paste0("plots/angsd_matrix/dapc/pie_charts/", gsub("\\s+", "", period),
                        ".", gsub("\\s+", "", region), ".stocks.png"),
      pie, width = 0.1, height = 0.1, dpi = 600, device = png
    )
    
    ## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
    # 2. plot NS inversion frequency
    nsdot <- plot_ns_freq(z, period, region)
    ggsave(
      filename = paste0("plots/angsd_matrix/dapc/pie_charts/", gsub("\\s+", "", period),
                        ".", gsub("\\s+", "", region), ".ns_freq.png"),
      nsdot, width = 0.1, height = 0.1, dpi = 600, device = png
    )
    
    ## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
    # 3. plot salinity pie
    salplot <- plot_salinity_pie(z, period, region)
    ggsave(
      filename = paste0("plots/angsd_matrix/dapc/pie_charts/", gsub("\\s+", "", period),
                        ".", gsub("\\s+", "", region), ".salinity.png"),
      salplot, width = 0.1, height = 0.1, dpi = 600, device = png
    )
    ## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
    # 4. plot N
    npdot <- plot_n(z, period, region)
    ggsave(
      filename = paste0("plots/angsd_matrix/dapc/pie_charts/", gsub("\\s+", "", period),
                        ".", gsub("\\s+", "", region), ".n_circle.png"),
      npdot, width = 0.1, height = 0.1, dpi = 600, device = png
    )
  }
}

## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
# 5. plot legends

ggsave(
  plot = plot_stocks_legend(z),
  filename = paste0("plots/angsd_matrix/dapc/pie_charts/Stocks_legend.pdf"),
  width = 90, height = 90, unit = "mm", dpi = 600, device = pdf
)

ggsave(
  plot = plot_ns_legend(z),
  filename = paste0("plots/angsd_matrix/dapc/pie_charts/NS_legend.pdf"),
  width = 90, height = 90, unit = "mm", dpi = 600, device = pdf
)

ggsave(
  plot = plot_salinity_legend(z),
  filename = paste0("plots/angsd_matrix/dapc/pie_charts/Salinity_legend.pdf"),
  width = 90, height = 90, unit = "mm", dpi = 600, device = pdf
)


## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
# 6. plot maps:
for (period in periods){
  mp <- plot_map(z, period)
  ggsave(
    filename = paste0("plots/angsd_matrix/dapc/pie_charts/", gsub("\\s+", "", period),
                      ".basemap.pdf"),
    mp, width = 90, height = 60, unit = "mm", device = pdf
  )
}


