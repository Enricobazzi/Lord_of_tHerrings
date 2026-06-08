library(tidyverse)
library(ggrepel)
library(sf)
library(cowplot)
library(scatterpie)

# palette of colors
colpal <- c(
  "Norwegian Spring-Spawner" = "#02d97c",
  "North Atlantic Autumn-Spawner" = "#FED439FF",
  "North Atlantic Spring-Spawner" = "#FFF5C8",
  "North Sea Autumn-Spawner" = "#ED3911",
  "North Sea Spring-Spawner" = "#FF7518",
  "Ireland & Britain Spring-Spawner" = "#91331FFF",
  "Ireland & Britain Autumn-Spawner" = "#a6721e",
  "Baltic Spring-Spawner" = "#71D0F5FF",
  "Baltic Autumn-Spawner" = "#1b639e",
  "Transition Spring-Spawner" = "#A55BA4",
  "Transition Spring-Spawner (Marine)" = "#733372",
  "Transition Spring-Spawner (Brackish)" = "#E0C7E0",
  "Transition Autumn-Spawner (Marine)" = "#C58AC4"
)

# colpal <- c(
#   "Norwegian Spring-Spawner" = "#02d97c",
#   "North Atlantic Autumn-Spawner" = "#FED439FF",
#   "North Sea Autumn-Spawner" = "#ED3911",
#   "North Atlantic Spring-Spawner" = "#3f704d",
#   "Britain & Ireland Spring-Spawner" = "#FF7518",
#   "Ireland & Britain Autumn-Spawner" = "#a6721e",
#   "Transition Spring-Spawner (Marine)" = "#FB9A99",
#   "Baltic Spring-Spawner" = "#1b639e",
#   "Transition Spring-Spawner" = "#b370b2",
#   "Baltic Autumn-Spawner" = "#440154FF",
#   "Transition Spring-Spawner (Brackish)" = "#71D0F5FF",
#   "Transition Autumn-Spawner (Marine)" = "#a6721e",#91331FFF
#   "North Sea Spring-Spawner" = "#FF7518"
# )

# name of the dataset used to generate matrix in PCANGSD
pcangsd_dataset <- "full_herr"
# name of the dataset used for DAPC analysis
dapc_dataset <- "wp1_final_bal"

eco_df <- read.table(
  file = paste0("data/angsd_matrix/dapc/ecotypes_table.",
                pcangsd_dataset, ".", dapc_dataset, ".csv"),
  header = T, sep = ","
)

summary_table <- eco_df |>
  count(Period, x, y, Ecotype, name = "n") |>
  pivot_wider(
    names_from = Ecotype,
    values_from = n,
    values_fill = 0
  )

pie_data <- summary_table |>
  pivot_longer(
    cols = -c(Period, x, y),
    names_to = "Ecotype",
    values_to = "Count"
  ) |> 
  group_by(Period, x, y) |> 
  group_walk(~ {

    p <- ggplot(.x, aes(x = 1, y = Count, fill = Ecotype)) +
      geom_col(width = 0.1, color = "grey20", linewidth = 0.5) +
      coord_polar(theta = "y", start = 0) +
      scale_fill_manual(values = c(colpal)) +
      theme_void() + theme(legend.position="none")
    
    ggsave(
      paste0(
        "pie_charts/Period_", .y$Period,
        "_X_", .y$x,
        "_Y_", .y$y,
        ".png"
      ),
      p,
      width = 0.5,
      height = 0.5,
      dpi = 300, device = png
    )
  }
  )

# load world map shapefile
world_map_sf <- sf::read_sf("data/maps/ne_50m_land/ne_50m_land.shp")
world_map_sf <- st_transform(world_map_sf, crs = 4326)
# set limits
xl <- c(-20, 25)
yl <- c(50, 70)

# Sea  : "#CFE8F3"
# Land : "#EAF4E1"
# Sea  : "#A9D6E5"
# Land : "#F1E9DA"
# Sea  : "#1F2937"
# Land : "#374151"
# Sea  : "#DCEEF2"
# Land : "#DCE8D0"

# Sea       "#DCEEF2"
# Land      "#F3EED9"
# Borders   "#A8A8A8"
# Coastline "#808080"

# starting map
basemap <- ggplot() +
  geom_sf(data = world_map_sf, fill = "#E3E3E3", color = "grey0", linewidth = 0.1) +
  xlab("") +
  ylab("") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "#F7F7F7"),
    panel.border = element_rect(color = "#A8A8A8", fill = NA, linewidth = 0.5)
  )

basemap +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    plot.margin = margin(30, 5, 0, 0)
  )

tp <- unique(eco_df$Period)[1] # iterate through time periods
for (tp in unique(eco_df$Period)){
  st <- summary_table |> filter(Period == tp) # subset table for period
  
  # new basemap for the period
  basemap2 <- basemap +
    geom_point(
      data = st, aes(x = x, y = y), shape = 1, size = 0.5
    ) +
    coord_sf(
      xlim = xl, ylim = yl, expand = FALSE, default_crs = NULL
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      plot.margin = margin(10, 5, 0, 0),
      axis.title = element_text(size = 6, face = "bold"),
      axis.text = element_text(size = 6),
      axis.ticks = element_line(linewidth = 0.1)
    )
  
  # save pdf of basemap for later manipulation :)
  ggsave(
    plot = basemap2,
    filename = paste0("plots/angsd_matrix/dapc/", tp, ".basemap.pdf"),
    width = 90, height = 60, unit = "mm", device = pdf
  )
    
  for (i in seq_len(nrow(st))) {
    x_val <- st$x[i]
    y_val <- st$y[i]
    # add pie charts
    filename <- paste0(
      "pie_charts/Period_", tp,
      "_X_", x_val,
      "_Y_", y_val,
      ".png"
    )
    basemap2 <- basemap2 +
      draw_image(filename, x = x_val, y = y_val, vjust = 0.5, hjust = 0.5, scale = 3) 
  }
  
  ggsave(
    plot = basemap2,
    filename = paste0("plots/angsd_matrix/dapc/", tp, ".piemap.png"),
    width = 90, height = 90, unit = "mm", dpi = 600, device = png
  )
}

# legend
df <- data.frame(
  x = seq_along(colpal),
  y = 1,
  ecotype = names(colpal)
)

p <- ggplot(df, aes(x, y, fill = ecotype)) +
  geom_point(shape = 21, size = 5) +
  scale_fill_manual(values = colpal) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

legend <- ggdraw(get_legend(p))
ggsave(
  plot = legend,
  filename = paste0("plots/angsd_matrix/dapc/Ecotypes_legend.png"),
  width = 90, height = 90, unit = "mm", dpi = 600, device = png
)


# time
library(ggplot2)

present <- as.numeric(format(Sys.Date(), "%Y"))

ggplot() +
  
  # faded pre-1600 continuation
  annotate(
    "rect",
    xmin = 1565, xmax = 1600,
    ymin = -Inf, ymax = Inf,
    fill = "grey85",
    alpha = 0.35
  ) +
  
  # full-height timeline
  annotate(
    "rect",
    xmin = 1600, xmax = present,
    ymin = -Inf, ymax = Inf,
    fill = "grey92"
  ) +
  
  # highlighted periods
  annotate(
    "rect",
    xmin = 1747, xmax = 1809,
    ymin = -Inf, ymax = Inf,
    fill = "#D8A24A",
    alpha = 0.85
  ) +
  annotate(
    "rect",
    xmin = 1877, xmax = 1906,
    ymin = -Inf, ymax = Inf,
    fill = "#D8A24A",
    alpha = 0.85
  ) +
  
  scale_x_continuous(
    limits = c(1560, present),
    breaks = seq(1600, present, 50),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  
  theme_void()  +
  theme(
    axis.text.x = element_text(size = 11),
    axis.ticks.x = element_line(),
    plot.margin = margin(0,0,0,0)
  )
ggsave("timeline.png", width = 8, height = 0.6)
