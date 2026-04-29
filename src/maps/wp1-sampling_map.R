library(tidyverse)
library(sf)
library(ggrepel)

get_samples_from_file <- function(file_path) {
  samples <- read.table(file_path)[, 1] |> as.character()
  return(samples)
}

plot_map <- function(world_map_sp, xl = c(-80, 30), yl = c(40, 70)){
  # Load shapefile as sp object
  # Convert to sf object
  world_map_sf <- st_as_sf(world_map_sp)
  # Plot
  wm <- ggplot() +
    geom_sf(data = world_map_sf, fill = "grey80", color = "grey0") +
    coord_sf(
      xlim = xl,
      ylim = yl,
      expand = FALSE,
      default_crs = st_crs(world_map_sf)
    )
  return(wm)
}

# define samples
subset_samples_name <- "wp1_final"
subset_samples_file <- paste0("data/angsd_matrix/bamlists/", subset_samples_name, ".sample_list.txt")
samples <- get_samples_from_file(subset_samples_file)
samples <- samples[samples != "HER135"] # Exclude sample HER135 - should be fixed in the future

# read samples metadata
tab <- read.table("data/samples_table.csv", sep = ",", header = TRUE, na.strings = c("UNKNOWN", "NA"))
# tab <- read.table("~/Documents/Silly-periods/data/samples_table.csv", sep = ",", header = TRUE, na.strings = c("UNKNOWN", "NA"))
tab <- tab %>%
  filter(sample_id %in% samples)

# read in world map
world_map_sp <- sf::read_sf("data/maps/ne_50m_land/ne_50m_land.shp")

modern <- tab %>%
  filter(
    time == "modern"
  ) %>%
  mutate(
    x = as.numeric(x),
    y = as.numeric(y)
  )

modern_summary <- modern |>
  count(x, y, spawn)

modern_map <- plot_map(world_map_sp, xl = c(-25, 20), yl = c(50, 70)) +
  geom_point(
    data = modern_summary,
    aes(x = x, y = y, fill = spawn), 
    size = 3, shape = 21
    ) +
  scale_fill_manual(values = c("autumn" = "darkcyan", "spring" = "gold2", "winter" = "darkblue")) +
  geom_label_repel(
    data = modern_summary |> filter(n > 1),
    aes(x = x, y = y, label = n),
    force = 1,
    size = 4,
    box.padding = 0.3,
    label.padding = 0.1,
    label.size = 0.2,       # border thickness
    fill = "white",         # background color
    color = "black"
  ) +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey20")
  ) +
  xlab("") + ylab("")
modern_map
ggsave("plots/maps/wp1-modern_map.png", modern_map, width = 7, height = 5)
# writeLines(modern$sample_id, "data/rdaforest/modern_samples.txt")

historical <- tab |>
  filter(
    time == "historical" | time == "sillperioder"
  ) |>
  mutate(
    x = as.numeric(x),
    y = as.numeric(y)
  ) |> 
  mutate(
    century = case_when(
      #is.na(year) ~ NA_character_,
      is.na(year) ~ "1800s",
      year < 1500 ~ "before 1500",
      year < 1800 ~ "1500 to 1800",
      TRUE ~ "1800s"
      )
    )
    
# paste0(floor(year / 100) + 1, "00s")
historical_summary <- historical %>%
  count(x, y, century)

historical_summary <- historical |>
  group_by(x, y, century, time) |>
  summarise(
    n = n(),
    mean_wg_depth = mean(wg.depth, na.rm = TRUE),
    min_wg_depth  = min(wg.depth, na.rm = TRUE),
    max_wg_depth  = max(wg.depth, na.rm = TRUE),
    .groups = "drop"
  )

historical_map <- plot_map(world_map_sp, xl = c(-10, 20), yl = c(50, 65)) +
  geom_point(
    data = historical_summary,
    aes(x = x, y = y, color = century, fill = century, shape = time),
    size = 3,
    # position = position_jitterdodge(jitter.width = 2, dodge.width = 2, seed = 1)
    ) +
  geom_point(
   data = historical_summary,
   aes(x = x, y = y, shape = time),
   size = 3, color = "black",
   # position = position_jitterdodge(jitter.width = 2, dodge.width = 2, seed = 1)
  ) +
  scale_shape_manual(values = c(21, 24)) +
  # scale_fill_binned(name = "year", n.breaks = 5) +
  geom_label_repel(
    data = historical_summary,
    aes(x = x, y = y,
        label = paste0(n, ": ~",
          round(mean_wg_depth, 1), "X")
        ),
    force = 20,             # any value > 0 works
    segment.color = "black",
    segment.size = 0.3,
    segment.alpha = 1,
    size = 3, max.overlaps = 500,
    box.padding = 0.1,
    label.padding = 0.1,
    label.size = 0.2,       # border thickness
    fill = "white",         # background color
    color = "black"         # text color
  ) + 
  xlab("") + ylab("") +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey20")
  ) 
historical_map
ggsave("plots/maps/wp1-historical_map.png", historical_map, width = 7, height = 5)

# Ggnewscale
# pca.plot <- autoplot(pc_ind,
#                      data = ind.level.df,
#                      loadings = TRUE,
#                      loadings.level = TRUE,
#                      loadings.label = TRUE, 
#                      loadings.label.size = 3, 
#                      loadings.colour = "grey20",
#                      loadings.label.colour = "grey20",
#                      color ="white") +
#   stat_ellipse(type = "norm", color = "grey20", level = 0.95, linewidth = 0.6) + 
#   geom_point(aes(color = plant_sp), alpha = 0.75) +
#   scale_color_manual(values = mycols) + 
#   new_scale_color()+
#   stat_ellipse(aes(color = continent), type = "norm", linetype = "dashed", level = 0.95) +
#   scale_color_manual(values = col.cont) +
#   theme(legend.position = "none")#

historical <- tab |>
  filter(
    species == "clupea_harengus",
    (region != "E-PACIFIC" & region != "NW-ATL" & region != "BALTIC" & region != "EUROPE") | is.na(region),
    year < 1930 | is.na(year),
    study == "current",
    wg.depth > 0
  ) |>
  mutate(
    x = as.numeric(x),
    y = as.numeric(y)
  ) |> 
  mutate(
    century = case_when(
      is.na(year) ~ NA_character_,
      year < 1500 ~ "before 1500",
      year < 1800 ~ "1500 to 1800",
      TRUE ~ "1800s"
    )
  )
modern |> filter(region == "TRANS") |> count(location)
modern |> filter(str_detect(location, "Røvær"))

modern |> 
  filter(region == "TRANS")

historical |> 
  filter(time == "historical" & region == "TRANS")

historical |> 
  filter(time == "sillperioder" & region == "TRANS" & (location == "Gullholmen" | location == "Masthugget"))

historical |> 
  filter(time == "historical" & region == "TRANS")

data.frame(dapc_df |>
  filter(time == "modern" & loc == "Norwegian Sea & Fjords" & group != 5))$old_id

dapc_df[startsWith(rownames(dapc_df), "fold"), ]
