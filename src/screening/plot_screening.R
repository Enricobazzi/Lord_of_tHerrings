library(tidyverse)

setwd("~/Documents/Lord_of_tHerrings/")

screen_df <- read.table("data/summary_screening/summary_screening.nocommanospace.tsv",
                        sep = "\t", header = TRUE)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### BOXPLOTS ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

box_museum <- ggplot() +
  geom_boxplot(data = screen_df, outlier.shape = NA,
               aes(x = museum, y = genome_recovery, fill = century)) 
box_museum <- box_museum +
  geom_jitter(data = screen_df, shape = 21,
              aes(x = museum, y = genome_recovery, fill = century),
              alpha = 1, width = 0.1, size = 3)

box_century <- ggplot() +
  geom_boxplot(data = screen_df, outlier.shape = NA,
               aes(x = century, y = genome_recovery, fill = region)) 
box_century <- box_century +
  geom_jitter(data = screen_df, shape = 21,
              aes(x = century, y = genome_recovery, fill = region),
              alpha = 1, width = 0.2, size = 3)

box_region <- ggplot() +
  geom_boxplot(data = screen_df, outlier.shape = NA,
               aes(x = region, y = genome_recovery, fill = century)) 

box_region <- box_region +
  geom_jitter(data = screen_df, shape = 21,
              aes(x = region, y = genome_recovery, fill = century),
              alpha = 1, width = 0.2, size = 3)

ggsave("plots/screening/box_museum.png",
       plot = box_museum,
       width = 12, height = 5, dpi = 300)
ggsave("plots/screening/box_century.png",
       plot = box_century,
       width = 12, height = 5, dpi = 300)
ggsave("plots/screening/box_region.png",
       plot = box_region,
       width = 12, height = 5, dpi = 300)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### CALCULATE LANES ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# yield cutoffs 
screen_gr15 <- screen_df %>% filter(genome_recovery > 0.15)
nlanes_gr15 <- sum(screen_gr15$prop_lane)

screen_gr20 <- screen_df %>% filter(genome_recovery > 0.20)
nlanes_gr20 <- sum(screen_gr20$prop_lane)

# calculate yield + location
screen_gr10_loc10 <- screen_df %>%
  filter(genome_recovery > 0.10) %>%
  group_by(location) %>%
  filter(n() >= 4) %>%
  slice_max(order_by = genome_recovery, n = 12, with_ties = FALSE) %>%
  ungroup()
nlanes_gr10_loc10 <- sum(screen_gr10_loc10$prop_lane)

# use 0.2 cutoff and add 9 samples from Högsten <- starting point
final_sp <- rbind(
  data.frame(
    screen_df %>% filter(genome_recovery > 0.20)
  ),
  data.frame(
    screen_df %>% filter(location == "Högsten" & genome_recovery > 0.1)
  )
)
sum(final_sp$prop_lane)

# remove samples to reach <6 lanes:
# 3 worst Hogsten samples: 301, 306 310
# 3 worst Stavanger samples: 179, 181, 185
# 1 worst Kämpinge/Höllviken sample: 406
final <- final_sp %>% filter(
  str_detect(ID, paste(c(
    "ND301", "ND310", "ND179", "ND181", "ND185", "ND406"
  ), collapse="|"), negate = TRUE)
)
nlanes_final <- sum(final$prop_lane)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### MAP SAMPLES ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

library(rgdal)
library(raster)
library(plotly)
library(htmlwidgets)
library(htmltools)
library(leaflet)
library(sf)
library(ggrepel)

century_color_dict <- c(
  "6th" = "#CC3333",
  "12th" = "#FF6600",
  "16th" = "#009900",
  "17th" = "#006666",
  "18th" = "#0066FF",
  "19th" = "#9933FF"
)

get_n_samples_from_loc <- function(loc, df){
  return(nrow(df %>% filter(location == loc)))
}

get_plot_df <- function(df){
  plot_df <- data.frame(
    x = df$x,
    y = df$y,
    location = df$location,
    color = as.character(century_color_dict[df$century]),
    n = as.numeric(lapply(df$location, get_n_samples_from_loc, df=df)),
    grr = df$genome_recovery,
    century = df$century,
    museum = df$museum
  ) %>%
    group_by(location) %>% summarize(avg_grr = mean(grr), across()) %>%
    select(-grr) %>% distinct() %>% 
    group_by(across(-museum)) %>% 
    summarise(museum = paste(unique(museum), collapse = "_"), .groups = "drop") %>%
    mutate(location = str_split_fixed(location, "__", 2)[,1]) %>%
    dplyr::ungroup()
  return(as.data.frame(plot_df))
}

plot_map <- function(){
  # Load shapefile as sp object
  world_map_sp <- rgdal::readOGR("~/Documents/Lynxtrogression_v2/plots/ms_figures/ne_50m_land/ne_50m_land.shp")
  # Convert to sf object
  world_map_sf <- st_as_sf(world_map_sp)
  # Plot
  wm <- ggplot() +
    geom_sf(data = world_map_sf, fill = "grey80", color = "grey0") +
    coord_sf(
      xlim = c(-8, 30),
      ylim = c(50, 67),
      expand = FALSE,
      default_crs = st_crs(world_map_sf)
    )
  return(wm)
}

plot_coords <- function(plot_df, title){
  cmap <- plot_map() + 
    geom_point(data = plot_df, aes(x = x, y = y, fill = century), 
               shape = 21, size = 3) +
    scale_color_manual(plot_df$color) +
    geom_label_repel(data = plot_df, aes(x = x, y = y,
                                         label = paste0(location, ": ", n)),
                     min.segment.length = unit(0.001, "mm"),
                     force_pull = 0.5, max.overlaps = 50, size = 3) +
    xlab("") + ylab("") +
    ggtitle(title)
  return(cmap)
}

leaf_coords <- function(plot_df){
  leaf_df <- leaflet(plot_df) %>%
    addTiles() %>%
    addCircleMarkers(
      ~x, ~y,
      radius = 3,
      popup = ~paste0(
        "<div style='font-size:18px;'>",
        "<b>", location, "</b><br>",
        "<i>", museum, "</i><br>",
        "n = ", n, "<br>",
        "average yield = ", round(avg_grr, digits = 5), "<br>",
        "age: ", century, " century"
      ),
      popupOptions = popupOptions(maxWidth = 300))
}

plot_gr15 <- get_plot_df(screen_gr15)
map_gr15 <- plot_coords(plot_gr15, paste0("Samples with 15% recovery\n", round(nlanes_gr15, digits = 3), " lanes"))
leaf_gr15 <- leaf_coords(plot_gr15)
ggsave("plots/screening/map_gr15.pdf", plot = map_gr15, height = 5, width = 5)
saveWidget(leaf_gr15, "plots/screening/map_gr15.html", selfcontained = TRUE)


plot_gr20 <- get_plot_df(screen_gr20)
map_gr20 <- plot_coords(plot_gr20, paste0("Samples with 20% recovery\n", round(nlanes_gr20, digits = 3), " lanes"))
leaf_gr20 <- leaf_coords(plot_gr20)
ggsave("plots/screening/map_gr20.pdf", plot = map_gr20, height = 5, width = 5)
saveWidget(leaf_gr20, "plots/screening/map_gr20.html", selfcontained = TRUE)

plot_gr10_loc10 <- get_plot_df(screen_gr10_loc10)
map_gr10_loc10 <- plot_coords(plot_gr10_loc10, paste0("Samples with 10% recovery\nminimum 4 maximum 12 per location\n", round(nlanes_gr10_loc10, digits = 3), " lanes"))
leaf_gr10_loc10 <- leaf_coords(plot_gr10_loc10)
ggsave("plots/screening/map_gr10_loc10.pdf", plot = map_gr10_loc10, height = 5, width = 5)
saveWidget(leaf_gr10_loc10, "plots/screening/map_gr10_loc10.html", selfcontained = TRUE)

plot_final <- get_plot_df(final)

map_final <- plot_map() + 
  geom_point(data = read.table("data/summary_screening/already_sequenced_data.txt", header = T),
             aes(x = x, y = y, fill = century), 
             shape = 23, size = 4) +
  geom_point(data = plot_final, aes(x = x, y = y, fill = century), 
             shape = 21, size = 3.5) +
  geom_label_repel(data = rbind(plot_final,
                                read.table("data/summary_screening/already_sequenced_data.txt", header = T)),
    aes(x = x, y = y, label = paste0(location, ": ", n)),
                   min.segment.length = unit(0.001, "mm"),
                   force_pull = 0.5, max.overlaps = 50, size = 2.5) +
  xlab("") + ylab("") +
  ggtitle(paste0("103 (124) samples - 17 (22) locations\n", round(nlanes_final, digits = 3), " lanes"))
ggsave("plots/screening/map_final.pdf", plot = map_final, height = 6, width = 8)
leaf_final <- leaf_coords(plot_final)
saveWidget(leaf_final, "plots/screening/map_final.html", selfcontained = TRUE)

