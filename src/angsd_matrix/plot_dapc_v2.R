library(tidyverse)
library(adegenet)
library(ggrepel)

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

# variance of LD1 & LD2
get_variance_label <- function(dapc_obj, axn) {
  return(paste0("LD", axn, " (", round(dapc_obj$eig[axn] / sum(dapc_obj$eig) * 100, 2), "%)"))
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

# plot dapc:
plot_dapc <- function(dapc_obj, metadata, color_by = "region") {
  
  dapc_df <- data.frame(
    id = metadata$new.id,
    Region = unlist(lapply(metadata$region, get_fancy_region)),
    Spawn = metadata$spawn,
    LD1 = dapc_obj$ind.coord[, 1],
    LD2 = dapc_obj$ind.coord[, 2],
    group = droplevels(dapc_obj$assign),
    old_id = metadata$sample_id,
    Period = unlist(lapply(metadata$period, get_fancy_period))
  )
  
  # dataframe for centroids with names:
  centroids <- dapc_df %>%
    group_by(group) %>%
    summarise(
      mean_LD1 = median(LD1, na.rm = TRUE),
      mean_LD2 = median(LD2, na.rm = TRUE),
      group_name = first(group)
    )
  
  if (color_by == "region") {
    dapc_plot <- ggplot (data = dapc_df,
                         aes(x = LD1, y = LD2, colour = Region, fill = Region, shape = Period))
  } else if (color_by == "spawn") {
    dapc_plot <- ggplot (data = dapc_df,
                         aes(x = LD1, y = LD2, colour = Spawn, fill = Spawn, shape = Period))
  }
  
  dapc_plot <- dapc_plot +
    # lines to group centroid
    geom_segment(
      data = dapc_df %>% left_join(centroids, by = "group"),
      aes(x = LD1, y = LD2, xend = mean_LD1, yend = mean_LD2),
      linewidth = 0.1,
      alpha = 0.8
    ) +
    # ellipse
    # stat_ellipse(aes(group = group), level = 0.5, linetype = 1, linewidth = 0.3) +
    # points
    geom_point(alpha = 0.7, size = 2) +
    # geom_point(
    #   data = dapc_df |> filter(Period == "Present-day"),
    #   aes(x = LD1, y = LD2), inherit.aes = F, shape = 1, size = 2
    # ) +
    geom_label_repel(data = centroids, aes(x = mean_LD1, y = mean_LD2, label = group),
                     fill = "white", size = 2, inherit.aes = F,
                     label.padding = unit(0.1, "lines"), segment.size = 0.05, force = 0.01) +
    # fancy plot:
    xlab(get_variance_label(dapc_obj, 1)) + ylab(get_variance_label(dapc_obj, 2)) +
    scale_fill_manual(values = c(colpal)) +
    scale_color_manual(values = c(colpal)) +
    scale_shape_manual(values = time_shapes) +
    theme_minimal(base_size = 7) # + theme(legend.position = "none")
  return(dapc_plot)
}


## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## CONSTANTS/VARIABLES:

# name of the dataset used to generate matrix in PCANGSD
pcangsd_dataset <- "full_herr"
# name of the dataset used for DAPC analysis
dapc_dataset <- "wp1_final_bal"
# sites
sites_name_lst <- c(
  "supplementary_file_7.v2",
  "baltic_v_atlantic.v2",
  "spring_v_autumn.v2"
)
# palette of colors
colpal <- c(
  "Bothnia" = "#440154FF",
  "Baltic" = "#1b639e",
  "Belt" = "#71D0F5FF",
  "Skagerrak & Kattegat" = "#b370b2",
  "North Sea" = "#ED3911",
  "Britain & Ireland" = "#91331FFF",
  "Norway" = "#02d97c",
  "North Atlantic" = "#FED439FF",
  "autumn" = "#91331FFF",
  "autumn/winter" = "#440154FF",
  "spring" = "#02d97c",
  "summer" = "#FED439FF",
  "winter" = "#1b639e"
)
# palette of shapes
time_shapes <- c(
  "Before 1700s" = 22,
  "1747–1805 Sillperiod" = 24,
  "Between Sillperiods" = 23,
  "1877–1906 Sillperiod" = 25,
  "Present-day" = 21
)


## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## PLOT DAPC RESULTS:

for (sites_name in sites_name_lst) {
  best_k <- read_best_k(dapc_dataset, sites_name)
  dapc_obj <- load_dapc(pcangsd_dataset, dapc_dataset, sites_name, best_k)
  metadata <- get_metadata(dapc_dataset)
  
  if (sites_name == "spring_v_autumn.v2") {
    plt <- plot_dapc(dapc_obj, metadata, color_by = "spawn")
  } else {
    plt <- plot_dapc(dapc_obj, metadata)
  }
  
  ggsave(
    filename = paste0("plots/angsd_matrix/dapc/", dapc_dataset, ".", sites_name, ".k", best_k, ".dapc.png"),
    plot = plt, width = 180, height = 90, unit = "mm", dpi = 300
  )
}

