library(tidyverse)
library(adegenet)

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
      "3" = "Transitional"
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

# assign ecotypes based on dapc results
assign_ecotype <- function(GeneticCluster, SpawnTime, Salinity) {
  
  if (
    GeneticCluster == "North Atlantic" &
    Salinity == "Marine" &
    SpawnTime == "Autumn"
  ) {
    ecotype <- "North Atlantic Autumn-Spawner"
  } else if (
    GeneticCluster == "Norwegian" &
    Salinity == "Marine" &
    SpawnTime == "Spring"
  ) {
    ecotype <- "Norwegian Spring-Spawner"
  } else if (
    GeneticCluster == "North Sea" &
    Salinity == "Marine" &
    SpawnTime == "Autumn"
  ) {
    ecotype <- "North Sea Autumn-Spawner"
  } else if (
    GeneticCluster == "Skagerrak & Kattegat" &
    Salinity == "Transitional" &
    SpawnTime == "Spring"
  ) {
    ecotype <- "Transition Spring-Spawner"
  } else if (
    GeneticCluster == "Baltic" &
    Salinity == "Brackish" &
    SpawnTime == "Spring"
  ) {
    ecotype <- "Baltic Spring-Spawner"
  } else if (
    GeneticCluster == "Britain & Ireland" &
    Salinity == "Marine" &
    SpawnTime == "Autumn"
  ) {
    ecotype <- "Ireland & Britain Autumn-Spawner"
  } else if (
    GeneticCluster == "Baltic" &
    Salinity == "Transitional" &
    SpawnTime == "Autumn"
  ) {
    ecotype <- "Baltic Autumn-Spawner"
  } else if (
    GeneticCluster == "Skagerrak & Kattegat" &
    Salinity == "Marine" &
    SpawnTime == "Spring"
  ) {
    ecotype <- "Transition Spring-Spawner (Marine)"
  } else if (
    GeneticCluster == "Skagerrak & Kattegat" &
    Salinity == "Brackish" &
    SpawnTime == "Spring"
  ) {
    ecotype <- "Transition Spring-Spawner (Brackish)"
  } else if (
    GeneticCluster == "Britain & Ireland" &
    Salinity == "Marine" &
    SpawnTime == "Spring"
  ) {
    ecotype <- "Ireland & Britain Spring-Spawner"
  } else if (
    GeneticCluster == "Skagerrak & Kattegat" &
    Salinity == "Transitional" &
    SpawnTime == "Autumn"
  ) {
    ecotype <- "Baltic Autumn-Spawner"
  } else if (
    GeneticCluster == "North Atlantic" &
    Salinity == "Marine" &
    SpawnTime == "Spring"
  ) {
    ecotype <- "North Atlantic Spring-Spawner"
  } else if (
    GeneticCluster == "North Sea" &
    Salinity == "Marine" &
    SpawnTime == "Spring"
  ) {
    ecotype <- "North Sea Spring-Spawner"
  } else if (
    GeneticCluster == "Skagerrak & Kattegat" &
    Salinity == "Marine" &
    SpawnTime == "Autumn"
  ) {
    ecotype <- "Transition Autumn-Spawner (Marine)"
  } else {
    ecotype <- "Other"
  }
  return (ecotype)
}

generate_ecotype_column <- function(eco_df) {
  mapply(
    function(x1, x2, x3) {
      # Example rule:
      assign_ecotype(x1, x2, x3)
    },
    eco_df$GeneticGroup, eco_df$SpawnTime, eco_df$Salinity,
    USE.NAMES = FALSE
  )
}

save_tables <- function(pcangsd_dataset, dapc_dataset) {
  metadata <- get_metadata(dapc_dataset)
  # table with one row per individual
  z <- data.frame(
    sample = metadata$sample_id,
    id = metadata$new.id,
    Period = metadata$period,
    Year = metadata$year,
    x = metadata$x,
    y = metadata$y,
    GeneticGroup = grp_assigns(pcangsd_dataset, dapc_dataset, "supplementary_file_7.v2"),
    SpawnTime = grp_assigns(pcangsd_dataset, dapc_dataset, "spring_v_autumn.v2"),
    Salinity = grp_assigns(pcangsd_dataset, dapc_dataset, "baltic_v_atlantic.v2")
  )
  # summary table
  fish_type_counts <- z |> 
    group_by(GeneticGroup, Salinity, SpawnTime) |> 
    summarise(
      Count = n(),
      .groups = "drop"
    ) |> 
    arrange(desc(Count))
  # add ecotype assignments
  z <- z |> mutate(Ecotype = generate_ecotype_column(z))
  fish_type_counts <- fish_type_counts |> 
    mutate(Ecotype = generate_ecotype_column(fish_type_counts))
  
  # save tables:
  write.table(
    z, file = paste0("data/angsd_matrix/dapc/ecotypes_table.",
                            pcangsd_dataset, ".", dapc_dataset, ".csv"),
    quote = F, row.names = F, sep = ","
  )
  write.table(
    fish_type_counts,
    file = paste0("data/angsd_matrix/dapc/summary_ecotypes.",
                         pcangsd_dataset, ".", dapc_dataset, ".csv"),
    quote = F, row.names = F, sep = ","
  )
}

## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## CONSTANTS/VARIABLES:

# name of the dataset used to generate matrix in PCANGSD
pcangsd_dataset <- "full_herr"
# name of the dataset used for DAPC analysis
dapc_dataset <- "wp1_final_bal"

## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## Assign Ecotypes and write table
save_tables(pcangsd_dataset, dapc_dataset)


################################################################################

library(dplyr)
metadata <- get_metadata(dapc_dataset)
# table with one row per individual
z <- data.frame(
  sample = metadata$sample_id,
  id = metadata$new.id,
  Region = metadata$region,
  Period = get_fancy_period(metadata$period),
  Year = metadata$year,
  x = metadata$x,
  y = metadata$y,
  GeneticGroup = grp_assigns(pcangsd_dataset, dapc_dataset, "supplementary_file_7.v2"),
  SpawnTime = grp_assigns(pcangsd_dataset, dapc_dataset, "spring_v_autumn.v2"),
  Salinity = grp_assigns(pcangsd_dataset, dapc_dataset, "baltic_v_atlantic.v2")
)
# summary table
fish_type_counts <- z |> 
  group_by(GeneticGroup, Salinity, SpawnTime) |> 
  summarise(
    Count = n(),
    .groups = "drop"
  ) |> 
  arrange(desc(Count))
# add ecotype assignments
z <- z |> mutate(Ecotype = generate_ecotype_column(z))
fish_type_counts <- fish_type_counts |> 
  mutate(Ecotype = generate_ecotype_column(fish_type_counts))

region <- "Skagerrak_&_Kattegat"
region <- "Norway"

salinity_freq <- z %>% filter(Region == region) |> 
  group_by(Period, Salinity) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Period) %>%
  mutate(prop = n / sum(n))

salinity_freq$Period <- factor(
  salinity_freq$Period,
  levels = c(
    "Before 1700s",
    "1747–1805 Sillperiod",
    "Between Sillperiods",
    "1877–1906 Sillperiod",
    "Present-day"
  ),
  ordered = TRUE
)

ggplot(
  salinity_freq,
  aes(x = Period,
      y = prop,
      color = Salinity,
      group = Salinity)
) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    x = NULL,
    y = "Proportion of observations"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(
  salinity_freq,
  aes(x = Period,
      y = prop,
      fill = Salinity)
) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = NULL,
    y = "Proportion of observations"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

