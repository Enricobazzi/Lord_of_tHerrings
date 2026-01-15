library(tidyverse)
library(ggplot2)
library(ggrepel)
library(adegenet)
library(cowplot)

get_samples_from_file <- function(file_path) {
  samples <- read.table(file_path)[, 1] |> as.character()
  return(samples)
}

get_samples_subset <- function(all_samples_file, subset_samples_file) {
  all_samples <- get_samples_from_file(all_samples_file)
  subset_samples <- get_samples_from_file(subset_samples_file)
  return(all_samples[all_samples %in% subset_samples])
}

read_matrix_from_file <- function(file_path) {
  mat <- as.matrix(read.table(file_path))
  return(mat)
}

get_samples_matrix <- function(matrix_file, samples, all_samples) {
  mat <- read_matrix_from_file(matrix_file)
  if (nrow(mat) == length(samples) && ncol(mat) == length(samples)) {
    rownames(mat) <- samples
    colnames(mat) <- samples
    return(mat)
  } else {
    sample_indices <- which(all_samples %in% samples)
    filtered_mat <- mat[sample_indices, sample_indices]
    rownames(filtered_mat) <- all_samples[sample_indices]
    colnames(filtered_mat) <- all_samples[sample_indices]
    return(filtered_mat)
  }
}

get_sample_metadata <- function(sample_data_file, samples) {
  sample_data <- read.table(sample_data_file, sep = ",",
                            header = TRUE, na.strings = "UNKNOWN")
  sample_data <- sample_data[sample_data$sample_id %in% samples, ]
  return(sample_data)
}

get_sample_ids <- function(sample_metadata, samples) {
  sample_ids <- sample_metadata$new.id[match(samples, sample_metadata$sample_id)]
  return(sample_ids)
}

get_locs <- function(sample_ids) {
  loc <- ifelse(
    stringr::str_split(sample_ids, "_", simplify = TRUE)[, 1] %in% c("iceland", 
                                                                     "faroe"),
    "Iceland & Faroe Islands",
    ifelse(
      stringr::str_split(sample_ids, "_", simplify = TRUE)[, 1] %in% c("more", "norwegian",
                                                                       "haugesund", "karmoy",
                                                                       "rovaer", "bergen",
                                                                       "stavanger", "foldfjorden",
                                                                       "norway", "unknown"),
      "Norwegian Sea & Fjords",
      ifelse(
        stringr::str_split(sample_ids, "_", simplify = TRUE)[, 1] %in% c("idefjord", "maseskar",
                                                                         "kaupang", "risor",
                                                                         "koster",
                                                                         "dynekilen",
                                                                         "nyalodose", "masthugget",
                                                                         "gullholmen",
                                                                         "kalvsund"),
        "Skagerrak-Kattegat",
        ifelse(
          stringr::str_split(sample_ids, "_", simplify = TRUE)[, 1] %in% c("northsea",
                                                                           "netherlands"),
          "North Sea",
          ifelse(
            stringr::str_split(sample_ids, "_", simplify = TRUE)[, 1] %in% c("celtic",
                                                                             "downs", "isleofman", "york",
                                                                             "scotland", "lyminge"),
            "United Kingdom",
            "Other"
          )
        )
      )
    )
  )
  return(loc)
}

get_times <- function(sample_metadata, samples) {
  times <- sample_metadata$time[match(samples, sample_metadata$sample_id)]
  return(times)
}

assign_group_name <- function(df) {
  group_names <- df %>%
    group_by(group) %>%
    mutate(group_name = names(sort(table(loc), decreasing = TRUE))[1]) %>%
    ungroup() %>%
    pull(group_name)
  return(group_names)
}


all_samples_name <- "wp1_all"
subset_samples_name <- "wp1_final"
sites_name <- "sf7_sites"
all_samples_file <- paste0("data/angsd_matrix/bamlists/", all_samples_name, ".sample_list.txt")
subset_samples_file <- paste0("data/angsd_matrix/bamlists/", subset_samples_name, ".sample_list.txt")
# subset_samples_file <- paste0("data/angsd_matrix/bamlists/", all_samples_name, ".sample_list.txt")
matrix_file <- paste0("data/angsd_matrix/", all_samples_name, ".", sites_name, ".pcangsd.cov")
sample_data_file <- "data/samples_table.csv"
all_samples <- get_samples_from_file(all_samples_file)
all_samples <- all_samples[all_samples != "HER135"] # Exclude sample HER135 - should be fixed in the future
samples <- get_samples_subset(all_samples_file, subset_samples_file)
samples <- samples[samples != "HER135"] # Exclude sample HER135 - should be fixed in the future
data_table <- get_sample_metadata(sample_data_file, samples)
sample_ids <- get_sample_ids(data_table, samples)
matrix <- get_samples_matrix(matrix_file, samples, all_samples)

# prior on first 4 (>90% of variance explained)
grp <- find.clusters(matrix, n.pca = 4, n.clust = 5)
grp
candidate_grps <- grp$grp
xval <- xvalDapc(matrix, candidate_grps, n.pca.max = 100,
                 result = "groupMean",
                 n.pca = NULL, n.rep = 50, xval.plot = FALSE)

# n.pca and n.da chosen based on xval results (10 and 4)
DAPC <- dapc(matrix, candidate_grps,
             n.pca = xval$DAPC$n.pca,
             n.da = xval$DAPC$n.da)

## if loading from saved DAPC object:
# DAPC <- readRDS("data/angsd_matrix/DAPC_wp1_final_sf7_sites.rds")
# scatter(DAPC)

# biplot of dapc - prepare data frame:
dapc_df <- data.frame(DAPC$ind.coord)
rownames(dapc_df) <- sample_ids
dapc_df$loc <- get_locs(sample_ids)
dapc_df$time <- get_times(data_table, samples)
dapc_df$group <- DAPC$assign
dapc_df$group_names <- assign_group_name(dapc_df)
# dataframe for centroids with names:
group_means <- dapc_df %>%
  group_by(group) %>%
  summarise(
    mean_LD1 = mean(LD1, na.rm = TRUE),
    mean_LD2 = mean(LD2, na.rm = TRUE),
    group_name = first(group_names)
  )
da_percent_variance <- round(100 * DAPC$eig / sum(DAPC$eig), 1)

# combine dfs for label repel
repel_df <- bind_rows(
  dapc_df %>%
    mutate(
      label = NA_character_
    ),
  group_means %>%
    transmute(
      LD1 = mean_LD1,
      LD2 = mean_LD2,
      label = group_name
    )
)
# new df for label positions
label_positions <- tibble(
  group = group_means$group,
  label_x = c(  -4, 2.5, 4.5,    5, -3.5),
  label_y = c(-0.5,  -6, 4.5, -4.5,  4.5)
)
group_means_labeled <- group_means %>%
  left_join(label_positions, by = "group")

# plot
p1 <- ggplot(dapc_df, aes(x = LD1, y = LD2,
                    fill = loc, color = time, shape = time,
                    size = time, stroke = time)) +
  ## lines to group centroid
  geom_segment(
    data = dapc_df %>% left_join(group_means, by = "group"),
    aes(x = LD1, y = LD2, xend = mean_LD1, yend = mean_LD2),
    inherit.aes = FALSE,
    color = "grey70",
    linewidth = 0.3,
    alpha = 0.5
  ) +
  ## lines to the group labels
  geom_segment(
    data = group_means_labeled,
    aes(
      x = mean_LD1,
      y = mean_LD2,
      xend = label_x,
      yend = label_y
    ),
    inherit.aes = FALSE,
    color = "grey70",
    linewidth = 0.3,
    alpha = 0.5
  ) + 
  ## centroid points
  geom_point(
    data = group_means,
    aes(x = mean_LD1, y = mean_LD2),
    inherit.aes = FALSE,
    shape = 1,
    size = 3,
    stroke = 0.2,
    color = "grey70"
  ) +
  ## ellipseS
  stat_ellipse(aes(group = group), level = 0.68, color = "grey25", linetype = 2, linewidth = 0.1) +
  ## sample points
  geom_point(alpha = 0.8) +
  ## group labels (custom positioned)
  geom_label(
    data = group_means_labeled,
    aes(
      x = label_x,
      y = label_y,
      label = group_name
    ),
    inherit.aes = FALSE,
    size = 3,
    color = "grey25",
    linewidth = 0.15,
    fill = "white"
  ) +
  theme_minimal() +
  scale_size_manual(values = c(
    "historical" = 3.2,
    "modern" = 2.5,
    "sillperioder" = 3
  )) +
  scale_fill_manual(values = c(
    "Iceland & Faroe Islands" = "#1b639e",
    "Norwegian Sea & Fjords" = "#02d97c",
    "Skagerrak-Kattegat" = "#b370b2",
    "North Sea" = "#e73f29",
    "United Kingdom" = "#a6721e",
    "Other" = "#7a7a7a"
  ))  +
  scale_color_manual(values = c(
    "historical" = "white",
    "modern" = "white",
    "sillperioder" = "grey25"
  ))  +
  scale_shape_manual(values = c(
    "historical" = 21,
    "modern" = 23,
    "sillperioder" = 21
  )) +
  scale_discrete_manual(
    aesthetics = "stroke",
    values = c(
      "historical" = 0.001,
      "modern" = 0.001,
      "sillperioder" = 0.6
    )
  ) +
  labs(title = "Discriminant Analysis of Principal Components (DAPC)",
       x = paste0("Discriminant Function 1 (", da_percent_variance[1], "%)"),
       y = paste0("Discriminant Function 2 (", da_percent_variance[2], "%)")) +
  theme(legend.position = "none")

p1

location_cols <- c(
  "Iceland & Faroe Islands" = "#1b639e",
  "Norwegian Sea & Fjords" = "#02d97c",
  "Skagerrak-Kattegat" = "#b370b2",
  "North Sea" = "#e73f29",
  "United Kingdom" = "#a6721e"
)

legend_locations <- ggplot(
  data.frame(
    location = names(location_cols),
    y = seq_along(location_cols)
  ),
  aes(x = 0.2, y = y)
) +
  geom_tile(
    aes(fill = location),
    width = 0.35,
    height = 0.35
  ) +
  geom_text(
    aes(x = 0.7, label = location),
    hjust = 0,
    size = 3
  ) +
  scale_fill_manual(values = location_cols) +
  expand_limits(x = c(0.5, 7)) +   # 🔑 reserve space
  labs(title = "Sampling Location") +
  coord_equal(clip = "off") +      # 🔑 square tiles, no squish
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.1, size = 10),
    plot.margin = margin(5, 5, 5, 5),
    plot.background = element_rect(fill = "white")
  )

# legend_locations

time_shapes <- c(
  "Historical" = 16,
  "Modern" = 18,
  "Sillperioder" = 21
)

legend_times <- ggplot(
  data.frame(
    time = names(time_shapes),
    y = seq_along(time_shapes)
  ),
  aes(x = 0.375, y = y)
) +
  geom_point(
    aes(shape = time),
    size = 2
  ) +
  geom_text(
    aes(x = 0.7, label = time),
    hjust = 0,
    size = 3
  ) +
  scale_shape_manual(values = time_shapes) +
  expand_limits(x = c(0.5, 7)) +   # 🔑 reserve space
  labs(title = "Sampling Time") +
  coord_equal(clip = "off") +      # 🔑 keeps vertical spacing consistent
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.1, size = 10),
    plot.margin = margin(5, 5, 5, 5),
    plot.background = element_rect(fill = "white")
  )

# legend_times

# 1️⃣ Determine number of rows in each legend
n_loc <- length(location_cols)
n_time <- length(time_shapes)

# 2️⃣ Add expand_limits to force the panels to full height
legend_locations <- legend_locations +
  expand_limits(y = c(0, n_loc + 1)) +
  theme(plot.margin = margin(5, 5, 5, 5))

legend_times <- legend_times +
  expand_limits(y = c(0, n_time + 1)) +
  theme(plot.margin = margin(5, 5, 5, 5))

# 3️⃣ Decide heights for stacking (normalized)
h1 <- 0.4   # top legend
h2 <- 0.3   # bottom legend
total_h <- h1 + h2
gap <- 0.05
bottom_y <- (1 - total_h - gap) / 2

legend_combined <- ggdraw() +
  draw_plot(legend_times, x = 0, y = bottom_y, width = 1, height = h2) +
  draw_plot(legend_locations, x = 0, y = bottom_y + h2 + gap, width = 1, height = h1)

# legend_combined

final_plot <- plot_grid(
  p1,
  legend_combined,
  ncol = 2,
  rel_widths = c(1, 0.3)
)

final_plot
ggsave("plots/angsd_matrix/dapc_wp1_final_sf7_sites.png",
       plot = final_plot,
       width = 8,
       height = 5,
       dpi = 300)
ggsave("plots/angsd_matrix/dapc_wp1_final_sf7_sites.pdf",
       plot = final_plot,
       width = 8,
       height = 5)

# Save DAPC object
saveRDS(DAPC, file = "data/angsd_matrix/DAPC_wp1_final_sf7_sites.rds")

####################################################################################

historical <- dapc_df %>%
  filter(time == "historical")

sillperioder <- dapc_df %>%
  filter(time == "sillperioder")

modern_uk <- dapc_df %>%
  filter(group_names == "United Kingdom" & time == "modern")

zz <- dapc_df |> 
  filter(time != "modern")

nosill_kattskag <- dapc_df %>%
  filter((time != "sillperioder" & loc == "Skagerrak-Kattegat"))

# https://mikkovihtakari.github.io/ggOceanMaps/index.html

post <- data.frame(DAPC$posterior)
rownames(post) <- sample_ids
post$loc <- get_locs(sample_ids)
post$time <- get_times(data_table, samples)
post$group <- DAPC$assign
post$group_names <- assign_group_name(post)
