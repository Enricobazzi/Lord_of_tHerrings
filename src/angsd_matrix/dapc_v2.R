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
            "Britain & Ireland",
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

# # prior on first 4 (>90% of variance explained)
# grp <- find.clusters(matrix, n.pca = 4, n.clust = 6)
# grp
# candidate_grps <- grp$grp
# xval <- xvalDapc(matrix, candidate_grps, n.pca.max = 100,
#                  result = "groupMean",
#                  n.pca = NULL, n.rep = 50, xval.plot = FALSE)
# 
# # n.pca and n.da chosen based on xval results (10 and 4)
# DAPC <- dapc(matrix, candidate_grps,
#              n.pca = xval$DAPC$n.pca,
#              n.da = xval$DAPC$n.da)
# # Save DAPC object
# saveRDS(DAPC, file = "data/angsd_matrix/DAPC_wp1_final_sf7_sites.rds")

## if loading from saved DAPC object:
DAPC <- readRDS("data/angsd_matrix/DAPC_wp1_final_sf7_sites.rds")
# scatter(DAPC)

# biplot of dapc - prepare data frame:
dapc_df <- data.frame(DAPC$ind.coord)
old_ids <- rownames(dapc_df)
rownames(dapc_df) <- sample_ids
dapc_df$loc <- get_locs(sample_ids)
dapc_df$time <- get_times(data_table, samples)
dapc_df$group <- DAPC$assign
dapc_df$group_names <- assign_group_name(dapc_df)
dapc_df$old_id <- old_ids

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
  label_x = c(  -4.2,   1.5, 4.5,    5, -3.5),
  label_y = c(  -3.5,  -6.5, 4.5,   -5,  4.5)
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
    text.color = "black",
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
    "Britain & Ireland" = "#a6721e",
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
  "Britain & Ireland" = "#a6721e"
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
ggsave(
  filename = paste0("plots/angsd_matrix/legend_locations.pdf"),
  plot = legend_locations,
  height = 3, width = 5
)
ggsave(
  filename = paste0("plots/angsd_matrix/legend_locations.png"),
  plot = legend_locations,
  height = 3, width = 5
)

legend_times <- legend_times +
  expand_limits(y = c(0, n_time + 1)) +
  theme(plot.margin = margin(5, 5, 5, 5))
ggsave(
  filename = paste0("plots/angsd_matrix/legend_times.pdf"),
  plot = legend_times,
  height = 3, width = 5
)
ggsave(
  filename = paste0("plots/angsd_matrix/legend_times.png"),
  plot = legend_times,
  height = 3, width = 5
)


# 3️⃣ Decide heights for stacking (normalized)
h1 <- 0.4   # top legend
h2 <- 0.3   # bottom legend
total_h <- h1 + h2
gap <- 0.05
bottom_y <- (1 - total_h - gap) / 2

legend_combined <- ggdraw() +
  draw_plot(legend_times, x = 0, y = bottom_y, width = 1, height = h2) +
  draw_plot(legend_locations, x = 0, y = bottom_y + h2 + gap, width = 1, height = h1)

final_plot <- plot_grid(
  p1,
  legend_combined,
  ncol = 2,
  rel_widths = c(1, 0.2)
)

final_plot
ggsave("plots/angsd_matrix/dapc_wp1_final_sf7_sites.png",
       plot = final_plot,
       width = 14*0.67,
       height = 6*0.67,
       dpi = 300)
ggsave("plots/angsd_matrix/dapc_wp1_final_sf7_sites.pdf",
       plot = final_plot,
       width = 14*0.67,
       height = 6*0.67)

# add posterior probability for groups to dapc_df
post <- data.frame(DAPC$posterior)
# transform number to group name
colnames(post) <- substr(colnames(post), 2, 2)
rownames(post) <- get_sample_ids(data_table, rownames(post))
# check if rownames of post and dapc_df are the same
all(rownames(post) == rownames(dapc_df))
dapc_df <- cbind(dapc_df, post)

# Save the dapc_df
write.table(
  dapc_df,
  file = "data/angsd_matrix/dapc_df_wp1_final_sf7_sites.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE
)

####################################################################################

# plot DAPC of:
# modern only
modern_dapc <- ggplot(
  dapc_df %>% filter(time == "modern"),
  aes(x = LD1, y = LD2,
      fill = loc, color = time, shape = time,
      size = time, stroke = time)) +
  ## lines from samples to group centroid
  geom_segment(
    data = dapc_df %>% filter(time == "modern") %>% left_join(group_means, by = "group"),
    aes(x = LD1, y = LD2, xend = mean_LD1, yend = mean_LD2),
    inherit.aes = FALSE,
    color = "grey70",
    linewidth = 0.3,
    alpha = 0.8
  ) +
  ## lines from centroid to the group labels
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
  geom_point(alpha = 1) +
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
    text.color = "black",
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
    "Britain & Ireland" = "#a6721e",
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
  labs(subtitle = "Modern", x = "LD1", y = "LD2") +
  theme(legend.position = "none")

modern_dapc

# historical + modern low-alpha
historical_dapc <- ggplot(
  dapc_df %>% filter(time == "modern" | time == "historical"),
  aes(x = LD1, y = LD2,
      fill = loc, color = time, shape = time,
      size = time, stroke = time)) +
  ## lines from samples 1 to group centroid
  geom_segment(
    data = dapc_df %>% filter(time == "historical") %>% left_join(group_means, by = "group"),
    aes(x = LD1, y = LD2, xend = mean_LD1, yend = mean_LD2),
    inherit.aes = FALSE,
    color = "grey70",
    linewidth = 0.3,
    alpha = 0.8
  ) +
  ## lines from samples 2 to group centroid
  geom_segment(
    data = dapc_df %>% filter(time == "modern") %>% left_join(group_means, by = "group"),
    aes(x = LD1, y = LD2, xend = mean_LD1, yend = mean_LD2),
    inherit.aes = FALSE,
    color = "grey70",
    linewidth = 0.3,
    alpha = 0.2
  ) +
  ## lines from centroid to the group labels
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
  ## sample points 2
  geom_point(
    data = dapc_df %>% filter(time == "modern"),
    aes(x = LD1, y = LD2,
        fill = loc, color = time, shape = time,
        size = time, stroke = time), alpha = 0.2) +
  ## sample points 1
  geom_point(
    data = dapc_df %>% filter(time == "historical"),
    aes(x = LD1, y = LD2,
        fill = loc, color = time, shape = time,
        size = time, stroke = time), alpha = 1) +
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
    text.color = "black",
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
    "Britain & Ireland" = "#a6721e",
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
  labs(subtitle = "Historical", x = "LD1", y = "LD2") +
  theme(legend.position = "none")

historical_dapc

# sillperioder + modern low-alpha

sillperioder_dapc <- ggplot(
  dapc_df %>% filter(time == "modern" | time == "sillperioder"),
  aes(x = LD1, y = LD2,
      fill = loc, color = time, shape = time,
      size = time, stroke = time)) +
  ## lines from samples 1 to group centroid
  geom_segment(
    data = dapc_df %>% filter(time == "sillperioder") %>% left_join(group_means, by = "group"),
    aes(x = LD1, y = LD2, xend = mean_LD1, yend = mean_LD2),
    inherit.aes = FALSE,
    color = "grey70",
    linewidth = 0.3,
    alpha = 0.8
  ) +
  ## lines from samples 2 to group centroid
  geom_segment(
    data = dapc_df %>% filter(time == "modern") %>% left_join(group_means, by = "group"),
    aes(x = LD1, y = LD2, xend = mean_LD1, yend = mean_LD2),
    inherit.aes = FALSE,
    color = "grey70",
    linewidth = 0.3,
    alpha = 0.2
  ) +
  ## lines from centroid to the group labels
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
  ## sample points 2
  geom_point(
    data = dapc_df %>% filter(time == "modern"),
    aes(x = LD1, y = LD2,
        fill = loc, color = time, shape = time,
        size = time, stroke = time), alpha = 0.2) +
  ## sample points 1
  geom_point(
    data = dapc_df %>% filter(time == "sillperioder"),
    aes(x = LD1, y = LD2,
        fill = loc, color = time, shape = time,
        size = time, stroke = time), alpha = 1) +
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
    text.color = "black",
    linewidth = 0.15,
    fill = "white"
  ) +
  theme_minimal() +
  scale_size_manual(values = c(
    "historical" = 3.2,
    "modern" = 2.5,
    "sillperioder" = 3.2
  )) +
  scale_fill_manual(values = c(
    "Iceland & Faroe Islands" = "#1b639e",
    "Norwegian Sea & Fjords" = "#02d97c",
    "Skagerrak-Kattegat" = "#b370b2",
    "North Sea" = "#e73f29",
    "Britain & Ireland" = "#a6721e",
    "Other" = "#7a7a7a"
  ))  +
  scale_color_manual(values = c(
    "historical" = "white",
    "modern" = "white",
    "sillperioder" = "white"
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
      "sillperioder" = 0.001
    )
  ) +
  labs(subtitle = "Sillperioder", x = "LD1", y = "LD2") +
  theme(legend.position = "none")

sillperioder_dapc

## add maps
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
  x = c(-11,   15,  10, -15,  18),
  y = c(64.5,  67,  54,  55, 58)
)
sillperioder_coords <- data.frame(
  loc = c("Iceland&FaroeIslands", "NorwegianSea&Fjords", "NorthSea", "Britain&Ireland", "Skagerrak-Kattegat"),
  x = c(-11,   0,   10, -15, 18),
  y = c(64.5,  65,  52,  55, 58)
)
historical_coords <- data.frame(
  loc = c("Iceland&FaroeIslands", "NorwegianSea&Fjords", "NorthSea", "Britain&Ireland", "Skagerrak-Kattegat"),
  x = c(-11,   0,  10,  -15,  18),
  y = c(64.5,  65,  50,  55, 58)
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

list_of_map_plots <- list()
i <- 1
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
    geom_sf(data = time_locations, color = "grey15", shape = 1, size = 2, stroke = 0.6) +
    geom_sf(data = loc_rectangles, fill = NA, color = "grey15", size = 0.3) +
    geom_sf(data = lines_to_edge,
            color = "grey15", linewidth = 0.3) +
    coord_sf(
      xlim = xl, ylim = yl, expand = FALSE, default_crs = NULL
    ) +
    # add title in the middle of the plot
    # ggtitle(label = paste0(str_to_title(time_period))) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      plot.margin = margin(30, 5, 0, 0)
    )
  
  for (location in unique(time_df$loc)){
    # add pie charts
    filename <- paste0("plots/angsd_matrix/pie_charts/", time_period, ".", location, ".png")
    popx <- loc_coords$x[loc_coords$loc == location]
    popy <- loc_coords$y[loc_coords$loc == location]
    basemap2 <- basemap2 +
      draw_image(filename, x = popx, y = popy, vjust = 0.5, hjust = 0.5, scale = 9) 
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
  list_of_map_plots[[i]] <- basemap2
  i <- i + 1 
}



## FINAL PLOT
final_plot <- plot_grid(
  ggdraw() +
    draw_plot(modern_dapc, x = 0, y = 0, width = 0.6, height = 1) +
    draw_plot(list_of_map_plots[[1]], x = 0.6, y = 0, width = 0.4, height = 1),
  ggdraw() +
    draw_plot(historical_dapc, x = 0, y = 0, width = 0.6, height = 1) +
    draw_plot(list_of_map_plots[[2]], x = 0.6, y = 0, width = 0.4, height = 1),
  ggdraw() +
    draw_plot(sillperioder_dapc, x = 0, y = 0, width = 0.6, height = 1) +
    draw_plot(list_of_map_plots[[3]], x = 0.6, y = 0, width = 0.4, height = 1),
  ncol = 1,
  labels = "AUTO"
)

ggsave(
  filename = paste0("plots/angsd_matrix/dapc.png"),
  plot = final_plot,
  height = 9, width = 8
)

ggsave(
  filename = paste0("plots/angsd_matrix/dapc.pdf"),
  plot = final_plot,
  height = 9, width = 8
)

