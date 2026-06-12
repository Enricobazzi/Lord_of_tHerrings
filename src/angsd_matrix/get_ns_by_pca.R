library(tidyverse)

## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## FUNCTIONS:

# get sample names for a dataset from the sample_list file:
get_samples_from_dataset <- function(dataset) {
  file_path <- paste0("data/angsd_matrix/bamlists/", dataset, ".sample_list.txt")
  samples <- read.table(file_path)[, 1] |> as.character()
  return (samples)
}

# get matrix of the pcangsd dataset + sites, and filter dapc_dataset samples:
get_matrix <- function(pcangsd_dataset, sites_name, dapc_dataset) {
  # read matrix
  file_path <- paste0("data/angsd_matrix/pcangsd/", pcangsd_dataset, ".", sites_name, ".pcangsd.cov")
  mat <- as.matrix(read.table(file_path))
  # decide which samples to keep
  pcangsd_samples <- get_samples_from_dataset(pcangsd_dataset)
  dapc_samples <- get_samples_from_dataset(dapc_dataset)
  sample_indices <- which(pcangsd_samples %in% dapc_samples)
  # filter samples
  filtered_mat <- mat[sample_indices, sample_indices]
  rownames(filtered_mat) <- pcangsd_samples[sample_indices]
  colnames(filtered_mat) <- pcangsd_samples[sample_indices]
  return (filtered_mat)
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

# run pca
run_pca <- function(cov_matrix) {
  e <- eigen(cov_matrix)
  return(e)
}

# get vectors of pca
get_eigenvects <- function(eigen_result, num_vectors = 10) {
  vectors <- eigen_result$vectors[, 1:num_vectors]
  return(vectors)
}

# get variance
get_percent_variance <- function(cov_matrix, npcs = 10) {
  eigen_result <- run_pca(cov_matrix)
  eigenvalues <- eigen_result$values[1:npcs]
  percent_variance <- (eigenvalues / sum(eigen_result$values)) * 100
  names(percent_variance) <- paste0("PC", seq_len(npcs))
  return(percent_variance)
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
  "winter" = "#1b639e",
  "North" = "#1b639e",
  "South" = "#ED3911",
  "Het" = "#FED439FF"
)

# palette of shapes
time_shapes <- c(
  "Before 1700s" = 22,
  "1747–1805 Sillperiod" = 24,
  "Between Sillperiods" = 23,
  "1877–1906 Sillperiod" = 25,
  "Present-day" = 21
)

# name of the dataset used to generate matrix in PCANGSD
pcangsd_dataset <- "full_herr"
# name of the dataset used for DAPC analysis
dapc_dataset <- "wp1_final_bal"
# sites = inversion
sites_name <- "ns_inversions.chr12"

metadata <- get_metadata(dapc_dataset)
matrix <- get_matrix(pcangsd_dataset, sites_name, dapc_dataset)
npcs <- 2

eigen_result <- run_pca(matrix)
eigenvectors <- get_eigenvects(eigen_result, npcs)
pca_df <- as.data.frame(eigenvectors)
colnames(pca_df) <- paste0("PC", seq_len(ncol(pca_df)))

plt_df <- cbind(
  data.frame(
    sample = metadata$sample_id,
    id = metadata$new.id,
    Region = unlist(lapply(metadata$region, get_fancy_region)),
    Spawn = metadata$spawn,
    old_id = metadata$sample_id,
    Period = unlist(lapply(metadata$period, get_fancy_period)),
    Year = metadata$year,
    Inversion = ifelse(pca_df$PC1 > 0.06, "South", 
                       ifelse(pca_df$PC1 < 0, "North", "Het"))
    ),
  pca_df
)

pcvar <- get_percent_variance(matrix)

# density <- ggplot(data = plt_df, aes(x = PC1, color = Inversion, fill = Inversion)) +
#   geom_density() +
#   geom_vline(xintercept = c(0, 0.06), linewidth = 0.3) +
#   scale_fill_manual(values = c(colpal)) +
#   scale_color_manual(values = c(colpal)) +
#   theme_bw()
# density


pca <- ggplot(data = plt_df, aes(x = PC1, y = PC2, color = Region, fill = Region, shape = Period)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(xintercept = c(0, 0.06), linewidth = 0.3) +
  scale_fill_manual(values = c(colpal)) +
  scale_color_manual(values = c(colpal)) +
  scale_shape_manual(values = time_shapes) +
  xlab(paste0("PC1 (", round(pcvar[1], 1), "%)")) + ylab(paste0("PC2 (", round(pcvar[2], 1), "%)")) +
  theme_bw() +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0, "cm")
  )


ggsave(
  filename = paste0("plots/angsd_matrix/dapc/", dapc_dataset, ".", sites_name, ".pca.png"),
  plot = pca, width = 180, height = 90, unit = "mm", dpi = 300
)

write.table(plt_df,
            paste0("data/angsd_matrix/dapc/ns_inv_table.", pcangsd_dataset, ".", dapc_dataset, ".csv"),
            sep = ",", row.names = F, quote = F)

################################################################################

plot_freq <- function(plt_df, region){
  plt_df$Period <- factor(
    plt_df$Period,
    levels = c(
      "Before 1700s",
      "1747–1805 Sillperiod",
      "Between Sillperiods",
      "1877–1906 Sillperiod",
      "Present-day"
    ),
    ordered = TRUE
  )
  
  period_freq <- plt_df %>% filter(Region == region) |> 
    mutate(
      north = case_when(
        Inversion == "North" ~ 1,
        Inversion == "South" ~ 0,
        Inversion == "Het" ~ 0.5
      )
    ) %>%
    group_by(Period) %>%
    summarise(
      north_freq = mean(north),
      south_freq = 1 - north_freq,
      n = n(),
      .groups = "drop"
    )
  
  plot_df <- period_freq %>%
    pivot_longer(
      c(north_freq, south_freq),
      names_to = "allele",
      values_to = "frequency"
    )
  
  #ggplot(plot_df,
  #       aes(Period, frequency,
  #           color = allele,
  #           group = allele)) +
  #  geom_point(size = 3) +
  #  geom_line() +
  #  theme_bw() +
  #  theme(
  #    axis.text.x = element_text(angle = 45, hjust = 1)
  #  )
  ggplot(
    plot_df,
    aes(x = Period,
        y = frequency,
        fill = allele)
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
  
}
plot_freq(plt_df, "Skagerrak & Kattegat")
