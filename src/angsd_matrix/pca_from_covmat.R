library(tidyverse)
library(ggplot2)
library(ggrepel)

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
    stringr::str_split(sample_ids, "_", simplify = TRUE)[, 1] %in% c("iceland", "faroe",
                                                             "norwegian", "more",
                                                             "haugesund", "karmoy",
                                                             "rovaer", "bergen",
                                                             "stavanger", "foldfjorden",
                                                             "norway"),
    "NE-Atlantic",
    ifelse(
      stringr::str_split(sample_ids, "_", simplify = TRUE)[, 1] %in% c("idefjord", "maseskar",
                                                               "kaupang", "risor",
                                                               "fehmarn", "koster",
                                                               "selso", "dynekilen",
                                                               "kampinge", "knastorp",
                                                               "nyalodose", "masthugget"),
      "Skagerrak-Kattegat",
      ifelse(
        stringr::str_split(sample_ids, "_", simplify = TRUE)[, 1] %in% c("celtic", "downs",
                                                                 "isleofman", "northsea",
                                                                 "netherlands", "york",
                                                                 "lyminge", "scotland"),
        "North Sea",
        "Other"
      )
    )
  )
  return(loc)
}

run_pca <- function(cov_matrix) {
  e <- eigen(cov_matrix)
  return(e)
}

get_eigenvects <- function(eigen_result, num_vectors = 10) {
  vectors <- eigen_result$vectors[, 1:num_vectors]
  return(vectors)
}

build_pca_df <- function(cov_matrix, sample_ids, npcs = 10) {
  eigen_result <- run_pca(cov_matrix)
  eigenvectors <- get_eigenvects(eigen_result, npcs)
  pca_df <- as.data.frame(eigenvectors)
  colnames(pca_df) <- paste0("PC", seq_len(ncol(pca_df)))
  pca_df$sample_id <- sample_ids
  pca_df$location <- get_locs(sample_ids)
  return(pca_df)
}

get_percent_variance <- function(cov_matrix, npcs = 10) {
  eigen_result <- run_pca(cov_matrix)
  eigenvalues <- eigen_result$values[1:npcs]
  percent_variance <- (eigenvalues / sum(eigen_result$values)) * 100
  names(percent_variance) <- paste0("PC", seq_len(npcs))
  return(percent_variance)
}

plot_pca <- function(pca_df, variance_explained, pc_x = 1, pc_y = 2) {
  ggplot(pca_df,
    aes(x = !!sym(paste0("PC", pc_x)), y = !!sym(paste0("PC", pc_y)),
        color = location, label = sample_id)) +
    geom_label_repel(
      size = 2, max.overlaps = 500, label.size = 0.1, box.padding = 0.1,
      label.padding = 0.1, force = 10,
    ) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(x = paste0("PC1 (", round(variance_explained[pc_x], 2), "%)"),
         y = paste0("PC2 (", round(variance_explained[pc_y], 2), "%)"),
         title = "PCA of samples based on covariance matrix") +
    scale_color_brewer(palette = "Set1")
}

all_samples_name <- "wp1_all"
subset_samples_name <- "wp1_modern_masthugget"
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
pca_df <- build_pca_df(matrix, sample_ids, npcs = 10)
variance_explained <- get_percent_variance(matrix, npcs = 10)

pca_plot <- plot_pca(pca_df, variance_explained, pc_x = 1, pc_y = 2)

ggsave(pca_plot,
  filename = paste0("plots/angsd_matrix/", all_samples_name, ".", subset_samples_name, ".", sites_name, ".pdf"),
  width = 12, height = 10
)
