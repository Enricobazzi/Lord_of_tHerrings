library(tidyverse)
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

all_samples_name <- "wp1_all"
subset_samples_name <- "wp1_modern"
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
rownames(matrix) <- sample_ids
colnames(matrix) <- sample_ids

library(adegenet)

grp <- find.clusters(matrix)

candidate_grps <- grp$grp

xval <- xvalDapc(matrix, candidate_grps, n.pca.max = 10,
                 result = "groupMean",
                 n.pca = NULL, n.rep = 10, xval.plot = TRUE)

DAPC <- dapc(matrix, candidate_grps,
             n.pca = xval$DAPC$n.pca,
             n.da = xval$DAPC$n.da)
scatter(DAPC)
dapc_matrix <- data.frame(DAPC$ind.coord)
dapc_post <- data.frame(DAPC$posterior)
