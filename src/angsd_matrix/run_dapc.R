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

# run dapc:
run_dapc <- function(matrix, k) {
  grp <- find.clusters(matrix, n.pca = 3, n.clust = k)
  candidate_grps <- grp$grp
  xval <- xvalDapc(matrix, candidate_grps, n.pca.max = 100,
                   result = "groupMean",
                   n.pca = 3, n.rep = 100, xval.plot = FALSE)
  DAPC <- dapc(matrix, candidate_grps,
               n.pca = xval$DAPC$n.pca,
               n.da = xval$DAPC$n.da)
  return (DAPC)
}

# save_dapc:
save_dapc <- function(pcangsd_dataset, dapc_dataset, sites_name, k) {
  matrix <- get_matrix(pcangsd_dataset, sites_name, dapc_dataset)
  DAPC <- run_dapc(matrix, k)
  saveRDS(DAPC, file = paste0("data/angsd_matrix/DAPC.", pcangsd_dataset, ".", dapc_dataset, ".", sites_name, ".k", k, ".rds"))
}

# load_dapc:
load_dapc <- function(pcangsd_dataset, dapc_dataset, sites_name, k) {
  DAPC <- readRDS(paste0("data/angsd_matrix/DAPC.", pcangsd_dataset, ".", dapc_dataset, ".", sites_name, ".k", k, ".rds"))
  return (DAPC)
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

# assign groups from DAPC object and a vector of group names:
get_dapc_groups <- function(dapc, group_names, label_col = "Group") {
  if (!all(dapc$assign %in% seq_along(group_names))) {
    warning("Some assignments do not match group_names and will be set to NA")
  }
  groups <- group_names[dapc$assign]
  posterior <- as.data.frame(dapc$posterior)
  colnames(posterior) <- group_names
  lds <- data.frame(dapc$ind.coord)
  colnames(lds) <- paste0(label_col, "_", colnames(lds))
  out <- cbind(data.frame(groups), lds, posterior)
  names(out)[1] <- label_col
  return (out)
}

# dataframe of groups (ecotypes):
get_groups_df <- function(pcangsd_dataset, dapc_dataset){
  groups <- cbind(
    get_dapc_groups(
      load_dapc(pcangsd_dataset, dapc_dataset, "supplementary_file_7.v2", 7),
      c("Norway_1", "Norway_2", "North_Atlantic", "North_Sea",
        "Skagerrak_&_Kattegat", "Baltic", "Britain_&_Ireland"),
      label_col = "GeneticCluster"
    ),
    get_dapc_groups(
      load_dapc(pcangsd_dataset, dapc_dataset, "ns_inversions.chr12", 3),
      c("N/S", "N/N", "S/S"),
      label_col = "Inversion"
    ),
    get_dapc_groups(
      load_dapc(pcangsd_dataset, dapc_dataset, "baltic_v_atlantic.v2", 3),
      c("Fresh", "Transition", "Salty"),
      label_col = "Salinity"
    ),
    get_dapc_groups(
      load_dapc(pcangsd_dataset, dapc_dataset, "spring_v_autumn.v2", 3),
      c("Autumn", "Baltic Spring", "Norwegian Spring"),
      label_col = "SpawnTime"
    )
  )
}


## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## RUN:

pcangsd_dataset <- "full_herr"
dapc_dataset <- "wp1_final_bal"

# Save DAPC results - RUN ONLY ONCE
# save_dapc(pcangsd_dataset, dapc_dataset, "supplementary_file_7.v2", 7)
# save_dapc(pcangsd_dataset, dapc_dataset, "ns_inversions.chr12", 3)
# save_dapc(pcangsd_dataset, dapc_dataset, "baltic_v_atlantic.v2", 3)
# save_dapc(pcangsd_dataset, dapc_dataset, "spring_v_autumn.v2", 3)

# Load the saved DAPC results if needed:
sf7_dapc <- load_dapc(pcangsd_dataset, dapc_dataset, "supplementary_file_7.v2", 7)
inv_dapc <- load_dapc(pcangsd_dataset, dapc_dataset, "ns_inversions.chr12", 3)
bva_dapc <- load_dapc(pcangsd_dataset, dapc_dataset, "baltic_v_atlantic.v2", 3)
sva_dapc <- load_dapc(pcangsd_dataset, dapc_dataset, "spring_v_autumn.v2", 3)

# Load metadata table
metadata <- get_metadata(dapc_dataset)

# Get dataframe of Ecotypes
groups <- get_groups_df(pcangsd_dataset, dapc_dataset)

## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## PLOTS:

ggplot(
  data = dapc_df,
  aes(x = LD1, color = region, fill = region)
) +
  geom_density(alpha = 0.5, bw = 0.5) +
  scale_color_manual(values = c(
    "North_Atlantic" = "#1b639e",
    "Norway" = "#02d97c",
    "Skagerrak_&_Kattegat" = "#b370b2",
    "North_Sea" = "#e73f29",
    "Britain_&_Ireland" = "#a6721e",
    "BALTIC" = "#440154FF",
    "Other" = "#7a7a7a"
  )) +
  scale_fill_manual(values = c(
    "North_Atlantic" = "#1b639e",
    "Norway" = "#02d97c",
    "Skagerrak_&_Kattegat" = "#b370b2",
    "North_Sea" = "#e73f29",
    "Britain_&_Ireland" = "#a6721e",
    "BALTIC" = "#440154FF",
    "Other" = "#7a7a7a"
  )) +
  theme_minimal()
