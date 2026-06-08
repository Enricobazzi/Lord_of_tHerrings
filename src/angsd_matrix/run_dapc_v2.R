library(tidyverse)
library(adegenet)
library(clue)

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

# get the table with results of 1000 find.clusters runs
calculate_k_table <- function(matrix, n.reps = 1000) {
  Ks <- replicate(n.reps, {
    set.seed(sample(1:1e6, 1))
    length(find.clusters(matrix, n.pca = 5, choose.n.clust = FALSE)$size)
  })
  print(data.frame(table(Ks)) |> 
          arrange(-Freq) |> 
          mutate(Freq = paste0(round(Freq/sum(Freq)*100, 1), "%")))
  return (data.frame(table(Ks)) |> arrange(-Freq))
}

# get optimal k and write table of best 1000 K searches
get_optimal_k <- function(matrix, dapc_dataset, sites_name, n.reps = 1000) {
  k_table <- calculate_k_table(matrix, n.reps)
  write.table(
    k_table,
    file = paste0("data/angsd_matrix/dapc/", dapc_dataset, ".", sites_name, ".k_table.txt"),
    quote = F, row.names = F
  )
  optimal_k <- as.numeric(as.character(k_table[which.max(k_table$Freq), ]$Ks))
  return (optimal_k)
}

# get consensus group assignments for selected K
get_consensus_groups <- function(matrix, optimal_k,  n.reps = 1000) {
  
  # Create a matrix to hold the assignments (Rows = individuals, Columns = replicates)
  assignment_matrix <- matrix(NA, nrow = nrow(matrix), ncol = n.reps)
  
  for (i in 1:n.reps) {
    # find.clusters uses K-means internally
    temp_clusters <- find.clusters(matrix, n.pca = 5, n.clust = optimal_k, choose.n.clust = FALSE)
    assignment_matrix[, i] <- temp_clusters$grp
  }
  
  # Convert assignments to a matrix of dummy variables
  # Then calculate co-occurrence probabilities
  coassign_array <- array(0, dim = c(nrow(matrix), nrow(matrix), n.reps))
  
  for (i in 1:n.reps) {
    # cl_membership creates a binary classification matrix
    mem_matrix <- as.cl_membership(as.factor(assignment_matrix[, i]))
    M <- unclass(mem_matrix)   # numeric membership matrix: individuals x clusters
    coassign_array[, , i] <- M %*% t(M)
  }
  consensus_matrix <- apply(coassign_array, c(1,2), mean)
  
  # Convert to a distance object (1 - consensus probability)
  dist_matrix <- as.dist(1 - consensus_matrix)
  # Hierarchical clustering with Average Linkage (UPGMA)
  hc_consensus <- hclust(dist_matrix, method = "average")
  # Cut the tree to get exactly K groups
  final.groups <- cutree(hc_consensus, k = optimal_k)
  return (final.groups)
}

# run DAPC:
run_dapc <- function(matrix, consensus_groups){
  xval <- xvalDapc(matrix, consensus_groups, n.pca.max = 100,
                   result = "groupMean",
                   n.rep = 100, xval.plot = FALSE)
  final.dapc <- dapc(matrix, consensus_groups,
                     n.pca = xval$DAPC$n.pca,
                     n.da = xval$DAPC$n.da)
  return (final.dapc)
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


## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## RUN DAPC and SAVE RESULTS:

for (sites_name in sites_name_lst) {
  print(paste(pcangsd_dataset, dapc_dataset, sites_name))
  
  # read covariance matrix
  matrix <- get_matrix(pcangsd_dataset, sites_name, dapc_dataset)
  
  # calculate optimal k
  optimal_k <- get_optimal_k(matrix, dapc_dataset, sites_name)
  
  # calculate consensus group assignments
  consensus_groups <- get_consensus_groups(matrix, optimal_k)

  # run dapc with cross-validation
  final.dapc <- run_dapc(matrix, consensus_groups)
  
  # save dapc object
  saveRDS(
    final.dapc,
    file = paste0("data/angsd_matrix/dapc/final_dapc.",
                  pcangsd_dataset, ".", dapc_dataset, ".", sites_name, ".k", optimal_k, ".rds")
  )
}
