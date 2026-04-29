library(tidyverse)
library(ggplot2)
library(ggrepel)
library(adegenet)
library(cowplot)

## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## CONSTANTS:


## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## FUNCTIONS:

# get sample names for a dataset from the sample_list file
get_samples_from_dataset <- function(dataset){
  file_path <- paste0("data/angsd_matrix/bamlists/", dataset, ".sample_list.txt")
  samples <- read.table(file_path)[, 1] |> as.character()
  return(samples)
}

# get matrix of the pcangsd dataset + sites, and filter dapc_dataset samples
get_matrix <- function(pcangsd_dataset, sites_name, dapc_dataset){
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
  return(filtered_mat)
}

# get samples metadata table
get_metadata <- function(samples){
  sample_data_file <- "~/Documents/Silly-periods/data/samples_table.csv"
  sample_data <- read.table(sample_data_file, sep = ",",
                            header = TRUE, na.strings = "UNKNOWN")
  sample_data <- sample_data[sample_data$sample_id %in% samples, ]
  return(sample_data)
}

# get new samples ids from the metadata table
get_sample_ids <- function(metadata, samples) {
  sample_ids <- metadata$new.id[match(samples, metadata$sample_id)]
  return(sample_ids)
}

# get group names from dapc dataframe
get_group_names <- function(df) {
  group_names <- df %>%
    group_by(group) %>%
    mutate(group_name = names(sort(table(region), decreasing = TRUE))[1]) %>%
    ungroup() %>%
    pull(group_name)
  return(group_names)
}

## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## VARIABLES/ARGUMENTS:

# dataset of full matrix generated with PCangsd
pcangsd_dataset <- "full_herr"
# dataset of samples to be analyzed in DAPC
dapc_dataset <- "wp1_final_bal"
# name of sites
sites_name <- "supplementary_file_7.v2"
# number of clusters (k)
k <- 7

# samples
samples <- get_samples_from_dataset(dapc_dataset)
# matrix
matrix <- get_matrix(pcangsd_dataset, sites_name, dapc_dataset)
# metadata
metadata <- get_metadata(samples)
# sample ids
sample_ids <- get_sample_ids(metadata, samples)

## DAPC
grp <- find.clusters(matrix, n.pca = 3, n.clust = k)
# grp <- find.clusters(matrix, n.pca = 3, choose.n.clust = F, criterion = "diffNgroup")
# usually k = 7 with this method
candidate_grps <- grp$grp
xval <- xvalDapc(matrix, candidate_grps, n.pca.max = 100,
                 result = "groupMean",
                 n.pca = 3, n.rep = 100, xval.plot = FALSE)
# n.pca and n.da chosen based on xval results (10 and 4)
DAPC <- dapc(matrix, candidate_grps,
             n.pca = xval$DAPC$n.pca,
             n.da = xval$DAPC$n.da)
# scatter(DAPC)
# # Save DAPC object
saveRDS(DAPC, file = paste0("data/angsd_matrix/DAPC.", pcangsd_dataset, ".", dapc_dataset, ".", sites_name, ".k", k, ".rds"))
# # Load the DAPC object
DAPC <- readRDS(paste0("data/angsd_matrix/DAPC.", pcangsd_dataset, ".", dapc_dataset, ".", sites_name, ".k", k, ".rds"))

# biplot of dapc - prepare data frame:
dapc_df <- data.frame(DAPC$ind.coord)
dapc_df$sample_id <- sample_ids
dapc_df$region <- metadata$region
dapc_df$period <- metadata$period

# add posterior probability for groups to dapc_df
post <- data.frame(DAPC$posterior)
# transform number to group name
# colnames(post) <- paste0("W", colnames(post))
dapc_df <- cbind(dapc_df, post)
# group
dapc_df$group <- DAPC$assign
dapc_df$groupname <- get_group_names(dapc_df)

# 
ggplot(
  data = dapc_df,
  aes(x = LD1, y = LD2, fill = region, shape = period)
  ) +
  geom_point(size = 3) +
  scale_fill_manual(values = c(
    "North_Atlantic" = "#1b639e",
    "Norway" = "#02d97c",
    "Skagerrak_&_Kattegat" = "#b370b2",
    "North_Sea" = "#e73f29",
    "Britain_&_Ireland" = "#a6721e",
    "BALTIC" = "#440154FF",
    "Other" = "#7a7a7a"
  )) +
  scale_shape_manual(values= c(
    "mh" = 23,
    "17sp" = 25,
    "18sp" = 24,
    "18rh" = 22,
    "ah" = 21
  )) +
  theme_minimal()


##

for (p in unique(dapc_df$period)) {
  for (r in unique(dapc_df$region)) {
    
    time_df <- dapc_df %>% filter(period == p, region == r)
    
    if (nrow(time_df) > 0) {
      
      a <- time_df |>
        gather(post_group, posterior, X1:X6) |>
        select(region, post_group, posterior)
      P <- a |>
        filter(region == r) |>
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
        geom_col(width = 0.1, color = "grey20", linewidth = 2, linejoin = "bevel") +
        coord_polar(theta = "y", start = 0) +
        scale_fill_manual(values = c(
          "X5" = "#1b639e",
          "X4" = "#02d97c",
          "X6" = "#3f704d",
          "X2" = "#b370b2",
          "X1" = "#e73f29",
          "X3" = "#a6721e",
          "X7" = "#440154FF",
          "Other" = "#7a7a7a"
        )) +
        theme_void() +
        theme(legend.position = "none")
    ggsave(
      filename = paste0("plots/angsd_matrix/pie_charts/", p, ".", r, ".png"),
      plot = PPP,
      height = 3, width = 3, device = png
    )
    }
  }
}

## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## ---------- ## ---------- ## ---------- ## ---------- ## ---------- ## -------
## BONUS SHIT:
## to get datasets for GONE2:

sk_mh <- dapc_df |> filter(X2 > 0.99, region == "Skagerrak_&_Kattegat", period == "mh")
print(rownames(sk_mh))
sk_18rh <- dapc_df |> filter(X2 > 0.99, region == "Skagerrak_&_Kattegat", period == "18rh")
print(rownames(sk_18rh))
