library(tidyverse)
library(ggrepel)

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

run_pca <- function(cov_matrix) {
  e <- eigen(cov_matrix)
  return(e)
}

get_eigenvects <- function(eigen_result, num_vectors = 10) {
  vectors <- eigen_result$vectors[, 1:num_vectors]
  return(vectors)
}

# palette of colors
colpal <- c(
  "Vasa Ship" = "black",
  "Bothnia" = "#440154FF",
  "Baltic" = "#1b639e",
  "South Baltic" = "navyblue",
  "Switzerland" = "grey",
  "Finland Gulf" = "blue",
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

# get metadata table
get_metadata <- function(dapc_dataset) {
  samples <- get_samples_from_dataset(dapc_dataset)
  sample_data_file <- "~/Documents/Silly-periods/data/samples_table.csv"
  sample_data <- read.table(sample_data_file, sep = ",",
                            header = TRUE, na.strings = "UNKNOWN")
  sample_data <- sample_data[sample_data$sample_id %in% samples, ]
  return (sample_data)
}

# get fancy region names
get_fancy_region <- function(region) {
  return(gsub("_", " ", region))
}

build_pca_df <- function(matrix, metadata, npcs = 10) {
  eigen_result <- run_pca(matrix)
  eigenvectors <- get_eigenvects(eigen_result, npcs)
  pca_df <- as.data.frame(eigenvectors)
  colnames(pca_df) <- paste0("PC", seq_len(ncol(pca_df)))
  # Use rownames from the matrix to ensure correct order
  pca_df$id <- metadata$new.id
  pca_df$sample <- metadata$sample_id
  pca_df$region <- unlist(lapply(metadata$region, get_fancy_region))
  pca_df$region[metadata$location == "Vasa Ship"] <- "Vasa Ship"
  pca_df$spawn <- metadata$spawn
  pca_df$time <- ifelse(metadata$year > 1905, "Modern", "Historical")
  return(pca_df)
}

get_percent_variance <- function(cov_matrix, npcs = 10) {
  eigen_result <- run_pca(cov_matrix)
  eigenvalues <- eigen_result$values[1:npcs]
  percent_variance <- (eigenvalues / sum(eigen_result$values)) * 100
  names(percent_variance) <- paste0("PC", seq_len(npcs))
  return(percent_variance)
}

plot_pca <- function(pca_df, variance_explained, pc_x = 1, pc_y = 2, by = "region") {
  if (by == "spawn"){
    return(
      ggplot(pca_df,
             aes(x = !!sym(paste0("PC", pc_x)), y = !!sym(paste0("PC", pc_y)),
                 color = spawn, fill = spawn, shape = time)) +
        geom_point(alpha = 0.7, size = 2) +
        # labels
        geom_label_repel(data = pca_df |> filter(id %in% c("vasa_1", "vasa_2", "vasa_3", "vasa_4", "vasa_5")),
                         aes(x = !!sym(paste0("PC", pc_x)), y = !!sym(paste0("PC", pc_y)),
                             label = sample),
                         fill = "white", size = 2.5, inherit.aes = F,
                         label.padding = unit(0.1, "lines"),
                         min.segment.length = 0,
                         segment.color = "black",
                         segment.size = 0.5,
                         segment.alpha = 0.8
        ) +
        theme_minimal() +
        labs(x = paste0(paste0("PC", pc_x), " (", round(variance_explained[pc_x], 2), "%)"),
             y = paste0(paste0("PC", pc_y), " (", round(variance_explained[pc_y], 2), "%)")) +
        #     title = "PCA of samples based on covariance matrix") +
        scale_color_manual(values = colpal) +
        scale_fill_manual(values = colpal) +
        scale_shape_manual(values = c("Historical" = 21, "Modern" = 22))
    )
  } else {
    return(
      ggplot(pca_df,
             aes(x = !!sym(paste0("PC", pc_x)), y = !!sym(paste0("PC", pc_y)),
                 color = region, fill = region, shape = time)) +
        geom_point(alpha = 0.7, size = 2) +
        # labels
        geom_label_repel(data = pca_df |> filter(id %in% c("vasa_1", "vasa_2", "vasa_3", "vasa_4", "vasa_5")),
                         aes(x = !!sym(paste0("PC", pc_x)), y = !!sym(paste0("PC", pc_y)),
                             label = sample),
                         fill = "white", size = 2.5, inherit.aes = F,
                         label.padding = unit(0.1, "lines"),
                         min.segment.length = 0,
                         segment.color = "black",
                         segment.size = 0.5,
                         segment.alpha = 0.8
        ) +
        theme_minimal() +
        labs(x = paste0(paste0("PC", pc_x), " (", round(variance_explained[pc_x], 2), "%)"),
             y = paste0(paste0("PC", pc_y), " (", round(variance_explained[pc_y], 2), "%)")) +
        #     title = "PCA of samples based on covariance matrix") +
        scale_color_manual(values = colpal) +
        scale_fill_manual(values = colpal) +
        scale_shape_manual(values = c("Historical" = 21, "Modern" = 22))
    )
  }
}

# name of the dataset used to generate matrix in PCANGSD
pcangsd_dataset <- "full_herr"
# name of the dataset used for DAPC analysis
# grep "harengus" ../Silly-periods/data/samples_table.csv | 
#   grep -vE "NA|NW-ATL|Lamichh|Barrio|Switz|HER135" | 
#   awk -F',' '$9 == "current" || $10 >= 0.1' | 
#   grep -Ev 'ND([^3]|3([^678]|6[^45]|7[^45]|8[^7])|$)' |
#   cut -d',' -f1 > data/angsd_matrix/bamlists/vasa_ship.sample_list.txt
dapc_dataset <- "vasa_ship"
# name of the dataset used for DAPC analysis
# grep "harengus" ../Silly-periods/data/samples_table.csv | 
#   grep -vE "NA|NW-ATL|Lamichh|Barrio|Switz|HER135" | 
#   awk -F',' '$9 == "current" || $10 >= 0.1' |
#   cut -d',' -f1 > data/angsd_matrix/bamlists/vasa_ship_plus.sample_list.txt
dapc_dataset <- "vasa_ship_plus"
# name of the dataset used for DAPC analysis
# grep "harengus" ../Silly-periods/data/samples_table.csv | 
#   grep -vE "NA|NW-ATL|Lamichh|Barrio|Switz|HER135" | 
#   awk -F',' '$9 == "current" || $10 >= 0.1' |
#   grep -E "Baltic|Finland|Belt|Bothnia" |
#   grep -Ev 'ND([^3]|3([^678]|6[^45]|7[^45]|8[^7])|$)' |
#   cut -d',' -f1 > data/angsd_matrix/bamlists/vasa_ship_baltic.sample_list.txt
dapc_dataset <- "vasa_ship_baltic"
# sites
sites_name_lst <- c(
  "supplementary_file_7.v2",
  "baltic_v_atlantic.v2",
  "spring_v_autumn.v2",
  "ns_inversions.chr12"
)

for (sites_name in sites_name_lst) {
  metadata <- get_metadata(dapc_dataset)
  matrix <- get_matrix(pcangsd_dataset, sites_name, dapc_dataset)
  pca_df <- build_pca_df(matrix, metadata)
  variance_explained <- get_percent_variance(matrix)
  if (sites_name == "spring_v_autumn.v2") {
    p1 <- plot_pca(pca_df, variance_explained, by = "spawn")
    p1 <- p1 +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size = 7),
        legend.title =  element_blank(),
        legend.key.size = unit(0.15, "cm"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.spacing.y = unit(0.5, "cm"),
        legend.margin = margin(0,0,0,0)
      ) +
      guides(
        color = guide_legend(nrow = 1, order = 1),
        shape = guide_legend(nrow = 1,order = 2),
        fill = "none"
      )
  } else if (dapc_dataset == "vasa_ship_baltic") {
    p1 <- plot_pca(pca_df, variance_explained)
    p1 <- p1 +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size = 7),
        legend.title =  element_blank(),
        legend.key.size = unit(0.15, "cm"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.spacing.y = unit(0.5, "cm"),
        legend.margin = margin(0,0,0,0)
      ) +
      guides(
        color = guide_legend(nrow = 1, order = 1),
        shape = guide_legend(nrow = 1,order = 2),
        fill = "none"
      )
  } else {
    p1 <- plot_pca(pca_df, variance_explained)
    p1 <- p1 +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size = 7),
        legend.title =  element_blank(),
        legend.key.size = unit(0.15, "cm"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.spacing.y = unit(0.5, "cm"),
        legend.margin = margin(0,0,0,0)
      ) +
      guides(
        color = guide_legend(nrow = 3, order = 1),
        shape = guide_legend(nrow = 1,order = 2),
        fill = "none"
      )
  }
  
  ggsave(
    filename = paste0("plots/vasa_ship/", dapc_dataset, ".", sites_name, ".pca.png"),
    plot = p1, width = 120, height = 90, unit = "mm", dpi = 300
  )
}

###

pca_ns <- build_pca_df(
  matrix = get_matrix(
    pcangsd_dataset = "full_herr",
    sites_name = "ns_inversions.chr12",
    dapc_dataset = "vasa_ship"
  ),
  metadata = get_metadata(dapc_dataset = "vasa_ship")
)

pca_ns_variance <- get_percent_variance(
  cov_matrix = get_matrix(
    pcangsd_dataset = "full_herr",
    sites_name = "ns_inversions.chr12",
    dapc_dataset = "vasa_ship"
  )
)

p1 <- plot_pca(pca_ns, pca_ns_variance, pc_x = 1, pc_y = 3)
p1 <- p1 +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 7),
    legend.title =  element_blank(),
    legend.key.size = unit(0.15, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.margin = margin(0,0,0,0)
  ) +
  guides(
    color = guide_legend(nrow = 3, order = 1),
    shape = guide_legend(nrow = 1,order = 2),
    fill = "none"
  )

ggsave(
  filename = paste0("plots/vasa_ship/vasa_ship.ns_inversions.chr12.pca.pc1pc3.png"),
  plot = p1, width = 120, height = 90, unit = "mm", dpi = 300
)

###

plot_pc1_density <- function(pca_df, variance_explained, by = "region") {
  return(
    if (by == "spawn") {
      ggplot(data = pca_df, aes(x = PC1, color = spawn, fill = spawn)) +
        # geom_density(alpha = 0.5, linewidth = 0.3, position = "identity") +
        geom_density(
          aes(y = after_stat(density * n)),
          alpha = 0.5, linewidth = 0.3,
          position = "identity"
        ) + 
        scale_color_manual(values = colpal) +
        scale_fill_manual(values = colpal) +
        labs(x = paste0(paste0("PC1 (", round(variance_explained[1], 2), "%)")),
             y = "Density") +
        theme_minimal()
    } else {
      ggplot(data = pca_df, aes(x = PC1, color = region, fill = region)) +
        # geom_density(alpha = 0.5, linewidth = 0.3, position = "identity") +
        geom_density(
          aes(y = after_stat(density * n)),
          alpha = 0.5, linewidth = 0.3,
          position = "identity"
        ) + 
        scale_color_manual(values = colpal) +
        scale_fill_manual(values = colpal) +
        labs(x = paste0(paste0("PC1 (", round(variance_explained[1], 2), "%)")),
             y = "Density") +
        theme_minimal()
    }
  )
}

pca_spawn <- build_pca_df(
  matrix = get_matrix(
    pcangsd_dataset = "full_herr",
    sites_name = "spring_v_autumn.v2",
    dapc_dataset = "vasa_ship_baltic"
  ),
  metadata = get_metadata(dapc_dataset = "vasa_ship_baltic")
  )

pca_spawn_variance <- get_percent_variance(
  cov_matrix = get_matrix(
    pcangsd_dataset = "full_herr",
    sites_name = "spring_v_autumn.v2",
    dapc_dataset = "vasa_ship_baltic"
  )
)

pca_salinity <- build_pca_df(
  matrix = get_matrix(
    pcangsd_dataset = "full_herr",
    sites_name = "baltic_v_atlantic.v2",
    dapc_dataset = "vasa_ship"
  ),
  metadata = get_metadata(dapc_dataset = "vasa_ship")
)

pca_salinity_variance <- get_percent_variance(
  cov_matrix = get_matrix(
    pcangsd_dataset = "full_herr",
    sites_name = "baltic_v_atlantic.v2",
    dapc_dataset = "vasa_ship"
  )
)


pca_ns <- build_pca_df(
  matrix = get_matrix(
    pcangsd_dataset = "full_herr",
    sites_name = "ns_inversions.chr12",
    dapc_dataset = "vasa_ship"
  ),
  metadata = get_metadata(dapc_dataset = "vasa_ship")
)

pca_ns_variance <- get_percent_variance(
  cov_matrix = get_matrix(
    pcangsd_dataset = "full_herr",
    sites_name = "ns_inversions.chr12",
    dapc_dataset = "vasa_ship"
  )
)

plot_pc1_density(pca_df = pca_salinity, variance_explained = pca_salinity_variance) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 7),
    legend.title =  element_blank(),
    legend.key.size = unit(0.15, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.margin = margin(0,0,0,0),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

