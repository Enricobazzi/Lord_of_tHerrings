library(tidyverse)
library(adegenet)
library(ggrepel)
library(cowplot)

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

# variance of LD1 & LD2
get_variance_label <- function(dapc_obj, axn) {
  return(paste0("LD", axn, " (", round(dapc1$eig[axn] / sum(dapc1$eig) * 100, 2), "%)"))
}

# plot dapc:
plot_dapc <- function(dapc_obj, metadata, color_by = "region") {

  dapc_df <- data.frame(
    id = metadata$new.id,
    region = metadata$region,
    spawn = metadata$spawn,
    LD1 = dapc_obj$ind.coord[, 1],
    LD2 = dapc_obj$ind.coord[, 2],
    group = droplevels(dapc_obj$assign),
    old_id = metadata$sample_id
  )
  
  # dataframe for centroids with names:
  centroids <- dapc_df %>%
    group_by(group) %>%
    summarise(
      mean_LD1 = median(LD1, na.rm = TRUE),
      mean_LD2 = median(LD2, na.rm = TRUE),
      group_name = first(group)
    )
  
  if (color_by == "region") {
    dapc_plot <- ggplot (data = dapc_df,
                         aes(x = LD1, y = LD2, colour = region, fill = region))
  } else if (color_by == "spawn") {
    dapc_plot <- ggplot (data = dapc_df,
                         aes(x = LD1, y = LD2, colour = spawn, fill = spawn))
  }
  
 dapc_plot <- dapc_plot +
    # lines to group centroid
    geom_segment(
      data = dapc_df %>% left_join(centroids, by = "group"),
      aes(x = LD1, y = LD2, xend = mean_LD1, yend = mean_LD2),
      linewidth = 0.3,
      alpha = 0.8
    ) +
    # ellipse
    # stat_ellipse(aes(group = group), level = 0.5, linetype = 1, linewidth = 0.3) +
    # points
    geom_point(alpha = 0.7, size = 1.8, shape = 19) +
    geom_point(
      data = dapc_df |> filter(id %in% c("vasa_1", "vasa_2", "vasa_3", "vasa_4", "vasa_5")),
      aes(x = LD1, y = LD2), inherit.aes = F, shape = 16
    ) +
    # labels
    geom_label_repel(data = dapc_df |> filter(id %in% c("vasa_1", "vasa_2", "vasa_3", "vasa_4", "vasa_5")),
                     aes(x = LD1, y = LD2, label = old_id),
                    fill = "white", size = 1, inherit.aes = F,
                    label.padding = unit(0.1, "lines"), linewidth = 0.05,
                    min.segment.length = 0,
                    segment.color = "black",
                    segment.size = 0.5,
                    segment.alpha = 0.8
    ) +
    geom_label_repel(data = centroids, aes(x = mean_LD1, y = mean_LD2, label = group),
                    fill = "white", size = 1.5, inherit.aes = F,
                    label.padding = unit(0.1, "lines"), linewidth = 0.05, force = 0.01) +
    # fancy plot:
    xlab(get_variance_label(dapc_obj, 1)) + ylab(get_variance_label(dapc_obj, 2)) +
    scale_fill_manual(values = c(
      '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000'
      
    )) +
    scale_color_manual(values = c(
      '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000'
    )) +
    theme_minimal(base_size = 7) # + theme(legend.position = "none")
  return(dapc_plot)
}

pcangsd_dataset <- "full_herr"
# cat ../Silly-periods/data/samples_table.csv | grep -E "mh|ND364|ND365|ND374|ND375|ND387" | \
# cut -d',' -f1 > data/angsd_matrix/bamlists/vasa_modern.sample_list.txt
# dapc_dataset <- "vasa_modern"
dapc_dataset <- "full_herr"
metadata <- get_metadata(dapc_dataset)

# 1. Genetic structure / location ?
sites_name <- "supplementary_file_7.v2"
k <- 8
mat <- get_matrix(pcangsd_dataset, sites_name, dapc_dataset)
dapc1 <- run_dapc(mat, k)
p1 <- plot_dapc(dapc1, metadata)
p1 <- p1 +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title =  element_blank(),
    legend.key.size = unit(0.15, "cm"),
    legend.spacing.x = unit(0.01, "cm"),
    legend.spacing.y = unit(0.01, "cm"),
    legend.margin = margin(0,0,0,0)
  ) +
  guides(color = guide_legend(nrow = 1))

ggsave(
  filename = paste0("plots/vasa_ship/", dapc_dataset, ".", sites_name, ".k", k, ".dapc.png"),
  plot = p1, width = 90, height = 90, unit = "mm", dpi = 300
)
# 2. Autumn or Spring?
sites_name <- "spring_v_autumn.v2"
k <- 3
mat <- get_matrix(pcangsd_dataset, sites_name, dapc_dataset)
dapc2 <- run_dapc(mat, k)
p2 <- plot_dapc(dapc2, metadata, "spawn")
p2
ggsave(
  filename = paste0("plots/vasa_ship/", dapc_dataset, ".", sites_name, ".k", k, ".dapc.png"),
  plot = p2, height = 4, width = 8
)


p1 +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title =  element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.spacing.x = unit(0.05, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(0,0,0,0)
  ) +
  guides(color = guide_legend(nrow = 1))
