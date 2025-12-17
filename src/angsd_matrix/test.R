library(tidyverse)
library(ggrepel)

samples <- read.table("data/angsd_matrix/bamlists/wp1_subset.sample_list.txt", header = FALSE)$V1
# remove HER135 from samples
samples <- samples[samples != "HER135"]
data_table <- read.table("data/samples_table.csv", sep = ",", header = TRUE, na.strings = "UNKNOWN")
data_table <- data_table[data_table$sample_id %in% samples, ]
sample_ids <- data_table$new.id[match(samples, data_table$sample_id)]

## PCA

name <- "data/angsd_matrix/output/wp1_subset.filtered.pcangsd.cov"
m <- as.matrix(read.table(name))
rownames(m) <- sample_ids
colnames(m) <- sample_ids
e <- eigen(m)
eigenvecs <- as.data.frame(e$vectors)
rownames(eigenvecs) <- sample_ids
eigenvecs$sample_id <- sample_ids
eigenvecs$pop <- str_split(sample_ids, "_", simplify = TRUE)[,1]
eigenvecs$loc <- ifelse(
  str_split(sample_ids, "_", simplify = TRUE)[,1] %in% c("iceland", "faroe", "norwegian"),
  "NE-Atlantic",
  ifelse(
    str_split(sample_ids, "_", simplify = TRUE)[,1] %in% c("idefjord", "maseskar", "risor"),
    "Skagerrak",
    ifelse(
      str_split(sample_ids, "_", simplify = TRUE)[,1] %in% c("celtic", "downs", "isleofman", "northsea"),
      "North Sea",
      "Other"
    )
  )
)

# plot pca color by loc
ggplot(eigenvecs, aes(x = V1, y = V2, color = loc, label = sample_id)) +
  geom_point(size = 3) +
  geom_label_repel(size = 2, max.overlaps = 200, label.size = 0.1, box.padding = 0.1, label.padding = 0.1, force = 20,
  ) +
  theme_minimal() +
  labs(x = paste0("PC1 (", round((e$values[1] / sum(e$values)) * 100, 2), "%)"),
       y = paste0("PC2 (", round((e$values[2] / sum(e$values)) * 100, 2), "%)"),
       title = "PCA of samples based on covariance matrix") +
  scale_color_brewer(palette = "Set1")


#neighbor joining
cor_mat <- cov2cor(m)
dist_mat <- as.dist(1 - cor_mat)
plot(ape::nj(dist_mat), use.edge.length=FALSE)
plot(hclust(dist(m), "ave"))

# root tree on node including norwegian 2 kalvsund 3 koster 4 koster 3 masthugget 2 norwegian 13 faroe 9 iceland 1 iceland 24 faroe 13 iceland 17 iceland 26 norwegian 28 norwegian 27 masthugget 8 masthugget 3 koster 2 gullholmen 2 masthugget 6 northsea 16 kalvsund 1 gullholmen 5 faroe 3 faroe 6 gullholmen 1 gullholmen 4 gullholmen 3 masthugget 7 masthugget 5 masthugget 4
tree <- ape::nj(dist_mat)
plot(tree, "u")
rooted_tree <- ape::root(tree, outgroup = c("norwegian_2", "kalvsund_3", "koster_4", "koster_3", "masthugget_2", "norwegian_13", "faroe_9", "iceland_1", "iceland_24", "faroe_13", "iceland_17", "iceland_26", "norwegian_28", "norwegian_27", "masthugget_8", "masthugget_3", "koster_2", "gullholmen_2", "masthugget_6", "northsea_16", "kalvsund_1", "gullholmen_5", "faroe_3", "faroe_6", "gullholmen_1", "gullholmen_4", "gullholmen_3", "masthugget_7", "masthugget_5", "masthugget_4"), resolve.root = TRUE)
plot(rooted_tree, use.edge.length=FALSE)
