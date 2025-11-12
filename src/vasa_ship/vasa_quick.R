library(tidyverse)
library(ggrepel)

##### PCA #####
sample <- "ND364"
vec <- paste0("data/vasa_quick/vasaship_quick.", sample,".eigenvec")
val <- paste0("data/vasa_quick/vasaship_quick.", sample,".eigenval")

# read eigenvalues
eigenvals <- read.table(val,
                        col.names = c("eigenval"))

# read eigenvectors
eigenvecs <- read.table(vec,
                        header = TRUE)[,-1]

percents <- eigenvals$eigenval/sum(eigenvals$eigenval)*100

# scree-plot of eigenvalues
scree <- ggplot() +
  geom_line(aes(x = as.numeric(row.names(eigenvals)), y = percents)) +
  geom_point(aes(x = as.numeric(row.names(eigenvals)), y = percents),
             fill = "grey", shape = 21, size = 2.5) +
  xlab(paste0("Number of Principal Component")) +
  ylab(paste0("Percentage of variance Explained")) +
  scale_x_continuous(n.breaks = 10) +
  theme_bw()
scree

strsplit(eigenvecs$IID, "")
eigenvecs$pop <- paste0(sapply(strsplit(eigenvecs$IID, ""), `[`, 1), sapply(strsplit(eigenvecs$IID, ""), `[`, 2))

pc1pc2 <- ggplot() +
  geom_point(data = eigenvecs,
             aes(x = PC1, y = PC2, fill = pop),
             shape = 21, size = 2.5) +
  xlab(paste0("PC1 - ", round(percents[1], 2), "%")) +
  ylab(paste0("PC2 - ", round(percents[2], 2), "%")) +
  theme_bw()
pc1pc2

pc1pc3 <- ggplot() +
  geom_point(data = eigenvecs,
             aes(x = PC1, y = PC3, fill = pop),
             shape = 21, size = 2.5) +
  xlab(paste0("PC1 - ", round(percents[1], 2), "%")) +
  ylab(paste0("PC3 - ", round(percents[3], 2), "%")) +
  theme_bw()
pc1pc3

pc2pc3 <- ggplot() +
  geom_point(data = eigenvecs,
             aes(x = PC2, y = PC3, fill = pop),
             shape = 21, size = 2.5) +
  xlab(paste0("PC2 - ", round(percents[2], 2), "%")) +
  ylab(paste0("PC3 - ", round(percents[3], 2), "%")) +
  theme_bw()
pc2pc3

pc3pc4 <- ggplot() +
  geom_point(data = eigenvecs,
             aes(x = PC3, y = PC4, fill = pop),
             shape = 21, size = 2.5) +
  xlab(paste0("PC3 - ", round(percents[3], 2), "%")) +
  ylab(paste0("PC4 - ", round(percents[4], 2), "%")) +
  theme_bw()
pc3pc4

ggplot() +
  geom_point(data = eigenvecs,
             aes(x = PC8, y = PC9, fill = pop),
             shape = 21, size = 2.5) +
  xlab(paste0("PC6 - ", round(percents[3], 2), "%")) +
  ylab(paste0("PC7 - ", round(percents[4], 2), "%")) +
  theme_bw()


#### distance from raw ####
library(adegenet)
sample <- "ND364"
gt_data <- read.PLINK(paste0("data/vasa_quick/vasaship_quick.", sample,".raw"))
gt_data <- as.matrix(gt_data)
hc <- hclust(as.dist(1-cor(t(gt_data), use = "pairwise.complete.obs")), method = "complete")
plot(hc)
a <-  t(na.omit(t(as.matrix(gt_data))))
PCA <- rda(a, scale=T)
PCs <- data.frame(scores(PCA, choices=1:4, display="sites", scaling=0))

# load data table with coordinates and samples (manually selected for now)
full_table <- read.table("data/samples_table.csv", sep = ",",
                         header = TRUE, na.strings = "UNKNOWN") %>%
  column_to_rownames("sample_id")
data_table <- full_table[rownames(PCs), ]

PCs$region <- data_table$region
PCs$spawn <- data_table$spawn
sample_ids <- rownames(PCs)
ggplot()+
  geom_point(data = PCs, aes(x = PC1, y = PC2, fill = region), shape = 21, size = 3) +
  geom_text_repel(data = PCs, aes(x = PC1, y = PC2, label = sample_ids),
                  family = "Verdana",
                  size = 2,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  max.overlaps = Inf)

