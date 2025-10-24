#### LOAD LIBRARIES ####

# install.packages("gradientForest", repos="http://R-Forge.R-project.org")
# install.packages("~/opt/RDA-forest-main/RDAforest_2.6.8.tar.gz")
library(RDAforest)
library(sdmpredictors)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(terra)
library(geosphere)
library(adegenet)
library(ggrepel)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### LISTS AND FUNCTIONS ####

# magma palette for map
map.pal <- rev(viridis(100, option = "magma"))

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### LOAD DATA ####

# load data table with coordinates and samples (manually selected for now)
data_table <- read.table("data/samples_table.csv", sep = ",", header = TRUE, na.strings = "UNKNOWN")
samples <- c("AAL1", "AAL2", "AAL3", "AF29", "AF30", "AF31", "AF8", "AK1", "AK2", "AK3", "AM27", "AM29", "AM33", "AM8", "BF16", "BF18", "BF19", "BF21", "BM14", "BM15", "BM16", "BM19", "Fehmarn3", "Fehmarn44", "Fehmarn6", "Gavle100", "Gavle54", "Gavle98", "MB1", "MB2", "MB3", "MHER001", "MHER003", "MHER006", "MHER007", "MHER008", "MHER009", "MHER010", "MHER011", "MHER012", "MHER013", "MHER014", "MHER015", "MHER016", "MHER017", "MHER018", "MHER019", "MHER020", "MHER021", "MHER022", "MHER023", "MHER024", "MHER025", "MHER026", "MHER027", "MHER028", "MHER029", "MHER030", "MHER034", "MHER035", "MHER036", "MHER037", "MHER038", "MHER039", "MHER044", "MHER045", "MHER046", "MHER052", "MHER053", "MHER054", "MHER055", "MHER056", "MHER057", "MHER058", "MHER059", "MHER061", "MHER062", "MHER063", "MHER065", "MHER066", "NorthSea13", "NorthSea19", "NorthSea34", "NSSH33", "NSSH34", "NSSH36", "Z12", "Z14", "Z4")
data_table <- data_table[data_table$sample_id %in% samples, ]

gt_data <- read.PLINK("data/rdaforest/ccc.raw")
gt_data <- as.matrix(gt_data)


PCA <- rda(gt_data, scale=T)
PCs <- data.frame(scores(PCA, choices=1:4, display="sites", scaling=0))
PCs$region <- data_table$region
PCs$sample <- data_table$sample_id
PCs$spawn <- data_table$spawn

ggplot()+
  geom_point(data = PCs, aes(x = PC1, y = PC2, fill = region), shape = 21)

