library(tidyverse)
library(adegenet)
library(vegan)
library(ggrepel)

cor_test_wrapper <- function(p_vec, env_vector){
  correlation_result <- cor.test(p_vec, env_vector, method = "kendall", exact = F)
  return(c(correlation_result$estimate,
           correlation_result$p.val))
}

samples <- c("AAL1", "AAL2", "AAL3", "AF29", "AF30", "AF31", "AF8", "AK1", "AK2", "AK3", "AM27", "AM29", "AM33", "AM8", "BF16", "BF18", "BF19", "BF21", "BM14", "BM15", "BM16", "BM19", "Fehmarn3", "Fehmarn44", "Fehmarn6", "Gavle100", "Gavle54", "Gavle98", "MHER001", "MHER003", "MHER006", "MHER007", "MHER008", "MHER009", "MHER010", "MHER011", "MHER012", "MHER013", "MHER014", "MHER015", "MHER016", "MHER017", "MHER018", "MHER019", "MHER020", "MHER021", "MHER022", "MHER023", "MHER024", "MHER025", "MHER026", "MHER027", "MHER028", "MHER029", "MHER030", "MHER034", "MHER035", "MHER036", "MHER037", "MHER038", "MHER039", "MHER044", "MHER045", "MHER046", "MHER052", "MHER053", "MHER054", "MHER055", "MHER056", "MHER057", "MHER058", "MHER059", "MHER061", "MHER062", "MHER063", "MHER065", "MHER066", "NorthSea13", "NorthSea19", "NorthSea34", "NSSH33", "NSSH34", "NSSH36", "Z12", "Z14", "Z4")

# data table
data_table <- read.table("data/samples_table.csv", sep = ",", header = TRUE, na.strings = "UNKNOWN")
data_table <- data_table[data_table$sample_id %in% samples, ]

# load genetic data
gt_data <- read.PLINK("data/rdaforest/raws/modern_samples.step3.biall_acanaf.02missmax_0fmiss_005maf.raw")
gt_data <- as.matrix(gt_data)
gt_data <- gt_data[rownames(gt_data) %in% samples, ]

# env: environmental values for samples
env <- read.table("data/rdaforest/etable.csv", sep=",")[, c(2,3)]
# eventually subset but for now only useful samples are in etable.csv
rownames(env) <- samples

# for each column (snp)
aaa <- data.frame()
for (i in 1:ncol(gt_data)){
  print(i)
  a <- cor_test_wrapper(gt_data[, i], env$sst_mean)
  aaa <- rbind(aaa, a)
}

chrom <- sapply(strsplit(colnames(gt_data), ":"), `[`, 1)
pos <- sapply(strsplit(colnames(gt_data), ":"), function(y) as.numeric(strsplit(y[2], "_")[[1]][1]))

zzz <- data.frame(
  CHROM = chrom,
  POS = pos,
  SNP = colnames(gt_data),
  kend = aaa[, 1],
  pval = aaa[, 2]
)

for (chr in unique(chrom)){
  if (grepl("^CM0", chr)){
    chr_df <- zzz[zzz$CHROM == chr, ]
    plot(chr_df$POS, -log10(chr_df$pval))
    title(paste0(chr))
  }
}

candtau <- zzz %>%
  filter(
    pval < quantile(zzz$pval, 0.01, na.rm = TRUE)
  )

# TRY SUBSETTING TO INVERSION
gt_df <- as.data.frame(t(gt_data))
gt_df$CHR <- as.character(chrom)
gt_df$POS <- as.numeric(pos)
bb <- gt_df %>%
  filter(
    CHR == "CM079356.1",
    POS > 22500000,
    POS < 25200000
  )
bb <- gt_df[candtau$SNP, ]
plot(hclust(dist(t(bb[, samples]))))

PCA <- rda(t(bb[, samples]), scale=T)
sample_df <- data.frame(scores(PCA, choices=c(1:4), display="sites", scaling="none"))
sample_df$population <- column_to_rownames(data_table  %>% `rownames<-`( NULL ), var = "sample_id")[rownames(sample_df), "region"]
col <- rep(NA, length(rownames(sample_df)))
col[grep("BALTIC", sample_df$population)] <- "#E31A1C"
col[grep("NE-ATL", sample_df$population)] <- "#3B528BFF"
col[grep("TRANS", sample_df$population)] <- "#A035AF"
sample_df$color <- col
sample_labels <- rownames(sample_df)
rda_ax_expl_constrain = round(x = (PCA$CA$eig / sum(PCA$CA$eig)) * 100,
                              digits = 2)
sample_df$sst <- env$sst_mean
sample_df$sss <- env$sss_mean

scale <- 3 + (log(sample_df$sst) - min(log(sample_df$sst))) * (8 - 3) / (max(log(sample_df$sst)) - min(log(sample_df$sst)))
scale <- 3 + (sample_df$sst - min(sample_df$sst)) * (8 - 3) / (max(sample_df$sst) - min(sample_df$sst))

# plot!
x = 1
y = 2
pca_plot <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.3) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.3) +
  geom_point(aes(x = sample_df[,x], y = sample_df[,y], fill = sample_df[, "sss"]),
             shape = 21, cex = 3, alpha = 0.8) +
  geom_text_repel(aes(x = sample_df[, x],
                      y = sample_df[, y],
                      label = sample_labels),
                  family = "Verdana",
                  size = 2,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  max.overlaps = Inf) +
  theme_bw(base_size = 14, base_family = "Verdana") +
  scale_fill_viridis() +
  theme(panel.grid = element_blank()) +
  xlab(paste0(colnames(sample_df)[x], " - ", rda_ax_expl_constrain[x], "%")) +
  ylab(paste0(colnames(sample_df)[y], " - ", rda_ax_expl_constrain[y], "%"))
pca_plot
