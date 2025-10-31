library(tidyverse)
library(vegan)
library(RColorBrewer)
library(adegenet)
library(ggrepel)

# functions
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]              # locus names in these tails
}

samples <- c("AAL1", "AAL2", "AAL3", "AF29", "AF30", "AF31", "AF8", "AK1", "AK2", "AK3", "AM27", "AM29", "AM33", "AM8", "BF16", "BF18", "BF19", "BF21", "BM14", "BM15", "BM16", "BM19", "Fehmarn3", "Fehmarn44", "Fehmarn6", "Gavle100", "Gavle54", "Gavle98", "MHER001", "MHER003", "MHER006", "MHER007", "MHER008", "MHER009", "MHER010", "MHER011", "MHER012", "MHER013", "MHER014", "MHER015", "MHER016", "MHER017", "MHER018", "MHER019", "MHER020", "MHER021", "MHER022", "MHER023", "MHER024", "MHER025", "MHER026", "MHER027", "MHER028", "MHER029", "MHER030", "MHER034", "MHER035", "MHER036", "MHER037", "MHER038", "MHER039", "MHER044", "MHER045", "MHER046", "MHER052", "MHER053", "MHER054", "MHER055", "MHER056", "MHER057", "MHER058", "MHER059", "MHER061", "MHER062", "MHER063", "MHER065", "MHER066", "NorthSea13", "NorthSea19", "NorthSea34", "NSSH33", "NSSH34", "NSSH36", "Z12", "Z14", "Z4")

# data table
data_table <- read.table("data/samples_table.csv", sep = ",", header = TRUE, na.strings = "UNKNOWN")
data_table <- data_table[data_table$sample_id %in% samples, ]


# env: environmental values for samples
etable <- read.table("data/rdaforest/etable.csv", sep=",")[, c(2,3)]
# eventually subset but for now only useful samples are in etable.csv
rownames(etable) <- samples

# geno: genetic data
gt_data <- read.PLINK("data/rdaforest/raws/modern_samples.step3.biall_acanaf.02missmax_0fmiss_005maf.raw")
gt_data <- as.matrix(gt_data)
gt_data <- gt_data[rownames(gt_data) %in% samples, ]
geno <- gt_data
# geno <- gt_data[, cand$snp]

# run RDA
rda <- rda(geno ~ ., data=etable, scale=TRUE)

# get candidates
load.rda <- scores(rda, choices=c(1, 2), display="species")
cand <- data.frame()
for(i in 1:2){
  candN <- outliers(load.rda[,i],3)
  candN <- cbind.data.frame(rep(i,times=length(candN)), names(candN), unname(candN))
  colnames(candN) <- c("axis","snp","loading")
  cand <- rbind(cand, candN)
}
cand$snp <- as.character(cand$snp)
cand <- cand[!duplicated(cand$snp),]
ncand <- NROW(cand)

# assign vars to candidates by correlation
foo <- matrix(nrow=(ncand), ncol=NCOL(etable))
colnames(foo) <- colnames(etable)
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- geno[,nam]
  foo[i,] <- apply(etable, 2, function(x) cor(x, snp.gen))
}
cand <- cbind.data.frame(cand,foo)  
# assign SNPs to predictor based on highest correlation value - might be worth exploring
predictors <- c()
correlations <- c()
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  preds <- c(abs(bar["sst_mean"]), abs(bar["sss_mean"]))
  predictor <- as.character(names(which.max(preds)))
  correlation <- as.numeric(preds[which.max(preds)])
  predictors <- append(predictors, predictor)
  correlations <- append(correlations, correlation)
}
cand$predictor <- predictors
cand$correlation <- correlations

# save candidate snp table
write.table(cand, "data/rdaforest/cand_snps.csv", sep = ",", col.names = TRUE, row.names = FALSE)

# prepare RDA plot
# samples portion:
sample_df <- data.frame(scores(rda, choices=c(1:2), display="sites", scaling="none"))
sample_df$population <- data_table$region
# colors
col <- rep(NA, length(rownames(sample_df)))
col[grep("BALTIC", sample_df$population)] <- "#E31A1C"
col[grep("NE-ATL", sample_df$population)] <- "#3B528BFF"
col[grep("TRANS", sample_df$population)] <- "#A035AF"
sample_df$color <- col
# labels
sample_labels <- rownames(sample_df)
# variance explained by rda axes
TAB_var <- data.frame(scores(rda, choices=c(1:2), display="bp"))
rda_ax_expl_constrain = round(x = (rda$CCA$eig / sum(rda$CCA$eig)) * 100,
                              digits = 2)
# plot!
x = 1
y = 2
rda_plot <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.3) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.3) +
  geom_segment(aes(xend=TAB_var[,x]/2, yend=TAB_var[,y]/2, x=0, y=0),
               colour="black", size=0.15, linetype=1,
               arrow=ggplot2::arrow(length = unit(0.01, "npc"))) +
  geom_text(aes(x=1.05*(TAB_var[,x]/2), y=1.05*(TAB_var[,y]/2),
                label = as.vector(row.names(TAB_var))),
            size = 3, family = "Verdana") +
  geom_point(aes(x = sample_df[,x], y = sample_df[,y]),
             shape = 21, cex = 3, fill = sample_df$color, alpha = 0.8) +
  geom_text_repel(aes(x = sample_df[, x],
                      y = sample_df[, y],
                      label = sample_labels),
                  family = "Verdana",
                  size = 2,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  max.overlaps = Inf) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme_bw(base_size = 14, base_family = "Verdana") +
  #theme(panel.background = element_blank(), panel.grid = element_blank(),
  #      plot.background = element_blank()) +
  theme(panel.grid = element_blank()) +
  xlab(paste0("RDA ", x, " - ", rda_ax_expl_constrain[x], "%")) +
  ylab(paste0("RDA ", y, " - ", rda_ax_expl_constrain[y], "%"))
rda_plot
ggsave("plots/rdaforest/rda/modern_samples.allsnps.png", plot = rda_plot, height =  10, width = 10)

####################################################################################
# "simple" pca of rda outlier snps
PCA <- rda(geno[, cand$snp], scale=T)
sample_df <- data.frame(scores(PCA, choices=c(1:4), display="sites", scaling="none"))
sample_df$population <- data_table$region
# colors
col <- rep(NA, length(rownames(sample_df)))
col[grep("BALTIC", sample_df$population)] <- "#E31A1C"
col[grep("NE-ATL", sample_df$population)] <- "#3B528BFF"
col[grep("TRANS", sample_df$population)] <- "#A035AF"
sample_df$color <- col
# labels
sample_labels <- rownames(sample_df)
# variance explained by rda axes
rda_ax_expl_constrain = round(x = (PCA$CA$eig / sum(PCA$CA$eig)) * 100,
                              digits = 2)
# plot!
x = 1
y = 2
pca_plot <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.3) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.3) +
  geom_point(aes(x = sample_df[,x], y = sample_df[,y]),
             shape = 21, cex = 3, fill = sample_df$color, alpha = 0.8) +
  geom_text_repel(aes(x = sample_df[, x],
                      y = sample_df[, y],
                      label = sample_labels),
                  family = "Verdana",
                  size = 2,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  max.overlaps = Inf) +
  theme_bw(base_size = 14, base_family = "Verdana") +
  theme(panel.grid = element_blank()) +
  xlab(paste0(colnames(sample_df)[x], " - ", rda_ax_expl_constrain[x], "%")) +
  ylab(paste0(colnames(sample_df)[y], " - ", rda_ax_expl_constrain[y], "%"))
pca_plot
