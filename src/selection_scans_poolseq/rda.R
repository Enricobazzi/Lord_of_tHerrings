library(tidyverse)
library(vegan)
library(RColorBrewer)
# functions
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]              # locus names in these tails
}

# locations
locations <- c('A_Kalix_Baltic_Spring', 'B_Vaxholm_Baltic_Spring',
               'DalBoB_Atlantic_Autumn', 'DalFB_Atlantic_Spring',
               'DalGeB_Atlantic_Autumn', 'DalInB_Atlantic_Spring',
               'DalNsF_Atlantic_Autumn', 'DalNsS_Atlantic_Spring',
               'G_Gamleby_Baltic_Spring', 'HGS10_Downs_EnglishChannel_Winter',
               'HGS11_RingkobingFjord_NorthSea_Spring',
               'HGS12_BornholmBasin_Baltic_Autumn', 'HGS15_NSSH_Atlantic_Spring',
               'HGS16_Orkney_NorthSea_Autumn', 'HGS17_IsleOfMan_IrishSea_Autumn',
               'HGS18_CelticSea_Atlantic_AutumnWinter',
               'HGS19_TeelinBay_Atlantic_Winter', 'HGS1_Riga_Baltic_Spring',
               'HGS20_CapeWrath_Atlantic_Spring', 'HGS21_Hebrides_Atlantic_Mixed',
               'HGS22_CapeWrath_Atlantic_Autumn', 'HGS23_Clyde_Atlantic_Spring',
               'HGS24_Landvik_Atlantic_Spring', 'HGS25_Lindas_Atlantic_Spring',
               'HGS26_Lusterfjorden_Atlantic_Spring', 'HGS27_Gloppen_Atlantic_Spring',
               'HGS2_Riga_Baltic_Spring', 'HGS3_Riga_Baltic_Autumn',
               'HGS4_Riga_Baltic_Autumn', 'HGS5_Schlei_Baltic_Autumn',
               'HGS6_Schlei_Baltic_Spring', 'HGS71_Rugen_Baltic_Spring',
               'HGS72_Rugen_Baltic_Spring', 'HGS8_KattegatNorth_Atlantic_Spring',
               'HGS9_Greenland_Atlantic_Spring', 'H_Fehmarn_Baltic_Autumn',
               'J_Traslovslage_Baltic_Spring', 'LandvikS17_Norway_Baltic_Spring',
               'N_NorthSea_Atlantic_Autumn', 'O_Hamburgsund_Atlantic_Spring',
               'PB10_Skagerrak_Atlantic_Spring', 'PB11_Kalmar_Baltic_Spring',
               'PB12_Karlskrona_Baltic_Spring', 'PB1_HastKar_Baltic_Spring',
               'PB2_Iceland_Atlantic_Spring', 'PB4_Hudiksvall_Baltic_Spring',
               'PB5_Galve_Baltic_Spring', 'PB6_Galve_Baltic_Summer',
               'PB7_Galve_Baltic_Autumn', 'PB9_Kattegat_Atlantic_Spring',
               'PN3_CentralBaltic_Baltic_Spring', 'Q_Norway_Atlantic_Atlantic_Spring',
               'TysklandS18_Germany_Baltic')

pacific_pops <- c("HWS1_Japan_SeaOfJapan", "HWS2_PechoraSea_BarentsSea",
                  "HWS3_WhiteSea_WhiteSea", "HWS4_KandalakshaBay_WhiteSea",
                  "HWS5_KandalakshaBay_WhiteSea", "HWS6_Balsfjord_Atlantic",
                  "PB8_Pacific_Pacific_Spring")

# load env data - t and s
tem <- read.table("data/selection_scans_poolseq/sst_mean.txt", col.names = locations)
sal <- read.table("data/selection_scans_poolseq/sss_mean.txt", col.names = locations)
etable <- pivot_longer(tem,
                    locations,
                    names_to = "location",
                    values_to = "sst_mean") %>% column_to_rownames("location")
etable$sss_mean <- as.numeric(sal[1,])

# load genetic data - use 001 subset for now
freq <- read.table("data/selection_scans_poolseq/subsets_freqs/60.Neff.001.freq")[, -(1:2)]
colnames(freq) <- strsplit(readLines("data/published_data/60.Neff.freq", n = 1), "\t")[[1]][-c(1:2)]
freq <- freq[ , !(names(freq) %in% pacific_pops)]
freq <- t(na.omit(freq))

# run rda
rda <- rda(freq ~ ., data=etable, scale=TRUE)

# get candidates
load.rda <- scores(rda, choices=c(1, 2), display="species")
cand <- data.frame()
for(i in 1:2){
  candN <- outliers(load.rda[,i],2.5)
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
  snp.gen <- freq[,nam]
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

# prepare RDA plot
# samples portion:
sample_df <- data.frame(scores(rda, choices=c(1:2), display="sites", scaling="none"))
loc <- c()
for (locname in strsplit(rownames(sample_df), "_")){
  loc <- append(loc,
                paste(locname[length(locname) - 1], locname[length(locname)], sep = "_")
                )
}
sample_df$population <- loc
# add colors
col <- rep(NA, length(rownames(sample_df)))
col[grep("Baltic_Spring", rownames(sample_df))] <- "#E31A1C"
col[grep("Baltic_Autumn", rownames(sample_df))] <- "#FF7F00"
col[grep("Baltic_Summer", rownames(sample_df))] <- "#B8860b"
col[grep("Germany_Baltic", rownames(sample_df))] <- "#FDBF6F"
col[grep("Atlantic_Mixed", rownames(sample_df))] <- "#FB9A99"
col[grep("Atlantic_Autumn", rownames(sample_df))] <- "#3B528BFF"
col[grep("Atlantic_Winter", rownames(sample_df))] <- "#0F4909"
col[grep("IrishSea_Autumn", rownames(sample_df))] <- "#B2DF8A"
col[grep("Atlantic_Spring", rownames(sample_df))] <- "#21908CFF"
col[grep("EnglishChannel_Winter", rownames(sample_df))] <- "#440154FF"
col[grep("NorthSea_Spring", rownames(sample_df))] <- "#CAB2D6"
col[grep("NorthSea_Autumn", rownames(sample_df))] <- "#A035AF"
sample_df$color <- col

# snps portion:
snps_df <- data.frame(scores(rda, choices=c(1:2), display="species", scaling="none"))
snps_df$names <- row.names(snps_df)
snps_df$type <- "Neutral"
snps_df$type[snps_df$names %in% cand$snp] <- "Candidate"
col.pred <- c()
for (n in 1:nrow(snps_df)){
  row <- snps_df[n,]
  if (row$type == "Candidate"){
    if (cand[cand$snp == row$names, "predictor"] == "sst_mean"){
      col.pred <- append(col.pred, "yellow")
    } else {
      col.pred <- append(col.pred, "green")
    }
  } else {
    col.pred <- append(col.pred, "#f1eef6")
  }
}
snps_df$color <- col.pred
# divide SNPs into 2 data.frames
neutral_snps_df <- snps_df %>% filter(type == "Neutral")
candidate_snps_df <- snps_df %>% filter(type == "Candidate")

# variables data frame
TAB_var <- data.frame(scores(rda, choices=c(1:2), display="bp"))
# variance explained by rda axes
rda_ax_expl_constrain = round(x = (rda$CCA$eig / sum(rda$CCA$eig)) * 100,
                              digits = 2)

# plot!
library(ggrepel)
sample_labels <- sapply(strsplit(rownames(sample_df), "_"), `[`, 2)

sample_labels <- sapply(strsplit(rownames(sample_df), "_"), function(x) {
  paste(head(x, 1), collapse = "_")
})

x = 1
y = 2
rda_plot <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.3) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.3) +
  #geom_point(aes(x = neutral_snps_df[,x] * 3, y = neutral_snps_df[,y] * 3),
  #           shape = 4, cex = 1, fill = neutral_snps_df$color, alpha = 0.5) +
  geom_segment(aes(xend=TAB_var[,x]/2, yend=TAB_var[,y]/2, x=0, y=0),
               colour="black", size=0.15, linetype=1,
               arrow=arrow(length = unit(0.01, "npc"))) +
  geom_text(aes(x=1.05*(TAB_var[,x]/2), y=1.05*(TAB_var[,y]/2),
                label = as.vector(row.names(TAB_var))),
            size = 4, family = "Verdana") +
  #geom_point(aes(x = candidate_snps_df[,x] * 3, y = candidate_snps_df[,y] * 3),
  #           shape = 21, cex = 1.5, fill = candidate_snps_df$color, alpha = 0.8) +
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
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        plot.background = element_blank()) +
  # xlim(min = min_lim, max = max_lim) +
  # ylim(min = min_lim, max = max_lim) +
  xlab(paste0("RDA ", x, " - ", rda_ax_expl_constrain[x], "%")) +
  ylab(paste0("RDA ", y, " - ", rda_ax_expl_constrain[y], "%"))
rda_plot

# pca of candidate snps only
library(FactoMineR)
cand_freq <- freq[ , as.character(cand$snp)]
cand_pca <- FactoMineR::PCA(cand_freq, graph = F, scale.unit=F)
pcdf <- as.data.frame(cand_pca$ind$coord)
col <- rep(NA, length(rownames(pcdf)))
col[grep("Baltic_Spring", rownames(pcdf))] <- "#E31A1C"
col[grep("Baltic_Autumn", rownames(pcdf))] <- "#FF7F00"
col[grep("Baltic_Summer", rownames(pcdf))] <- "#B8860b"
col[grep("Germany_Baltic", rownames(pcdf))] <- "#FDBF6F"
col[grep("Atlantic_Mixed", rownames(pcdf))] <- "#FB9A99"
col[grep("Atlantic_Autumn", rownames(pcdf))] <- "#3B528BFF"
col[grep("Atlantic_Winter", rownames(pcdf))] <- "#0F4909"
col[grep("IrishSea_Autumn", rownames(pcdf))] <- "#B2DF8A"
col[grep("Atlantic_Spring", rownames(pcdf))] <- "#21908CFF"
col[grep("EnglishChannel_Winter", rownames(pcdf))] <- "#440154FF"
col[grep("NorthSea_Spring", rownames(pcdf))] <- "#CAB2D6"
col[grep("NorthSea_Autumn", rownames(pcdf))] <- "#A035AF"
pcdf$color <- col
ggplot() +
  geom_point(data = pcdf, aes(x=Dim.1, y=Dim.2), shape=21,
             fill=pcdf$color, size = 2.5)
