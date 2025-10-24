library(tidyverse)

# locations
locations <- c('CHROM', 'POS', 'A_Kalix_Baltic_Spring', 'B_Vaxholm_Baltic_Spring',
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

cor_test_wrapper <- function(p_vec, env_vector){
  correlation_result <- cor.test(p_vec, env_vector, method = "kendall", exact = F)
  return(c(correlation_result$estimate,
           correlation_result$p.val))
}

# load genetic data - use 001 subset for now
freq <- read.table("data/selection_scans_poolseq/subsets_freqs/60.Neff.001.freq")#[, -(1:2)]
colnames(freq) <- strsplit(readLines("data/published_data/60.Neff.freq", n = 1), "\t")[[1]]#[-c(1:2)]
freq <- freq[ , (names(freq) %in% locations)]
freq <- freq[rowSums(is.na(freq)) <= 0.2 * 53, ]

# env data
tem <- read.table("data/selection_scans_poolseq/sst_mean.txt", col.names = locations[-c(1:2)])
sal <- read.table("data/selection_scans_poolseq/sss_mean.txt", col.names = locations[-c(1:2)])
env <- pivot_longer(tem,
                    locations[-c(1:2)],
                    values_to = "sst_mean") %>% column_to_rownames("name")
env$sss_mean <- as.numeric(sal[1,])


# first row
aaa <- data.frame()
for (i in 1:nrow(freq)){
  print(i)
  a <- cor_test_wrapper(t(freq[i,-(1:2)]), env$sst_mean)
  aaa <- rbind(aaa, a)
}
zzz <- cbind(freq$CHROM, freq$POS, aaa)
colnames(zzz) <- c("CHROM", "POS", "kend", "pval")
for (i in 1:26){
  chr <- zzz[zzz$CHROM == paste0("chr", i), ]
  plot(chr$POS, -log10(chr$pval))
  title(paste0("chr", i))
}

zzz <- drop_na(zzz)
zzz[(zzz$CHROM == "chr12") & (-log10(zzz$pval) > 8), ]

candi <- zzz[(-log10(zzz$pval) > 5), ]
nnn <- freq %>%
  semi_join(candi, by = c("CHROM", "POS"))
