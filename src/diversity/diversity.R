library(tidyverse)
library(ggrepel)
library(adegenet)
library(cowplot)

snps_bed <- "data/angsd_matrix/sites/supplementary_file_7.v2.bed"
snps <- read.table(snps_bed, col.names = c("Chr", "St", "En"))

chromosomes <- c("NC_045152.1", "NC_045153.1", "NC_045154.1", "NC_045155.1", 
                 "NC_045156.1", "NC_045157.1", "NC_045158.1", "NC_045159.1",
                 "NC_045160.1", "NC_045161.1", "NC_045162.1", "NC_045163.1",
                 "NC_045164.1", "NC_045165.1", "NC_045166.1", "NC_045167.1",
                 "NC_045168.1", "NC_045169.1", "NC_045170.1", "NC_045171.1",
                 "NC_045172.1", "NC_045173.1", "NC_045174.1", "NC_045175.1",
                 "NC_045176.1", "NC_045177.1")

datasets <- c("sk-17sp", "sk-18rh", "sk-18sp", "sk-mh", "ns-ah", "ns-18rh", "ns-18sp", "ns-mh")

read_thetas_chr <- function(dataset, chr) {
  aa <- read.table(paste0("data/diversity/output/", dataset, ".folded.thetas.W10000.S5000.pestPG"),
                   comment.char = "", header = T)
  aaa <- aa |> 
    filter(Chr == chr) |>
    mutate(Pi = tP/nSites, tWat = tW/nSites, TD = Tajima, DataSet = dataset) |>
    select(DataSet, Chr, WinCenter, Pi, tWat, TD)
  return (aaa)
}

baseplot <- function() {
  bp <- ggplot()
  return (bp)
}

plot_vari <- function(thetas, vari) {
  p <- geom_point(
    data = thetas,
    aes(x = WinCenter, y = .data[[vari]], colour = DataSet),
    alpha = 0.3, size = 0.1
    )
  return(p)
}

plot_vari_smooth <- function(thetas, vari) {
  p <- geom_smooth(
    data = thetas,
    aes(x = WinCenter, y = .data[[vari]], colour = DataSet),
    se = FALSE, method = "gam", formula = y ~ s(x)
  )
  return(p)
}

for (sea in c("sk", "ns")) {
  
  for (chr in chromosomes) {
    
    for (vari in c("Pi", "tWat", "TD")) {
      
      plot <- baseplot()
      for (dataset in datasets){
        if (startsWith(dataset, sea)){
          thetas <- read_thetas_chr(dataset, chr)
          plot <- plot + plot_vari(thetas, vari)
        }
      }
      for (dataset in datasets){
        if (startsWith(dataset, sea)){
          thetas <- read_thetas_chr(dataset, chr)
          plot <- plot + plot_vari_smooth(thetas, vari)
        }
      }
      
      if (vari == "Pi" || vari == "tWat"){
        plot <- plot + geom_vline(data = snps |> filter(Chr == chr) , aes(xintercept = St), linewidth = 0.1) + 
          ylim(c(0,0.008)) + theme_minimal()
      } else {
        plot <- plot + geom_vline(data = snps |> filter(Chr == chr), aes(xintercept = St), linewidth = 0.1) + 
          theme_minimal()
      }
      
      ggsave(
        filename = paste0("plots/diversity/", sea, "/", vari, ".", chr, ".pdf"),
        plot = plot, height = 4, width = 8
      )
    }
  }
}

####

datasets <- c("sk-17sp", "sk-18rh", "sk-18sp", "sk-mh", "ns-ah", "ns-18rh", "ns-18sp", "ns-mh")
datasets <- c("sk-17sp", "sk-18rh", "sk-18sp", "sk-mh")

for (dataset in datasets){
  aa <- read.table(paste0("data/diversity/output/", dataset, ".folded.thetas.W10000.S5000.pestPG"),
                   comment.char = "", header = T)
  print(dataset)
  # print(mean(aa$tP/aa$nSites, na.rm =T))
  # print(exp(mean(log(aa$tP/aa$nSites), na.rm = T)))
  print(exp(weighted.mean(log(aa$tP / aa$nSites), w = aa$nSites, na.rm = TRUE)))
  # print(sum(aa$tP, na.rm = TRUE) / sum(aa$nSites, na.rm = TRUE))
  print(exp(weighted.mean(log(aa$tW / aa$nSites), w = aa$nSites, na.rm = TRUE)))
  # print(sum(aa$tW, na.rm = TRUE) / sum(aa$nSites, na.rm = TRUE))
  # print(weighted.mean(aa$Tajima, w = aa$nSites, na.rm = TRUE))
}

####

library(tidyverse)

pcangsd_dataset <- "full_herr"
dataset <- "wp1_final_bal"

# get sample names for a dataset from the sample_list file:
get_samples_from_dataset <- function(dataset) {
  file_path <- paste0("data/angsd_matrix/bamlists/", dataset, ".sample_list.txt")
  samples <- read.table(file_path)[, 1] |> as.character()
  return (samples)
}

# get metadata table
get_metadata <- function(dataset) {
  samples <- get_samples_from_dataset(dataset)
  sample_data_file <- "~/Documents/Silly-periods/data/samples_table.csv"
  sample_data <- read.table(sample_data_file, sep = ",",
                            header = TRUE, na.strings = "UNKNOWN")
  sample_data <- sample_data[sample_data$sample_id %in% samples, ]
  return (sample_data)
}

# get sample het
get_het <- function(sample) {
  # from https://www.popgen.dk/angsd/index.php/Heterozygosity
  a <- scan(paste0("data/diversity/output/", sample, ".het_notrans.ml"))
  return (a[2] / sum(a))
}

# get stock dataframe
get_stock_df <- function(pcangsd_dataset, dataset) {
  return (
    read.table(paste0("data/angsd_matrix/dapc/stocks_table.",
                      pcangsd_dataset, ".", dataset, ".csv"),
               sep = ",", header = T)
    )
}

# get sample stock 
get_sample_stock <- function(sample, stock_df){
  stock <- stock_df$Stock[which(stock_df$sample == sample)]
  return (stock)
}

samples <- get_samples_from_dataset(dataset)
metadata <- get_metadata(dataset)

het_df <- data.frame(
  ID = metadata$new.id,
  Sample = metadata$sample_id,
  Region = metadata$region,
  Period = metadata$period,
  Year = metadata$year,
  Stock = sapply(metadata$sample_id, get_sample_stock,
                 stock_df = get_stock_df(pcangsd_dataset, dataset),
                 USE.NAMES = FALSE),
  Heterozygosity = unlist(lapply(metadata$sample_id, get_het))
)
stock <- "Skagerrak & Kattegat Spring-Spawner"
stock <- "Norwegian Spring-Spawner"
stock <- "North Sea Autumn-Spawner"
stock <- "Ireland & Britain Autumn-Spawner"
z <- het_df |> filter(Stock == stock)
for (period in c("17sp", "18rh", "18sp", "mh")){
  print(paste(stock, period))
  zz <- z |> filter(Period == period)
  print(mean(zz$Heterozygosity))
}

region <- "North_Sea"
z <- het_df |> filter(Region == region)
for (period in c("mh")){
  print(paste(region, period))
  zz <- z |> filter(Period == period)
  print(mean(zz$Heterozygosity))
}
