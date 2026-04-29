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
dataset <- datasets[4]
for (dataset in datasets){
  aa <- read.table(paste0("data/diversity/output/", dataset, ".folded.thetas.W10000.S5000.pestPG"),
                   comment.char = "", header = T)
  print(dataset)
  # print(mean(aa$tP/aa$nSites, na.rm =T))
  print(exp(mean(log(aa$tP/aa$nSites), na.rm = T)))
}
