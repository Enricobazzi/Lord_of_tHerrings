# from https://github.com/aalbrechtsen/NGSadmix/blob/master/TUTORIAL.md
# pop <- read.table("Demo/Data/Demo1pop.info", as.is = TRUE)[, 1]
# q <- read.table("Demo/Results/Demo1NGSadmix.qopt")
# ord <- order(pop)
# par(mar = c(7, 4, 1, 1))
# barplot(
#   t(q)[, ord],
#   col = c(2, 1, 3),
#   names.arg = pop[ord],
#   las = 2,
#   ylab = "Demo1 admixture proportions",
#   cex.names = 0.75
# )

# get sample names for a dataset from the sample_list file:
get_samples_from_dataset <- function(dataset) {
  file_path <- paste0("data/angsd_matrix/bamlists/", dataset, ".sample_list.txt")
  samples <- read.table(file_path)[, 1] |> as.character()
  return (samples)
}

# get metadata table
get_newids_from_dataset <- function(dataset) {
  samples <- get_samples_from_dataset(dataset)
  sample_data_file <- "~/Documents/Silly-periods/data/samples_table.csv"
  sample_data <- read.table(sample_data_file, sep = ",",
                            header = TRUE, na.strings = "UNKNOWN")
  sample_data <- sample_data[sample_data$sample_id %in% samples, ]
  return (sample_data$new.id)
}

dataset <- "wp1_final_bal"
samples <- get_newids_from_dataset(dataset)
q <- read.table(paste0("data/ngsadmix/wp1_final_bal.admix.k_4.seed_1.qopt"), header=FALSE)

pdf("admix.pdf", width = 60)
barplot(
  t(q),
  col = c(2:20),
  names.arg = samples,
  las = 2,
  ylab = "admixt",
  cex.names = 0.75
)
dev.off()
