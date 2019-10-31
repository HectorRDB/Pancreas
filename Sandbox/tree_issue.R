library(SummarizedExperiment)
library(parallel)
library(matrixStats)
library(clusterExperiment)
library(tidyverse)
library(Dune)
library(mclust)
loc <- "/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/segerstolpe"
rsec <- "/scratch/users/singlecell/Pancreas/ProcessedData/segerstolpe_Rsec.rds"
out <- "/accounts/projects/epurdom/singlecell/Pancreas/Data/Dune/segerstolpe_comp2"
sc3_p <- "0"
seurat_p <- "2.1.30"
monocle_p <- "k_15"

sc3 <- read.csv(paste0(loc, "_SC3.csv"))[, -1]
colnames(sc3) <- str_remove(colnames(sc3), "^X") %>% str_replace("\\.", "-")
Names <- sc3$cells
sc3 <- sc3[, sc3_p] %>% as.numeric()
Monocle <- read.csv(paste0(loc, "_Monocle.csv"))[, -1]
Monocle <- as.data.frame(Monocle)[, monocle_p] %>% as.numeric()
Rsec <- readRDS(rsec)

for (clustering in c("sc3", "Monocle", "Seurat")) {
  Rsec <- addClusterings(Rsec, get(clustering), clusterLabels = clustering)
}

cutoffs <- seq(from = 0, to = 1, by = .05)
res <- list()
clustering <- "sc3"
print(clustering)
Rsec2 <- makeDendrogram(Rsec, whichCluster = clustering)