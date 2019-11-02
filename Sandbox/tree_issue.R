library(SummarizedExperiment)
library(parallel)
library(matrixStats)
library(clusterExperiment)
library(tidyverse)
library(Dune)
library(mclust)
loc <- "/accounts/projects/epurdom/singlecell/Pancreas/Data/singleMethod/segerstolpe"
rsec <- "/scratch/users/singlecell/Pancreas/ProcessedData/segerstolpe_Rsec.rds"

sc3 <- read.csv(paste0(loc, "_SC3.csv"))[, -1]
colnames(sc3) <- str_remove(colnames(sc3), "^X") %>% str_replace("\\.", "-")
sc3 <- sc3[, "0"] %>% as.numeric()
Rsec <- readRDS(rsec)

Rsec <- addClusterings(Rsec, sc3, clusterLabels = "sc3")
Rsec <- makeDendrogram(Rsec, whichCluster = "sc3")