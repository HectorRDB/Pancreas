suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(tidyverse)
    library(here)
})


source(here("Scripts", "Replicability", "variable_genes.R"))
source(here("Scripts", "Replicability", "datasets.R"))


create_data = function() {
    dataset = list(
        baron = readRDS("/scratch/users/singlecell/Pancreas/RawData/baron.rds"),
        # muraro = readRDS("/scratch/users/singlecell/Pancreas/RawData/muraro.rds"),
        # xin = readRDS("/scratch/users/singlecell/Pancreas/RawData/xin.rds"),
        segerstolpe = readRDS("/scratch/users/singlecell/Pancreas/RawData/segerstolpe.rds")
    )

    col_data <- fuse_coldata(dataset)
    dataset <- fuse_datasets(dataset)
    dataset <- SingleCellExperiment(
      assays = list("counts" = as.matrix(dataset)),
      colData = as(col_data, "DataFrame")
      )
    hvg <- variable_genes(dataset)
    rowData(dataset)$is_hvg <- rownames(dataset) %in% hvg

    saveRDS(dataset, here("Data", "full_data.rds"))
}