suppressWarnings(library(optparse))

# Arguments for R Script ----
option_list <- list(
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Where to store the output"
  ),
  make_option(c("-l", "--location"),
              action = "store", default = NA, type = "character",
              help = "The location of the data"
  ),
  make_option(c("-p", "--preproc"),
              action = "store", default = NA, type = "character",
              help = "Where to store the normalized data"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$l)) {
  loc <- opt$l
  cat("The selected dataset is located at", loc)
} else {
  stop("Missing l argument")
}

if (!is.na(opt$o)) {
  output <- opt$o
} else {
  stop("Missing o argument")
}

# Load data ----
library(Seurat)
library(dplyr)
library(stringr)
library(SingleCellExperiment)

sce <- readRDS(file = loc)

# Setup ----
sSeurat <- CreateSeuratObject(counts = assays(sce)$counts, project = 'allen40K')
sSeurat <- NormalizeData(object = sSeurat, normalization.method = "LogNormalize")
sSeurat <- FindVariableFeatures(object = sSeurat, mean.function = ExpMean,
                                dispersion.function = LogVMR, do.plot = F)
sSeurat <- ScaleData(object = sSeurat, vars.to.regress = c("nCount_RNA", "human"))
sce <- as.SingleCellExperiment(sSeurat)
sSeurat <- RunPCA(object = sSeurat, ndims.print = 1, npcs = 100)

# Run clustering ----
clusterMatrix <- NULL
for (RESOLUTION in seq(from = 0.3, to = 2.5, by = .1)) {
  print(RESOLUTION)
  sSeurat_star <- FindNeighbors(sSeurat, dims = 1:50)
  sSeurat_star <- FindClusters(sSeurat_star, resolution = RESOLUTION)
  cluster_labels <- Idents(sSeurat_star)
  clusters <- unique(Idents(sSeurat_star))
  for (cluster in clusters) {
    print(paste0("...", cluster))
    sSeurat_loc <- sSeurat[, Idents(sSeurat_star) == cluster]
    sSeurat_loc <- FindNeighbors(sSeurat_loc, dims = 1:50)
    sSeurat_loc <- FindClusters(sSeurat_loc, resolution = RESOLUTION)
    cluster_labels[Idents(sSeurat_star) == cluster] <- 
      paste0(cluster, "_", Idents(sSeurat_loc))
  }
  clusterMatrix <- cbind(clusterMatrix, cluster_labels)
  colnames(clusterMatrix)[ncol(clusterMatrix)] <- RESOLUTION
}

clusterMatrix <- as.data.frame(clusterMatrix)
clusterMatrix$cells <- colnames(sce)
write.csv(clusterMatrix, file = output)
