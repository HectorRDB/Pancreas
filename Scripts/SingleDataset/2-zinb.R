suppressWarnings(library(optparse))

# Arguments for R Script ----
option_list <- list(
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Where to store the reduced dim object"
  ),
  make_option(c("-l", "--location"),
              action = "store", default = NA, type = "character",
              help = "The location of the data"
  ),
  make_option(c("-n", "--nCores"),
              action = "store", default = 1,
              help = "Number of cores to use [default %default]"
  ),
  make_option(c("-i", "--initial"),
              action = "store", default = 10,
              help = "The smallest value of K [default %default]"
  ),
  make_option(c("-f", "--final"),
              action = "store", default = 50,
              help = "The largest value of K [default %default]"
  ),
  make_option(c("-d", "--dims"),
              action = "store", default = 5,
              help = "How many values of K between i and f [default %default]"
  ),
  make_option(c("-p", "--plots"),
              action = "store", default = NA, type = "character",
              help = "Location of the visual output. Default to [default %default], no output"
  ),
  make_option(c("-m", "--meta"),
              action = "store", default = NA, type = "character",
              help = "Location of metadata"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$l)) {
  loc <- opt$l
  cat("The selected dataset is located at", loc, "\n")
} else {
  stop("Missing l argument\n")
}
if (!is.na(opt$o)) {
  output_r <- opt$o
} else {
  stop("Missing o argument\n")
}

library(Rtsne)
library(stringr)
library(clusterExperiment)
library(dplyr)
library(BiocParallel)
library(zinbwave)
library(matrixStats)
library(ggplot2)

# Load data ----
sce <- readRDS(file = loc)

# Run ZinbWave ----
NCORES <- as.numeric(opt$n)
BiocParallel::register(MulticoreParam(NCORES))

vars <- matrixStats::rowVars(logcounts(sce))
cat("Running with K = 0 on the full data\n")
cat("Number of cores:", NCORES, "\n")
cat("Time to run zinbwave (seconds):\n")
if (ncol(sce) < 5000) {
  print(system.time(zinb0 <- zinbwave(sce, X = "~human")))  
}

ind <- vars > sort(vars,decreasing = TRUE)[1000]
whichGenes <- rownames(sce)[ind]
zinbDims <- floor(seq(from = opt$i, to = opt$f, length.out = opt$d))
cat("Using the following values for K :", zinbDims)
sceVar <- sce[ind,]

zinbWs <- lapply(zinbDims, function(zinbDim) {
  cat("Running with K = ", zinbDim, " on the filtered data\n")
  cat("Number of cores:", NCORES, "\n")
  cat("Time to run zinbwave (seconds):\n")
  print(system.time(zinb <- zinbwave(sceVar, K = zinbDim, X = "~human")))
  return(zinb)
})

ref_clusters <- read.csv(opt$m)$cell_type1

for (i in 1:length(zinbWs)) {
  type <- paste0("zinb-K", zinbDims[i])
  print(type)
  print("....Saving data")
  reducedDim(sce, type = type) <- zinbW <- reducedDim(zinbWs[[i]])
  if (!is.na(opt$p)) {
    print("....t-SNE")
    TNSE <- Rtsne(zinbW, initial_dims = min(50, zinbDims[i]))
    df <- data.frame(x = TNSE$Y[, 1], y = TNSE$Y[, 2],
                     col = ref_clusters)
    p <- ggplot(df, aes(x = x, y = y, col = col)) +
      geom_point(size = .4, alpha = .3) +
      theme_classic() +
      labs(x = "dim1", y = "dim2")
    ggsave(paste0(opt$p, "_K_", zinbDims[i], ".png"), p)
    print("....Saving plot")
  }
}

saveRDS(sce, file = output_r)
