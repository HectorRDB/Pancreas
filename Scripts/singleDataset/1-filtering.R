library(optparse)

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
  make_option(c("-c", "--cutoff"),
              action = "store", default = 50, type = "character",
              help = "The cutoff for filtering"
  ),
  make_option(c("-r", "--reads"),
              action = "store", default = 50, type = "character",
              help = "The read cutoff for filtering"
  ),
  make_option(c("-m", "--meta"),
              action = "store", default = NA, type = "character",
              help = "Where to store the metadata"
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

library(dplyr)
library(stringr)
library(SingleCellExperiment)

# Load data per se ----

Sce <- readRDS(loc)
Sce$human <- word(colnames(Sce), 1, sep = "_") %>% unlist()
if (str_detect(loc, "segerstolpe")) {
  Sce <- Sce[, Sce$cell_quality == "OK"]
  Sce <- Sce[, Sce$human != "AZ"]
}
write.csv(colData(Sce), opt$m)
cat("Preparing the data", "\n")
counts <- counts(Sce)
counts <- as.matrix(counts)
print(dim(counts))
filt <- rowSums(counts >= opt$r) >= opt$c
print(sum(filt))
print(sum(counts[filt, ]) / sum(counts))
Sce <- Sce[filt, ]

cat("Saving output to ", output)
saveRDS(Sce, file = output)
