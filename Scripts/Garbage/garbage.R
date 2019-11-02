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
  make_option(c("-n", "--nCores"),
              action = "store", default = 1, type = "integer",
              help = "Number of cores to use"
  ),
  make_option(c("-S", "--SeuratParam"),
              action = "store", default = NA, type = "character",
              help = "Parameter to use for Seurat"
  ),
  make_option(c("-C", "--C3"),
              action = "store", default = NA, type = "character",
              help = "SC3 parameter"
  ),
  make_option(c("-m", "--monocle"),
              action = "store", default = NA, type = "character",
              help = "Monocle parameter"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$l)) {
  loc <- opt$l
  cat("The selected dataset is located at ", loc, "\n")
} else {
  stop("Missing l argument")
}

if (!is.na(opt$o)) {
  output <- opt$o
} else {
  stop("Missing o argument")
  cat("The output will be stored at ", output, "\n")
}

if (!is.na(opt$C)) {
  sc3_p <- opt$C
} else {
  stop("Missing C argument")
}

if (!is.na(opt$C)) {
  seurat_p <- opt$S
} else {
  stop("Missing S argument")
}

if (!is.na(opt$m)) {
  monocle_p <- opt$m
} else {
  stop("Missing m argument")
}

library(SummarizedExperiment)
library(parallel)
library(matrixStats)
library(clusterExperiment)
library(tidyverse)
library(Dune)
library(mclust)

# Load Data ----
# Load sc3 clustering results
sc3 <- read.csv(paste0(loc, "_SC3.csv"))[, -1]
colnames(sc3) <- str_remove(colnames(sc3), "^X") %>% str_replace("\\.", "-")
Names <- sc3$cells
sc3 <- sc3[, sc3_p] %>% as.numeric()

# Load Seurat clustering results
Seurat <- read.csv(paste0(loc, "_Seurat.csv"))[, -1]
colnames(Seurat) <- str_remove(colnames(Seurat), "^X")
Seurat <- Seurat[, seurat_p] %>% as.numeric()

# Load Monocle clustering results
Monocle <- read.csv(paste0(loc, "_Monocle.csv"))[, -1]
Monocle <- as.data.frame(Monocle)[, monocle_p] %>% as.numeric()

# Do the consensus clustering with one additional garbage clustering ----
clusMat <- data.frame("sc3" = sc3, "Monocle" = Monocle, "Seurat" = Seurat)
rownames(clusMat) <- Names  
n_clus <- mean(n_distinct(sc3), n_distinct(Seurat), n_distinct(Monocle))
clusMat$garbage <- sample(n_clus, size = nrow(clusMat), replace = TRUE)
print(paste0("Number of cores: ", opt$n))
print(system.time(
  merger <- Dune(clusMat = clusMat, nCores = opt$n)
))
cat("Finished Consensus Merge\n")
saveRDS(merger,  paste0(output, "_1-bad_merger.rds"))

# Save the matrix with all the consensus steps ----
print("...Initial")
initialMat <- merger$initialMat
initialMat <- as.matrix(initialMat) 

print("...Final consensus")
currentMat <- merger$currentMat
currentMat <- as.matrix(currentMat) 

print("...Intermediary consensus at 33.3%")
stopMatrix_33 <- intermediateMat(merger = merger,
                                 p = 1/3)
stopMatrix_33 <- as.matrix(stopMatrix_33)

print("...Intermediary consensus at 66.7%")
stopMatrix_66 <- intermediateMat(merger = merger, p = 2/3)
stopMatrix_66 <- as.matrix(stopMatrix_66)

print("...Intermediary consensus at 90%")
stopMatrix_90 <- intermediateMat(merger = merger, p = .9)
stopMatrix_90 <- as.matrix(stopMatrix_90)

print("...Full matrix")
mat <- cbind(as.character(Names),
             initialMat,  stopMatrix_33, stopMatrix_66,  stopMatrix_90,
             currentMat)

chars <- c("sc3", "Monocle", "Seurat", "garbage")

colnames(mat) <- c("cells",
                   paste(chars, "Initial", sep = "-"), paste(chars, "33", sep = "-"),
                   paste(chars, "66", sep = "-"), paste(chars, "90", sep = "-"),
                   paste(chars, "Final", sep = "-")
)

write_csv(x = as.data.frame(mat), path = paste0(output, "_1-bad.csv"))

# Do the consensus clustering with two additional garbage clustering ----
clusMat <- data.frame("sc3" = sc3, "Monocle" = Monocle, "Seurat" = Seurat)
rownames(clusMat) <- Names  
n_clus <- mean(n_distinct(sc3), n_distinct(Seurat), n_distinct(Monocle))
clusMat$garbage1 <- sample(n_clus, size = nrow(clusMat), replace = TRUE)
clusMat$garbage2 <- sample(n_clus, size = nrow(clusMat), replace = TRUE)
print(paste0("Number of cores: ", opt$n))
print(system.time(
  merger <- Dune(clusMat = clusMat, nCores = opt$n)
))
cat("Finished Consensus Merge\n")
saveRDS(merger,  paste0(output, "_2-bad_merger.rds"))

# Save the matrix with all the consensus steps ----
print("...Initial")
initialMat <- merger$initialMat
initialMat <- as.matrix(initialMat) 

print("...Final consensus")
currentMat <- merger$currentMat
currentMat <- as.matrix(currentMat) 

print("...Intermediary consensus at 33.3%")
stopMatrix_33 <- intermediateMat(merger = merger,
                                 p = 1/3)
stopMatrix_33 <- as.matrix(stopMatrix_33)

print("...Intermediary consensus at 66.7%")
stopMatrix_66 <- intermediateMat(merger = merger, p = 2/3)
stopMatrix_66 <- as.matrix(stopMatrix_66)

print("...Intermediary consensus at 90%")
stopMatrix_90 <- intermediateMat(merger = merger, p = .9)
stopMatrix_90 <- as.matrix(stopMatrix_90)

print("...Full matrix")
mat <- cbind(as.character(Names),
             initialMat,  stopMatrix_33, stopMatrix_66,  stopMatrix_90,
             currentMat)

chars <- c("sc3", "Monocle", "Seurat", "garbage1", "garbage2")

colnames(mat) <- c("cells",
                   paste(chars, "Initial", sep = "-"), paste(chars, "33", sep = "-"),
                   paste(chars, "66", sep = "-"), paste(chars, "90", sep = "-"),
                   paste(chars, "Final", sep = "-")
)

write_csv(x = as.data.frame(mat), path = paste0(output, "_2-bad.csv"))

# Do the consensus clustering with three additional garbage clustering ----
clusMat <- data.frame("sc3" = sc3, "Monocle" = Monocle, "Seurat" = Seurat)
rownames(clusMat) <- Names  
n_clus <- mean(n_distinct(sc3), n_distinct(Seurat), n_distinct(Monocle))
clusMat$garbage1 <- sample(n_clus, size = nrow(clusMat), replace = TRUE)
clusMat$garbage2 <- sample(n_clus, size = nrow(clusMat), replace = TRUE)
clusMat$garbage3 <- sample(n_clus, size = nrow(clusMat), replace = TRUE)
print(paste0("Number of cores: ", opt$n))
print(system.time(
  merger <- Dune(clusMat = clusMat, nCores = opt$n)
))
cat("Finished Consensus Merge\n")
saveRDS(merger,  paste0(output, "_3-bad_merger.rds"))

# Save the matrix with all the consensus steps ----
print("...Initial")
initialMat <- merger$initialMat
initialMat <- as.matrix(initialMat) 

print("...Final consensus")
currentMat <- merger$currentMat
currentMat <- as.matrix(currentMat) 

print("...Intermediary consensus at 33.3%")
stopMatrix_33 <- intermediateMat(merger = merger,
                                 p = 1/3)
stopMatrix_33 <- as.matrix(stopMatrix_33)

print("...Intermediary consensus at 66.7%")
stopMatrix_66 <- intermediateMat(merger = merger, p = 2/3)
stopMatrix_66 <- as.matrix(stopMatrix_66)

print("...Intermediary consensus at 90%")
stopMatrix_90 <- intermediateMat(merger = merger, p = .9)
stopMatrix_90 <- as.matrix(stopMatrix_90)

print("...Full matrix")
mat <- cbind(as.character(Names),
             initialMat,  stopMatrix_33, stopMatrix_66,  stopMatrix_90,
             currentMat)

chars <- c("sc3", "Monocle", "Seurat", "garbage1", "garbage2", "garbage3")

colnames(mat) <- c("cells",
                   paste(chars, "Initial", sep = "-"), paste(chars, "33", sep = "-"),
                   paste(chars, "66", sep = "-"), paste(chars, "90", sep = "-"),
                   paste(chars, "Final", sep = "-")
)

write_csv(x = as.data.frame(mat), path = paste0(output, "_3-bad.csv"))

