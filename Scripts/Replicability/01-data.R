suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
})

# Helper functions ----
load_data <- function() {
  readRDS(here("Data", "full_data.rds"))
}

# Load dune data ----
load_Dune_labels <- function(cell_names, data_path = here("Data"),
                             size = "") {

  input_dir <- file.path(data_path, "Dune")
    label_matrix <- bind_rows(
      baron = read.csv(
        file.path(input_dir, paste0("baron_", size, "_Dune.csv"))),
      segerstolpe = read.csv(
        file.path(input_dir, paste0("segerstolpe_", size, "_Dune.csv"))),
      .id = "dataset"
    )
  
  # reorder cells to match data
  row_match <- match(cell_names, label_matrix$cells)
  label_matrix <- label_matrix[row_match, ]
  
  return(label_matrix)
}

load_Dune_NMI_labels <- function(cell_names, data_path = here("Data"),
                                 size = "") {
  
  input_dir <- file.path(data_path, "Dune")
  label_matrix <- bind_rows(
    baron = read.csv(
      file.path(input_dir, paste0("baron_", size, "_Dune_NMI.csv"))),
    segerstolpe = read.csv(
      file.path(input_dir, paste0("segerstolpe_", size, "_Dune_NMI.csv"))),
    .id = "dataset"
  )
  
  # reorder cells to match data
  row_match <- match(cell_names, label_matrix$cells)
  label_matrix <- label_matrix[row_match, ]
  
  return(label_matrix)
}

# Load hierarchical ----

load_single_merge_labels <- function(cell_names, data_path = here("Data"),
                                     size = "", type = "DE") {
  input_dir <- file.path(data_path, "Dune")
  result <- bind_rows(
    baron = read.csv(
      file.path(input_dir, 
                paste0("baron_", size , "_hierarchical_", type,
                       ".csv"))),
    segerstolpe = read.csv(
      file.path(input_dir,
                paste0("segerstolpe_", size , "_hierarchical_", type,
                       ".csv"))),
    .id = "dataset"
  )

  # NOTE: no cell ids, we assume that the order of cells is same as data

  # restrict to steps where both datasets have at least 2 clusters
  n_clusters_baron <- apply(result[result$dataset == "baron", ], 2,
                            function(x) length(table(x)))
  n_clusters_segerstolpe <- apply(result[result$dataset == "segerstolpe", ], 2,
                             function(x) length(table(x)))
  keep <- n_clusters_baron > 1 & n_clusters_segerstolpe > 1
  keep[1] <- TRUE
  result <- result[, keep]
  result <- result %>% select(dataset, cells, colnames(result))
  
  # rename columns for compatibility with analysis/visualization modules
  # (format METHOD_NAME.MERGING_STEP)
  methods <- colnames(result)[-(1:2)]
  if (type == "DE")  {
    method_name <- stringr::word(methods, 1, sep = stringr::fixed("."))
    method_level <- stringr::word(methods, 3, sep = stringr::fixed("."))
    method_level[is.na(method_level)] <- "00"
    colnames(result)[-(1:2)] <- paste(method_name, method_level, sep = ".")
  } else {
    methods <- stringr::word(methods, 2, sep = stringr::fixed("."))
    colnames(result)[-(1:2)] <- methods
  }
  
  # reorder cells to match data
  row_match <- match(cell_names, result$cells)
  result <- result[row_match, ]
  
  return(result)
}

# Load single Method ----
read_single_method <- function(filename) {
  result <- read.csv(filename, check.names = FALSE, row.names = 1)
  if ("cells" %in% colnames(result)) {
    result <- as_tibble(result)
  } else {
    result <- as_tibble(result, rownames = "cells")  
  }
  return(result %>% select(cells, colnames(result)))
}

load_single_method <- function(cell_names, data_path = here("Data"), method) {
  input_dir <- file.path(data_path, "singleMethod")
  result <- bind_rows(
    baron = read_single_method(
      file.path(input_dir, paste0("baron_", method, ".csv"))),
    segerstolpe = read_single_method(
      file.path(input_dir, paste0("segerstolpe_", method, ".csv"))),
    .id = "dataset"
  )
  
  # restrict to steps where both datasets have at least 2 clusters
  n_clusters_baron <- apply(result[result$dataset == "baron", ], 2,
                            function(x) length(table(x)))
  n_clusters_segerstolpe <- apply(result[result$dataset == "segerstolpe", ], 2,
                             function(x) length(table(x)))
  keep <- n_clusters_baron > 1 & n_clusters_segerstolpe > 1
  keep[1] <- TRUE
  result <- result[, keep]
  
  # rename columns for compatibility with analysis/visualization modules
  # (format METHOD_NAME.PARAMETER)
  colnames(result)[-(1:2)] <- paste(method,
                                    colnames(result)[-(1:2)], sep = ".")
  
  # reorder cells to match data
  row_match <- match(cell_names, result$cells)
  result <- result[row_match, ]
  
  return(as.data.frame(result))
}

load_single_seurat_labels <- function(cell_names, data_path = here("Data")) {
  return(load_single_method(cell_names, data_path, method = "Seurat"))
}

load_single_sc3_labels <- function(cell_names, data_path = here("Data")) {
  return(load_single_method(cell_names, data_path, method = "SC3"))
}

load_single_monocle_labels <- function(cell_names, data_path = here("Data")) {
  return(load_single_method(cell_names, data_path, method = "Monocle"))
}

# Load garbage ----
load_garbage_labels <- function(cell_names, data_path = here("Data"),
                                comp = 1, rep = 1, size = 1) {
  
  input_dir <- file.path(data_path, "Garbage")
  label_matrix <- bind_rows(
    baron = read.csv(
      file.path(input_dir, paste0("Garb", comp, "_baron_", rep, "_bad", size, ".csv"))),
    segerstolpe = read.csv(
      file.path(input_dir, paste0("Garb", comp, "_segerstolpe_", rep, "_bad", size, ".csv"))),
    .id = "dataset"
  )
  
  # reorder cells to match data
  row_match <- match(cell_names, label_matrix$cells)
  label_matrix <- label_matrix[row_match, ]
  
  return(label_matrix)
}

# Load downsampling ----
load_down_labels <- function(cell_names, data_path = here("Data"),
                             comp = 1, fraction = .05) {
  
  input_dir <- file.path(data_path, "Downsampling")
  label_matrix <- bind_rows(
    baron = read.csv(
      file.path(input_dir, paste0("baron_comp", comp, "_", 100 * fraction, ".csv"))),
    segerstolpe = read.csv(
      file.path(input_dir, paste0("segerstolpe_comp", comp, "_", 100 * fraction, ".csv"))),
    .id = "dataset"
  )
  
  # reorder cells to match data
  row_match <- match(cell_names, label_matrix$cells)
  label_matrix <- label_matrix[row_match, ]
  
  return(label_matrix)
}