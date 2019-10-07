# Packages and helper scripts ----
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(here)
})

source(here("Scripts", "Replicability", "meta_components.R"))
source(here("Scripts", "Replicability", "graph_visualization.R"))
source(here("Scripts", "Replicability", "01-data.R"))

# Helper functions ---- 
compute_replicability <- function(dataset, label_matrix, output_dir) {
  if (!file.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  label_matrix <- dplyr::select(label_matrix, -dataset, -cells)
  label_sets <- colnames(label_matrix)
  
  is_not_noise <- dataset$class_label != "Noise"
  is_hvg <- rowData(dataset)$is_hvg
  
  my_dataset <- dataset[is_hvg, is_not_noise]
  my_labels <- label_matrix[is_not_noise, ]
  for (current_set in label_sets) {
    labels <- my_labels[, current_set]
    is_not_na <- !is.na(labels)
    if (sum(is_not_na) == 0) next
    stats <- analyze_components(my_dataset[, is_not_na], labels[is_not_na])
    export_components(stats, file.path(output_dir, current_set))
  }
}

analyze_components <- function(dataset, labels) {
  is_nonzero_cell <- Matrix::colSums(assay(dataset)) != 0
  dataset <- dataset[, is_nonzero_cell]
  labels <- labels[is_nonzero_cell]
  
  best_hits <- compute_best_hits(dataset, labels)
  components <- extract_components(best_hits, 0.6)
  
  return(list(
    best_hits = best_hits,
    components = components,
    labels = paste(dataset$study_id, labels, sep = "|")
  ))
}

export_components <- function(component_obj, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  write.table(component_obj$best_hits, file.path(output_dir, "best_hits.txt"))
  plot_components(component_obj$best_hits, component_obj$components$modules,
                  output_dir)
  write_component_summary(component_obj$components, output_dir)
  export_filtered_best_hits(component_obj$best_hits,
                            component_obj$components$modules, output_dir)
  
  graph <- make_directed_graph(component_obj$best_hits, 0.6, 1)
  graph <- color_graph(graph, component_obj$labels)
  pdf(file.path(output_dir, "graph_visualization.pdf"))
  plot_directed_graph(graph, 1)
  dev.off()
}

analyze_data <- function(dataset, label_matrix, output_dir) {
  compute_replicability(dataset, label_matrix, file.path(output_dir))
}

# Main functions ----
## Dune ----
analyze_Dunes <- function(data_path = here("Data"),
                          output_dir = here("Data", "Replicability",
                                            "mn_results", "Dune")) {
  dataset <- load_smart_data()
  # Dune comp1
  labels <- load_Dune_labels(colnames(dataset), data_path)
  analyze_data(dataset, labels, paste0(output_dir, "/comp1"))
  # Dune comp2
  labels <- load_Dune_labels(colnames(dataset), data_path, size = "comp2")
  analyze_data(dataset, labels, paste0(output_dir, "/comp2"))
  # Dune comp3
  labels <- load_Dune_labels(colnames(dataset), data_path, size = "comp3")
  analyze_data(dataset, labels, paste0(output_dir, "/comp3"))
}

## Hierarchical ----
analyze_single_merge <- function(data_path = here("Data"),
                                 output_dir = here("Data", "Replicability",
                                                   "mn_results", "singleTree")) {
  dataset <- load_smart_data()
  # DE
  ## comp1 single Merge
  labels <- load_single_merge_labels(colnames(dataset), data_path)
  analyze_data(dataset, labels, paste0(output_dir, "/comp1_DE"))
  ## Comparison 2
  labels <- load_single_merge_labels(colnames(dataset), data_path,
                                     size = "_comp2")
  analyze_data(dataset, labels, paste0(output_dir, "/comp2_DE"))
  ## Comparison 3
  labels <- load_single_merge_labels(colnames(dataset), data_path,
                                     size = "_comp3")
  analyze_data(dataset, labels, paste0(output_dir, "/comp3_DE"))
  # Dist
  ## comp1 single Merge
  labels <- load_single_merge_labels(colnames(dataset), data_path, type = "Dist")
  analyze_data(dataset, labels, paste0(output_dir, "/comp1_Dist"))
  ## Comparison 2
  labels <- load_single_merge_labels(colnames(dataset), data_path,
                                     size = "_comp2", type = "Dist")
  analyze_data(dataset, labels, paste0(output_dir, "/comp2_Dist"))
  ## Comparison 3
  labels <- load_single_merge_labels(colnames(dataset), data_path,
                                     size = "_comp3", type = "Dist")
  analyze_data(dataset, labels, paste0(output_dir, "/comp3_Dist"))
}

## single Method ----
analyze_single_methods <- function(
  data_path = here("Data"),
  output_dir = here("Data", "Replicability", "mn_results", "SingleMethod")) 
  {
  dataset <- load_smart_data()
  # Seurat
  labels <- load_single_seurat_labels(colnames(dataset), data_path)
  analyze_data(dataset, labels, output_dir)
  # SC3
  labels <- load_single_sc3_labels(colnames(dataset), data_path)
  analyze_data(dataset, labels, output_dir)
  # Monocle
  labels <- load_single_monocle_labels(colnames(dataset), data_path)
  analyze_data(dataset, labels, output_dir)
}

## To run ----
main <- function() {
  print("single Method smart")
  analyze_single_methods()
  print("single merge")
  analyze_single_merge()
  print("single All Dunes")
  analyze_Dunes()
}

if (!interactive()) {
  main()
}