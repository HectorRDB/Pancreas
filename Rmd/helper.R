libs <- c("here", "dplyr", "ggplot2", "tidyr", "stringr", "readr", "mclust",
          "RColorBrewer", "purrr", "Dune")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)
# utils functions ----
toRank <- function(i) {
  case_when(
    i == "Initial" ~ 1,
    i == 33 ~ 2,
    i == 66 ~ 3,
    i == 90 ~ 4,
    i == "Final" ~ 5
  )
}

comp_with_ref <- function(x, y) {
  return(c("clusters" = n_distinct(x), "ARI" = adjustedRandIndex(x, y)))
}

# With hierarchical ----
load_dune <- function(dataset, comp, ref) {
  df <- read.csv(here("Data", "Dune",
                        paste0(dataset, "_", comp, "_Dune.csv"))) %>%
    arrange(cells) %>%
    select(-cells) %>%
    map_df(., comp_with_ref, y = ref) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename("n_clus" = V1, "ari" = V2) %>%
    mutate(clustering_method = rownames(.),
           level = word(clustering_method, 2, sep = "\\."),
           clustering_method = word(clustering_method, 1, sep = "\\."))
  df$level <- lapply(df$level, toRank) %>% unlist()
  df$level <- as.numeric(df$level)
  return(df)
}

load_DE_tree <- function(dataset, comp, ref) {
  df <- read.csv(here("Data", "Dune",
                        paste0(dataset, "_", comp, "_hierarchical_DE.csv"))) %>%
    arrange(cells) %>%
    select(-cells) %>%
    map_df(., comp_with_ref, y = ref) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename("n_clus" = V1, "ari" = V2) %>%
    mutate(clustering_method = rownames(.),
           level = word(clustering_method, 2, sep = "\\_"),
           level = as.numeric(level),
           clustering_method = word(clustering_method, 1, sep = "\\."))
  return(df)
}

load_Dist_tree <- function(dataset, comp, ref) {
  df <- read.csv(here("Data", "Dune",
                      paste0(dataset, "_", comp, "_hierarchical_Dist.csv"))) %>%
    arrange(cells) %>%
    select(-cells) %>%
    map_df(., comp_with_ref, y = ref) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename("n_clus" = V1, "ari" = V2) %>%
    mutate(clustering_method = rownames(.),
           level = word(clustering_method, 2, sep = "\\_"),
           level = as.numeric(level),
           clustering_method = word(clustering_method, 1, sep = "\\."))
  return(df)
}

# With single methods ----
load_sc3 <- function(dataset, dune) {
  sc3_init <- dune[, "sc3.Initial"]
  SC3 <- read.csv(here("Data", "singleMethod",
                       paste0(dataset, "_SC3.csv"))) %>%
    select(-X) %>%
    arrange(cells) %>%
    select(-cells)
  colnames(SC3) <- str_remove(colnames(SC3), "X") %>%
    str_replace("\\.", "-") %>% unlist
  sc3_param <- colnames(SC3)[which.max(colSums(SC3 == sc3_init))] %>%
    as.numeric()
  keep <- as.numeric(colnames(SC3)) <= sc3_param
  SC3 <- SC3[, keep]
  colnames(SC3) <- paste("sc3", colnames(SC3), sep = "_")
  return(SC3)
}

load_monocle <- function(dataset, dune) {
  Monocle_init <- dune[, "Monocle.Initial"]
  Monocle <- read.csv(here("Data", "singleMethod",
                       paste0(dataset, "_Monocle.csv"))) %>%
    select(-X) %>%
    arrange(cells) %>%
    select(-cells)
  colnames(Monocle) <- str_remove(colnames(Monocle), "k_") %>% unlist()
  Monocle_param <- colnames(Monocle)[which.max(colSums(Monocle == Monocle_init))] %>%
    as.numeric()
  keep <- as.numeric(colnames(Monocle)) >= Monocle_param
  Monocle <- Monocle[, keep]
  colnames(Monocle) <- paste("Monocle", colnames(Monocle), sep = "_")
  return(Monocle)
}

load_seurat <- function(dataset, dune) {
  Seurat_init <- dune[, "Seurat.Initial"]
  Seurat <- read.csv(here("Data", "singleMethod",
                           paste0(dataset, "_Seurat.csv"))) %>%
    select(-X) %>%
    arrange(cells) %>%
    select(-cells)
  colnames(Seurat) <- str_remove(colnames(Seurat), "X") %>% unlist()
  Seurat_param <- colnames(Seurat)[which.max(colSums(Seurat == Seurat_init))]
  k_seurat <- word(Seurat_param, 3, sep = "\\.")
  keep <- str_detect(colnames(Seurat), k_seurat) %>% unlist()
  Seurat <- Seurat[, keep]
  colnames(Seurat) <- str_remove(colnames(Seurat), paste0(".", k_seurat))
  Seurat_param <- colnames(Seurat)[which.max(colSums(Seurat == Seurat_init))] %>%
    as.numeric()
  keep <- as.numeric(colnames(Seurat)) <= Seurat_param
  Seurat <- Seurat[, keep]
  colnames(Seurat) <- paste("Seurat", colnames(Seurat), sep = "_")
  return(Seurat)
}

load_singles <- function(dataset, comp, ref) {
  dune <- df <- read.csv(here("Data", "Dune",
                              paste0(dataset, "_", comp, "_Dune.csv"))) %>%
    arrange(cells)
  Seurat <- load_seurat(dataset, dune)
  SC3 <- load_sc3(dataset, dune)
  Monocle <- load_monocle(dataset, dune)
  df <- cbind(Seurat, SC3, Monocle) %>%
    as.data.frame() %>%
    map_df(., comp_with_ref, y = ref) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename("n_clus" = V1, "ari" = V2) %>%
    mutate(clustering_method = rownames(.),
           level = word(clustering_method, 2, sep = "_"),
           level = as.numeric(level),
           clustering_method = word(clustering_method, 1, sep = "_"))
  
}
# Main reports ----

plot_all <- function(dataset, comp, ref) {
  Dune <- readRDS(here("Data", "Dune", paste0(dataset, "_", comp, ".rds")))
  p1 <- plotPrePost(Dune)
  p2 <- plot_grid(
    plotARIs(clusMat = Dune$initialMat) + ggtitle("Before Merging"),
    plotARIs(clusMat = Dune$currentMat) + ggtitle("After Merging"),
    ncol = 2
  )
  p3 <- ARItrend(Dune)
  
  p4 <- plotComp(dataset, comp, ref)
  return(list(p1, p2, p3, p4))
}

plotComp <- function(dataset, comp, ref) {
  df1 <- bind_rows("Dune" = load_dune(dataset, comp, ref),
                   "DE" = load_DE_tree(dataset, comp, ref),
                   .id = "Method") %>%
    mutate(comp = "Dune versus Hierarchical DE")
  df2 <- bind_rows("Dune" = load_dune(dataset, comp, ref),
                   "Dist" = load_Dist_tree(dataset, comp, ref),
                   .id = "Method") %>%
    mutate(comp = "Dune versus Hierarchical Dist")
  df3 <- bind_rows("Dune" = load_dune(dataset, comp, ref),
                   "Param" = load_singles(dataset, comp, ref),
                   .id = "Method") %>%
    mutate(comp = "Dune versus new parameters")
  df <- bind_rows(df1, df2, df3) %>%
    arrange(level)
  
  p <- ggplot(df %>% filter(n_clus > 10),
              aes(x = n_clus,
                  y = ari,
                  linetype = clustering_method,
                  col = Method,
                  group = interaction(clustering_method, Method))) +
    geom_path(size = 1.8) +
    theme_classic() +
    # scale_linetype_manual(values = linetypes) +
    scale_color_brewer(type = "qual") +
    labs(x = "Resolution",
         y = "ARI with glod standard",
         linetype = "Clustering\nmethod",
         col = "Method of\nmerging") +
    facet_wrap(~comp) +
  NULL
  return(p)
}
