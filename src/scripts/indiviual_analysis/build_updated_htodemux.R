# Look into EB. CH and CP also don't seem quite right.

# Look into how you see doublets that are Ins Ins

library(DoubletFinder)
library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

#sample <- args[[1]]
#sample <- gsub("__.*", "", sample)
sample <- "JH310-12_EB"

#results_dir <- args[[2]]
results_dir <- here("results")

#sample_info <- args[[4]]
sample_info <- here("files/sample_info.tsv")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident


if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
save_dir <- file.path(results_dir, "R_analysis", sample)

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_start.rds"))

HTODemuxUpdate <- function(object, assay = "HTO", positive.quantile = 0.99,
                           init = NULL, nstarts = 100, kfunc = "clara", 
                           nsamples = 100, seed = 42, verbose = TRUE){
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay)
  counts <- GetAssayData(object = object, assay = assay, slot = "counts")[, 
                                                                          colnames(x = object)]
  counts <- as.matrix(x = counts)
  ncenters <- init %||% (nrow(x = data) + 1)
  switch(EXPR = kfunc, kmeans = {
    init.clusters <- kmeans(x = t(x = GetAssayData(object = object, 
                                                   assay = assay)), centers = ncenters, nstart = nstarts)
    Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
  }, clara = {
    init.clusters <- clara(x = t(x = GetAssayData(object = object, 
                                                  assay = assay)), k = ncenters, samples = nsamples)
    Idents(object = object, cells = names(x = init.clusters$clustering), 
           drop = TRUE) <- init.clusters$clustering
  }, stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'"))
  average.expression <- AverageExpression(object = object, 
                                          assays = assay, verbose = FALSE)[[assay]]
  if (sum(average.expression == 0) > 0) {
    stop("Cells with zero counts exist as a cluster.")
  }
  discrete <- GetAssayData(object = object, assay = assay)
  discrete[discrete > 0] <- 0
  
  proportion <- discrete
  for (iter in rownames(x = data)) {
    values <- counts[iter, colnames(object)]
    values.use <- values[WhichCells(object = object, idents = levels(x = Idents(object = object))[[which.min(x = average.expression[iter, 
    ])]])]
    fit <- suppressWarnings(expr = fitdistrplus::fitdist(data = values.use, 
                                                         distr = "nbinom"))
    cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
    discrete[iter, names(x = which(x = values > cutoff))] <- 1
    if (verbose) {
      message(paste0("Cutoff for ", iter, " : ", cutoff, 
                     " reads"))
    }
    
    proportion[iter,] <- values / cutoff
  }
  npositive <- colSums(x = discrete)
  classification.global <- npositive
  classification.global[npositive == 0] <- "Negative"
  classification.global[npositive == 1] <- "Singlet"
  classification.global[npositive > 1] <- "Doublet"
  
  # Here I do it based on the proportion (count / cutoff). This way, if
  # the value was higher than the cutoff it will be greater than 1 and if it
  # is lower it will be less than 1.
  donor.id <- rownames(x = proportion)
  hash.max <- apply(X = proportion, MARGIN = 2, FUN = max)
  hash.maxID <- sapply(X = 1:ncol(x = proportion), FUN = function(x) {
    max_indices <- which(proportion[, x] == hash.max[x])
    return(donor.id[min(max_indices)])
  })
  
  hash.second <- apply(X = proportion, MARGIN = 2, FUN = function(x) {
    sorted_values <- sort(x, decreasing = TRUE)
    if (length(sorted_values) > 1) {
      second_max <- sorted_values[2] # Get the second highest value
      return(second_max)
    } else {
      return(NA) # Handle cases where there is no second highest value
    }
  })
  
  hash.secondID <- sapply(X = 1:ncol(x = proportion), FUN = function(x) {
    idx <- which(proportion[, x] == hash.second[x])
    if (length(idx) > 1) {
      # If there are multiple indices for the second highest value, choose the one not equal to hash.maxID
      idx <- idx[which(names(idx) != hash.maxID[x])]
    }
    return(donor.id[idx[1]])
  })
  
  
  hash.margin <- hash.max - hash.second
  doublet_id <- sapply(X = 1:length(x = hash.maxID), FUN = function(x) {
    return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), 
                 collapse = "_"))
  })
  classification <- classification.global
  classification[classification.global == "Negative"] <- "Negative"
  classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == 
                                                                           "Singlet")]
  classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == 
                                                                           "Doublet")]
  classification.metadata <- data.frame(hash.maxID, hash.secondID, 
                                        hash.margin, classification, classification.global)
  colnames(x = classification.metadata) <- paste(assay, c("maxID", 
                                                          "secondID", "margin", "classification", "classification.global"), 
                                                 sep = "_")
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste0(assay, "_classification")
  doublets <- rownames(x = object[[]])[which(object[[paste0(assay, 
                                                            "_classification.global")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- "Doublet"
  object$hash.ID <- Idents(object = object)
  
  
  return(object)
}


seurat_object <- HTODemuxUpdate(seurat_data, assay = "SCAR_TET", 
                                positive.quantile = 0.90,
                                kfunc = "kmeans")




all_counts <- GetAssayData(seurat_object, assay = "SCAR_TET",
                           slot = "counts") %>%
  as.matrix() %>%
  t() %>%
  data.frame

meta_data <- seurat_object[[]] %>%
  dplyr::select(SCAR_TET_classification, hash.ID) %>%
  merge(all_counts, by = "row.names")

table(meta_data$SCAR_TET_classification, meta_data$INS.tet > 54)
table(meta_data$SCAR_TET_classification, meta_data$TET.tet > 13)
table(meta_data$SCAR_TET_classification, meta_data$IA2.tet > 7)
table(meta_data$SCAR_TET_classification, meta_data$DNA.tet > 6)
table(meta_data$SCAR_TET_classification, meta_data$GAD.tet > 13)


meta_data[meta_data$SCAR_TET_classification == "DNA-tet_GAD-tet" & 
            meta_data$INS.tet > 54,]


meta_data[meta_data$SCAR_TET_classification == "INS-tet_TET-tet" & 
            meta_data$INS.tet <= 54,]


