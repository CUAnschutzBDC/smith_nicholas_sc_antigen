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

seurat_object <- HTODemux(seurat_data, assay = "SCAR_TET", 
                         positive.quantile = 0.90,
                         kfunc = "kmeans")


all_counts <- GetAssayData(seurat_data, assay = "SCAR_TET",
                           slot = "counts") %>%
  as.matrix() %>%
  t() %>%
  data.frame

meta_data <- seurat_data[[]] %>%
  dplyr::select(SCAR_TET_classification, hash.ID) %>%
  merge(all_counts, by = "row.names")

object <- seurat_data
assay <- "SCAR_TET"
positive.quantile <- 0.90
init <- NULL
nstarts <- 100
kfunc <- "kmeans"
nsamples <- 100
seed <- 42
verbose <- TRUE


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
}
npositive <- colSums(x = discrete)
classification.global <- npositive
classification.global[npositive == 0] <- "Negative"
classification.global[npositive == 1] <- "Singlet"
classification.global[npositive > 1] <- "Doublet"

# This cell is wrong AAACCTGCACGACGAA-1

donor.id = rownames(x = data)
hash.max <- apply(X = data, MARGIN = 2, FUN = max)
hash.maxID <- apply(X = data, MARGIN = 2, FUN = which.max)
hash.second <- apply(X = data, MARGIN = 2, FUN = Seurat:::MaxN, N = 2)
hash.maxID <- as.character(x = donor.id[sapply(X = 1:ncol(x = data), 
                                               FUN = function(x) {
                                                 return(which(x = data[, x] == hash.max[x])[1])
                                               })])
hash.secondID <- as.character(x = donor.id[sapply(X = 1:ncol(x = data), 
                                                  FUN = function(x) {
                                                    return(which(x = data[, x] == hash.second[x])[1])
                                                  })])

# The problem is that the "max" is determined regardless of the background
# levels

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
