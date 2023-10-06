library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(clustree)
library(harmony)
library(singlecellmethods)
library(batchelor)


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

remove_ambience <- TRUE

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

#sample <- args[[1]]
#sample <- gsub("__.*", "", sample)
sample <- "merged"

#sample_info <- args[[4]]
sample_info <- here("files/sample_info.tsv")

#results_dir <- args[[2]]
results_dir <- here("results")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

samples_use <- sample_info[sample_info$sample != sample, ]$sample

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident


ADT_pca <- sample_info$adt_pca
run_adt_umap <- sample_info$adt_umap

RNA_pcs <- sample_info$PCs
resolution <- sample_info$resolution
batch_correction <- sample_info$batch_correction



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
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

pca_list <- list("rna" = c(assay = seurat_assay, reduction_name = "pca"))

if(remove_ambience){
  pca_list <- c(pca_list,
                list("ambience" = c(assay = "AMBRNA", 
                                    reduction_name = "ambpca")))
} 


# UMAP -------------------------------------------------------------------------
for(assay_type in names(pca_list)){
  seurat_assay <- pca_list[[assay_type]][["assay"]]
  reduction_name <- pca_list[[assay_type]][["reduction_name"]]

  reduction_save <- paste(assay_type, batch_correction, sep = "_")
  
  # Remove previous clustering
  remove_cols <- colnames(seurat_data[[]])[grepl("res\\.[0-9]",
                                                 colnames(seurat_data[[]]))]
  
  for (i in remove_cols){
    seurat_data[[i]] <- NULL
  }
  
  set.seed(0)
  umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                           resolution = 0.8, assay = seurat_assay, HTO = HTO,
                           reduction = reduction_save)
  ## Cluster --------------------------------------------------------------------
  seurat_data <- FindClusters(seurat_data, resolution = c(0.2, 0.5, 0.6,
                                                          0.8, 1, 1.2))
  
  clustering_columns <- colnames(seurat_data[[]])[grepl("RNA_snn_res",
                                                        colnames(seurat_data[[]]))]
  
  quality_columns <- c("nCount_RNA", "nFeature_RNA", "nCount_ADT",
                       "nFeature_ADT", "percent.mt", "Doublet_finder",
                       "nCount_TET", "nFeature_TET")
  
  plot_columns <- c(clustering_columns, quality_columns)
  
  all_plots <- plotDimRed(seurat_data, col_by = plot_columns,
                          plot_type = paste0(reduction_save, ".umap"),
                          ggrastr = TRUE)
  
  # Save these resolutions to a file
  pdf(file.path(save_dir, "images", 
                paste0(assay_type, "clustering_resolution.pdf")))
  print(clustree(seurat_data))
  print(all_plots)
  
  dev.off()
}