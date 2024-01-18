library(here)
library(scAnalysisR)
library(pheatmap)
library(tidyverse)
library(splitstackshape)
library(viridis)
library(Seurat)

normalization_method <- "log" # can be SCT or log
# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

#sample <- args[[1]]
#sample <- gsub("__.*", "", sample)
sample <- "merged"

#results_dir <- args[[2]]
results_dir <- here("results")

#sample_info <- args[[4]]
sample_info <- here("files/sample_info.tsv")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

sample_info <- sample_info[sample_info$sample == sample,]
vars.to.regress <- NULL
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
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))

# Plots to make
# 1 frequency of cells that are "doublets"
scAnalysisR::stacked_barplots(seurat_data, meta_col = scar_hash_id)

# 2 Doublets broken down by type
# 3 Number and or frequency of GAD-insulin cells in each type
# 4 Expression of insulin receptor


# Expression by sample
all_plots <- lapply(unique(seurat_data$sample), function(x){
  seurat_sub <- subset(seurat_data, subset = sample == x)
  return_plot <- featDistPlot(seurat_sub, "INS-tet", assay = "SCAR_TET_LOG",
                              sep_by = "scar_hash_id", combine = FALSE)
  
  return(return_plot)
  
})

featDistPlot(seurat_data, "INS-tet", assay = "SCAR_TET_LOG",
             sep_by = "sample", combine = FALSE, col_by = "scar_hash_id")

featDistPlot(seurat_data, "INS-tet", assay = "SCAR_TET",
             sep_by = "sample", combine = FALSE, col_by = "scar_hash_id")

featDistPlot(seurat_data, "INS-tet", assay = "SCAR_TET",
             sep_by = "SCAR_TET_classification", combine = FALSE,
             col_by = "SCAR_TET_classification")

# Look into EB. CH and CP also don't seem quite right.

# Look into how you see doublets that are Ins Ins