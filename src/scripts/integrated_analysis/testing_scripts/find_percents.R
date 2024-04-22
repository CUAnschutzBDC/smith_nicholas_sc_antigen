library(here)
library(scAnalysisR)
library(pheatmap)
library(tidyverse)
library(splitstackshape)
library(circlize)
library(viridis)
library(ggalluvial)
library(Seurat)
library(djvdj)

normalization_method <- "log" # can be SCT or log
# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

#args <- commandArgs(trailingOnly = TRUE)

args <- c("merged", here("results"), "", here("files/sample_info.tsv"))

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "merged"

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- args[[4]]
#sample_info <- here("files/sample_info.tsv")

# Set directories
save_dir <- file.path(results_dir, "R_analysis", sample)

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))

seurat_data_new <- readRDS(file.path(save_dir, "rda_obj", "seurat_adtnorm.rds"))


cm <- confusionMatrix(seurat_data$Status, seurat_data$tet_name_cutoff)
cm <- cm / rowSums(cm) * 100


percent_total <- table(seurat_data$tet_name_cutoff) / ncol(seurat_data) * 100

cm_new <- confusionMatrix(seurat_data_new$Status,
                          seurat_data_new$tet_name_cutoff)
cm_new <- cm_new / rowSums(cm_new) * 100

percent_total_new <- table(seurat_data_new$tet_name_cutoff) / ncol(seurat_data_new) * 100
