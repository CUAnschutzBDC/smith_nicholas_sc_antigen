library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)
library(here)
library(scAnalysisR)
library(viridis)
library(clustifyr)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

remove_ambience <- TRUE

normalization_method <- "log" # can be SCT or log

cor_cutoff <- 0.5

all_ref_dir <-
  "/beevol/home/wellskri/Analysis/references/single_cell_references"

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "merged"

sample_info <- args[[4]]
#sample_info <- here("files/sample_info.tsv")

results_dir <- args[[2]]
#results_dir <- here("results")

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

#-------------------------------------------------------------------------------
for(assay_type in names(pca_list)){
  seurat_assay <- pca_list[[assay_type]][["assay"]]
  reduction_name <- pca_list[[assay_type]][["reduction_name"]]

  plot_type <- paste0(assay_type, "_", batch_correction, ".umap")
  
  # NOTE update the section below to include any references pertanent to your
  # sample.
  
  ################
  # Seurat pbmc #
  ###############
  
  # Information for cell mapping
  ref_dir <- file.path(all_ref_dir, "pbmc/seurat")
  
  ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference_l2.csv"),
                      header = TRUE, row.names = 1)
  
  pdf(file.path(save_dir, "images", 
                paste0(assay_type, "_celltype_mapping.pdf")))
  
  cluster_res <- name_clusters(seurat_data, ref_mat,
                               save_dir = save_dir,
                               save_name = "comb_celltype_seurat", ADT = FALSE,
                               assay = seurat_assay,
                               nfeatures = 1000, 
                               clusters = paste0(seurat_assay, "_cluster"),
                               plot_type = plot_type,
                               features = VariableFeatures(seurat_data))
  
  seurat_data <- cluster_res$object
  
  seurat_res_seurat <- cluster_res$RNA
  
  grid::grid.newpage()
  
  print(plotDimRed(seurat_data, 
                   col_by = paste0(seurat_assay, "_comb_celltype_seurat"),
                   plot_type = plot_type, ggrastr = TRUE))
  
  cm <- confusionMatrix(seurat_data[[paste0(seurat_assay,
                                            "_comb_celltype_seurat")]][[1]],
                        seurat_data$RNA_celltype_seurat)
  cm <- cm / rowSums(cm)
  
  pheatmap::pheatmap(cm)
  
  write.csv(seurat_res_seurat, 
            file.path(save_dir, "files", 
                      paste0(assay_type, "_celltype_mapping_seurat.csv")))
  
  
  #-------------------------------------------------------------------------------
  
  ####################
  # Published object #
  ####################
  
  ref_dir <- here("files/210825_object")
  
  
  ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference.csv"),
                      header = TRUE, row.names = 1)
  
  grid::grid.newpage()
  
  cluster_res <- name_clusters(seurat_data, ref_mat,
                               save_dir = save_dir,
                               save_name = "comb_celltype_bnd", ADT = FALSE,
                               assay = seurat_assay,
                               features = VariableFeatures(seurat_data),
                               clusters = paste0(seurat_assay, "_cluster"),
                               plot_type = plot_type)
  
  seurat_data <- cluster_res$object
  
  seurat_res_bnd <- cluster_res$RNA
  
  grid::grid.newpage()
  
  new_colors <- c("Resting_memory" = "#924bdb", # Resting memory
                  "Naive_1" = "#69ba3d", # Naive 1
                  "Naive_2" = "#9a43a4", # Naive 2
                  "Memory_IgE_IgG" = "#bf9b31", # Memory IgE/IgG1
                  "Naive_3" = "#6477ce", # Naive 3
                  "Memory_IgA" = "#d15131", # Memory IA
                  "Early_memory" = "#4c9e8e", # Early Memory
                  "BND2" = "#cc4570", #Bnd2
                  "DN2" = "#648d4f", # DN2
                  "Activated_memory" = "#985978", # Activated memory
                  "Activated_naive" = "#a06846") # Activated naive
  
  print(plotDimRed(seurat_data,
                   col_by = paste0(seurat_assay, "_comb_celltype_bnd"),
                   plot_type = plot_type, color = new_colors,
                   ggrastr = TRUE))
  
  cm <- confusionMatrix(seurat_data[[paste0(seurat_assay,
                                            "_comb_celltype_bnd")]][[1]],
                        seurat_data$RNA_celltype_bnd)
  cm <- cm / rowSums(cm)
  
  pheatmap::pheatmap(cm)
  
  
  write.csv(seurat_res_bnd,
            file.path(save_dir, "files",
                      paste0(assay_type, "_celltype_mapping_bnd.csv")))
  
  
  #-------------------------------------------------------------------------------
  celltype_colors <- readRDS(file = file.path(save_dir, "color_palette.rds"))
  
  # Merge dfs
  
  seurat_res <- cbind(seurat_res_seurat, seurat_res_bnd)
  
  seurat_cluster <- cor_to_call(seurat_res) %>% 
    mutate(type = ifelse(r < cor_cutoff, "undetermined", type))
  
  new_clusters <- seurat_cluster$type
  names(new_clusters) <- seurat_cluster$cluster
  seurat_data[[paste0(seurat_assay, "_combined_celltype")]] <- 
    new_clusters[seurat_data[[paste0(seurat_assay, "_cluster")]][[1]]]
  
  print(plotDimRed(seurat_data,
                   col_by = paste0(seurat_assay, "_combined_celltype"),
                   plot_type = plot_type,
                   color = celltype_colors, ggrastr = TRUE))
  
  cm <- confusionMatrix(seurat_data[[paste0(seurat_assay, "_combined_celltype")]][[1]],
                        seurat_data$RNA_celltype)
  cm <- cm / rowSums(cm)
  
  pheatmap::pheatmap(cm)
  
  dev.off()
}

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))