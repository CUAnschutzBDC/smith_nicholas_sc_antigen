library(here)
library(scAnalysisR)
library(pheatmap)
library(tidyverse)
library(MAST)
library(SingleCellExperiment)
library(Seurat)
library(cowplot)
library(SCOPfunctions)
library(KEGGREST)
library(org.Hs.eg.db)

source(here("src/scripts/muscat_plotting_functions.R"))

all_colors <- readRDS(file = file.path("files/all_colors.rds"))


final_colors <- all_colors$cell_type_colors

final_colors <- final_colors[names(final_colors) != "Translational_Intermediate"]

tetramer_colors <- all_colors$tetramer_colors

sample_colors <- all_colors$sample_colors


status_colors <- all_colors$status_colors

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
DefaultAssay(seurat_data) <- "RNA"

seurat_data <- DietSeurat(object = seurat_data, assays = "RNA")

seurat_data$Status <- factor(seurat_data$Status, levels = c("ND", "AAB", "T1D"))

# Function to run DE
run_de <- function(sisters, seurat_data, save_name){
  seurat_subset <- subset(seurat_data, subset = sample %in% sisters)
  
  print(unique(seurat_subset$sample))
  
  test_mapping <- c("DNA.tet" = "Other",
                    "GAD.tet" = "Islet_Reactive",
                    "IA2.tet" = "Islet_Reactive",
                    "INS.tet" = "Islet_Reactive",
                    "Islet_Multi_Reactive" = "Islet_Reactive",
                    "Negative" = "Other",
                    "Other_Multi_Reactive" = "Other",
                    "TET.tet" = "Other")
  
  seurat_subset$test_id <- test_mapping[seurat_subset$tet_name_cutoff]
  
  # Name columns for muscat sample_id, group_id, cluster_id
  seurat_subset$sample_id <- paste(seurat_subset$sample, seurat_subset$test_id,
                                 sep = "_")
  
  seurat_subset$test_full_new_status <- paste(seurat_subset$Status, 
                                            seurat_subset$test_id,
                                            sep = "_")
  
  
  Idents(seurat_subset) <- "test_full_new_status"
  
  tests_run <- list(c("T1D_Islet_Reactive", "ND_Islet_Reactive"))
  
  
  all_markers <- lapply(tests_run, function(x){
    markers <- FindMarkers(seurat_subset, test.use = "MAST", 
                           ident.1 = x[[1]],
                           ident.2 = x[[2]])
    
    markers$ident.1 <- x[[1]]
    markers$ident.2 <- x[[2]]
    markers$gene <- rownames(markers)
    
    return(markers)
    
  })
  
  all_markers <- do.call(rbind, all_markers)
  
  markers_sig <- all_markers[all_markers$p_val_adj < 0.05,]
  
  print(table(markers_sig$ident.1, markers_sig$ident.2))
  
  markers_sig$cluster <- paste(markers_sig$ident.1,
                               markers_sig$ident.2,
                               sep = "_")
  # save markers
  write.csv(markers_sig, file.path(save_dir, "files",
                                   paste0(save_name, "_mast_de.csv")))
  
}

#seurat_data$sample <- as.character(seurat_data$sample)

sisters <- c("109", "108", "107")

twins <- c("109", "108")
run_de(sisters = sisters, seurat_data = seurat_data, save_name = "sisters")
run_de(sisters = twins, seurat_data = seurat_data, save_name = "twins")
