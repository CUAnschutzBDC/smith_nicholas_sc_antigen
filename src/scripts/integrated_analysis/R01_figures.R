library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)
library(djvdj)
library(UpSetR)
library(muscat)
library(SingleCellExperiment)

source(here("src/scripts/muscat_plotting_functions.R"))

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "merged"

all_ref_dir <-
  "/beevol/home/wellskri/Analysis/references/single_cell_references"

HTO <- FALSE
ADT <- FALSE

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
base_dir <- here()

save_dir <- file.path(base_dir, "results", "R_analysis", sample)
fig_dir <- file.path(save_dir, "images", "figures")

ifelse(!dir.exists(fig_dir), dir.create(fig_dir), FALSE)

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))

plot_tetramers <- c("INS-tet", "TET-tet", "IA2-tet", "GAD-tet", "DNA-tet")

antigen_colors<- c("Negative" = "#7fc97f",
                   "INS-tet" = "#beaed4",
                   "TET-tet" = "#fdc086",
                   "IA2-tet" = "#ffff99",
                   "GAD-tet" = "#386cb0",
                   "DNA-tet" = "#f0027f")

b_cells <- c("Activated_memory", "B.intermediate", "BND2", "Memory_IgA",
             "Naive_1", "Naive_3", "Plasmablast", "Resting_memory")

seurat_data <- subset(seurat_data, subset = RNA_combined_celltype %in% b_cells)

`%notin%` <- Negate(`%in%`)

seurat_data_tet <- subset(seurat_data, subset = tet_hash_id %notin% 
                        c("Negative", "Doublet"))

seurat_data_scar <- subset(seurat_data, subset = scar_hash_id %notin%
                             c("Negative", "Doublet"))

p1 <- scAnalysisR::stacked_barplots(seurat_data_tet, meta_col = "tet_hash_id",
                              color = antigen_colors,
                              split_by = "RNA_combined_celltype") + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))


p2 <- scAnalysisR::stacked_barplots(seurat_data_scar, meta_col = "scar_hash_id",
                              color = antigen_colors,
                              split_by = "RNA_combined_celltype") + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))


p3 <- scAnalysisR::stacked_barplots(seurat_data_tet, meta_col = "tet_hash_id",
                              color = antigen_colors,
                              split_by = "Status") + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))


p4 <- scAnalysisR::stacked_barplots(seurat_data_scar, meta_col = "scar_hash_id",
                              color = antigen_colors,
                              split_by = "Status") + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))


pdf(file.path(save_dir, "images", "tetramer_barplots.pdf"),
    height = 6, width = 8)

print(p2)
print(p4)

dev.off()
