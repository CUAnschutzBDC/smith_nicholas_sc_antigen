library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(muscat)
library(scran)
library(DESeq2)


source(here("src/scripts/muscat_plotting_functions.R"))

calc_logfc <- 1

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

#args <- commandArgs(trailingOnly = TRUE)

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

seurat_data$doublet_chains <- paste(seurat_data$Doublet_finder, seurat_data$n_chains, sep = "_")

plot1 <- featDistPlot(seurat_data, geneset = "nFeature_RNA", 
             sep_by = "Doublet_finder", combine = FALSE)
plot2 <- featDistPlot(seurat_data, geneset = "nCount_RNA", 
             sep_by = "Doublet_finder", combine = FALSE)

plot3 <- featDistPlot(seurat_data, geneset = "nFeature_RNA", 
             sep_by = "RNA_combined_celltype", 
             col_by = "Doublet_finder", combine = FALSE)
plot4 <- featDistPlot(seurat_data, geneset = "nCount_RNA", 
             sep_by = "RNA_combined_celltype", 
             col_by = "Doublet_finder", combine = FALSE)


plot5 <- featDistPlot(seurat_data, geneset = "n_chains",
             sep_by = "Doublet_finder", combine = FALSE)

plot6 <- featDistPlot(seurat_data, geneset = "nFeature_RNA",
             sep_by = "doublet_chains", combine = FALSE)
plot7 <- featDistPlot(seurat_data, geneset = "nCount_RNA",
             sep_by = "doublet_chains", combine = FALSE)


plot8 <- featDistPlot(seurat_data, geneset = "nFeature_RNA", 
             sep_by = "RNA_combined_celltype", 
             col_by = "doublet_chains", combine = FALSE)
plot9 <- featDistPlot(seurat_data, geneset = "nCount_RNA", 
             sep_by = "RNA_combined_celltype", 
             col_by = "doublet_chains", combine = FALSE)

pdf(file.path(save_dir, "images", "doublet_exploring.pdf"))
print(plot1)
print(plot2)
print(plot3)
print(plot4)
print(plot5)
print(plot6)
print(plot7)
print(plot8)
print(plot9)

cm <- confusionMatrix(seurat_data$Doublet_finder,
                      seurat_data$n_chains)
cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

cm <- confusionMatrix(seurat_data$Doublet_finder,
                      seurat_data$RNA_combined_celltype)
cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

cm <- confusionMatrix(seurat_data$RNA_combined_celltype,
                      seurat_data$Doublet_finder)
cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)




seurat_no_doublet <- subset(seurat_data, subset = Doublet_finder == "Singlet")

Idents(seurat_no_doublet) <- "RNA_combined_celltype"


top_markers <- FindAllMarkers(seurat_no_doublet, logfc.threshold = calc_logfc)

top_fifty_markers <- top_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 50, wt = avg_log2FC)


seurat_data$cluster_doublet <- paste(seurat_data$RNA_combined_celltype,
                                     seurat_data$Doublet_finder,
                                     sep = "_")


meta_data <- seurat_data[[]] %>%
  dplyr::select(RNA_combined_celltype, Doublet_finder, 
                cluster_doublet) %>%
  dplyr::distinct() %>%
  tibble::remove_rownames() 

rownames(meta_data) <- meta_data$cluster_doublet


celltype_colors <- readRDS(file.path(save_dir, "color_palette.rds"))
doublet_colors <- RColorBrewer::brewer.pal(n = 3, name = "Set1")
doublet_colors <- doublet_colors[1:2]
names(doublet_colors) <- c("Doublet", "Singlet")
  
all_colors <- list("RNA_combined_celltype" = celltype_colors,
                   "Doublet_finder" = doublet_colors)

grid::grid.newpage()
plot_heatmap(seurat_object = seurat_data,
             gene_list = unique(top_fifty_markers$gene), 
             meta_col = "cluster_doublet", colors = NULL, 
             meta_df = meta_data, color_list = all_colors, max_val = 2.5, 
             min_val = -2.5,  cluster_rows = TRUE, 
             cluster_cols = TRUE, average_expression = TRUE, 
             plot_meta_col = FALSE, plot_rownames = FALSE,
             cell_order = NULL, return_data = FALSE)



seurat_data$nchains <- seurat_data$n_chains
seurat_data$nchains[is.na(seurat_data$nchains)] <- 0
seurat_data$nchains <- factor(seurat_data$nchains)

seurat_data$celltype_nchains <- paste(seurat_data$RNA_combined_celltype,
                                     seurat_data$nchains,
                                     sep = "_")


meta_data <- seurat_data[[]] %>%
  dplyr::select(RNA_combined_celltype, nchains, 
                celltype_nchains) %>%
  dplyr::distinct() %>%
  tibble::remove_rownames() 

rownames(meta_data) <- meta_data$celltype_nchains


celltype_colors <- readRDS(file.path(save_dir, "color_palette.rds"))
n_chains_colors <- RColorBrewer::brewer.pal(n = 5, name = "Set1")
names(n_chains_colors) <- c("0", "1", "2", "3", "4")

all_colors <- list("RNA_combined_celltype" = celltype_colors,
                   "nchains" = n_chains_colors)
grid::grid.newpage()
plot_heatmap(seurat_object = seurat_data,
             gene_list = unique(top_fifty_markers$gene), 
             meta_col = "celltype_nchains", colors = NULL, 
             meta_df = meta_data, color_list = all_colors, max_val = 2.5, 
             min_val = -2.5,  cluster_rows = TRUE, 
             cluster_cols = TRUE, average_expression = TRUE, 
             plot_meta_col = FALSE, plot_rownames = FALSE,
             cell_order = NULL, return_data = FALSE)

# Repeat with no T cells
keep_celltypes <- c("Activated_memory", "B.intermediate", "BND2", "Memory_IgA",
                    "Naive_1", "Naive3", "Plasmablast", "Resting_memory")

seurat_b <- subset(seurat_data, subset = RNA_combined_celltype %in% keep_celltypes)

seurat_no_doublet <- subset(seurat_b, subset = Doublet_finder == "Singlet")

Idents(seurat_no_doublet) <- "RNA_combined_celltype"


top_markers <- FindAllMarkers(seurat_no_doublet, logfc.threshold = calc_logfc)

top_fifty_markers <- top_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 50, wt = avg_log2FC)


seurat_b$cluster_doublet <- paste(seurat_b$RNA_combined_celltype,
                                  seurat_b$Doublet_finder,
                                     sep = "_")

seurat_b$RNA_combined_celltype <- factor(seurat_b$RNA_combined_celltype)

meta_data <- seurat_b[[]] %>%
  dplyr::select(RNA_combined_celltype, Doublet_finder, 
                cluster_doublet) %>%
  dplyr::distinct() %>%
  tibble::remove_rownames() 

rownames(meta_data) <- meta_data$cluster_doublet


celltype_colors <- readRDS(file.path(save_dir, "color_palette.rds"))
doublet_colors <- RColorBrewer::brewer.pal(n = 3, name = "Set1")
doublet_colors <- doublet_colors[1:2]
names(doublet_colors) <- c("Doublet", "Singlet")

all_colors <- list("RNA_combined_celltype" = celltype_colors[names(celltype_colors) %in%
                                                                     keep_celltypes],
                   "Doublet_finder" = doublet_colors)

grid::grid.newpage()
plot_heatmap(seurat_object = seurat_b,
             gene_list = unique(top_fifty_markers$gene), 
             meta_col = "cluster_doublet", colors = NULL, 
             meta_df = meta_data, color_list = all_colors, max_val = 2.5, 
             min_val = -2.5,  cluster_rows = TRUE, 
             cluster_cols = TRUE, average_expression = TRUE, 
             plot_meta_col = FALSE, plot_rownames = FALSE,
             cell_order = NULL, return_data = FALSE)



seurat_b$nchains <- seurat_b$n_chains
seurat_b$nchains[is.na(seurat_b$nchains)] <- 0
seurat_b$nchains <- factor(seurat_b$nchains)

seurat_b$celltype_nchains <- paste(seurat_b$RNA_combined_celltype,
                                   seurat_b$nchains,
                                      sep = "_")


meta_data <- seurat_b[[]] %>%
  dplyr::select(RNA_combined_celltype, nchains, 
                celltype_nchains) %>%
  dplyr::distinct() %>%
  tibble::remove_rownames() 

rownames(meta_data) <- meta_data$celltype_nchains


celltype_colors <- readRDS(file.path(save_dir, "color_palette.rds"))
n_chains_colors <- RColorBrewer::brewer.pal(n = 5, name = "Set1")
names(n_chains_colors) <- c("0", "1", "2", "3", "4")

all_colors <- list("RNA_combined_celltype" = celltype_colors[names(celltype_colors) %in%
                                                                     keep_celltypes],
                   "nchains" = n_chains_colors)
grid::grid.newpage()
plot_heatmap(seurat_object = seurat_b,
             gene_list = unique(top_fifty_markers$gene), 
             meta_col = "celltype_nchains", colors = NULL, 
             meta_df = meta_data, color_list = all_colors, max_val = 2.5, 
             min_val = -2.5,  cluster_rows = TRUE, 
             cluster_cols = TRUE, average_expression = TRUE, 
             plot_meta_col = FALSE, plot_rownames = FALSE,
             cell_order = NULL, return_data = FALSE)


seurat_b_doublet <- subset(seurat_b, subset = Doublet_finder == "Doublet")

grid::grid.newpage()
plot_heatmap2(seurat_object = seurat_b_doublet,
             gene_list = unique(top_fifty_markers$gene), 
             meta_col = "RNA_combined_celltype", 
             colors = celltype_colors[names(celltype_colors) %in% keep_celltypes], 
             max_val = 2.5, min_val = -2.5,  cluster_rows = TRUE, 
             cluster_cols = TRUE, average_expression = FALSE, 
             plot_meta_col = TRUE, plot_rownames = FALSE,
             cell_order = NULL, return_data = FALSE)

dev.off()
