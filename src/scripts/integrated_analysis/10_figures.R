library(here)
library(scAnalysisR)
library(pheatmap)
library(tidyverse)
library(splitstackshape)
library(circlize)
library(viridis)
library(ggalluvial)
library(Seurat)

normalization_method <- "log" # can be SCT or log
# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

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

# Set up factors ---------------------------------------------------------------
seurat_data$Status <- factor(seurat_data$Status, levels = c("ND", "AAB", "T1D"))

plotting_levels <- c("Negative", "DNA-tet", "TET-tet",
                     "Islet_Multi_Reactive", "IA2-tet", 
                     "GAD-tet", "INS-tet")

seurat_data$tet_hash_id <- factor(seurat_data$tet_hash_id, 
                                  levels = plotting_levels)

image_dir <- file.path(save_dir, "images", "final_figures")

ifelse(!dir.exists(image_dir), dir.create(image_dir), FALSE)

# Colors -----------------------------------------------------------------------
all_colors <- readRDS(file = file.path("files/all_colors.rds"))


final_colors <- all_colors$cell_type_colors

tetramer_colors <- all_colors$tetramer_colors

sample_colors <- all_colors$sample_colors


status_colors <- all_colors$status_colors

# Plots ------------------------------------------------------------------------

## QC plots --------------------------------------------------------------------
# 1.	Guidance on which QC plots are necessary to show
# Violin plots of mito, features, genes across samples
quality_plot <- featDistPlot(seurat_data, geneset = c("nCount_RNA",
                                                      "nFeature_RNA", 
                                                      "percent.mt"),
                             sep_by = "sample", col_by = "Status",
                             color = status_colors,
                             combine = FALSE)

quality_plot$percent.mt <- quality_plot$percent.mt +
  ggplot2::ylim(0, 20)

final_quality <- cowplot::plot_grid(plotlist = quality_plot, nrow = 3)

pdf(file.path(image_dir, "quality_plots.pdf"))
print(final_quality)
dev.off()

# UMAP of samples before and after correction

## Cell type information -------------------------------------------------------
# 2.	Updated UMAP of all cells clustered according to mRNA for cell type
# (Naïve, Transitional (intermediate), Resting Memory, Activated Memory, Plasmablast)

pdf(file.path(image_dir, "rna_umap"), width = 8, height = 8)

print(plotDimRed(seurat_data, col_by = "RNA_combined_celltype", 
                 plot_type = "rna_mnn.umap",
           color = final_colors, ggrastr = TRUE))

dev.off()

# 3.	List of top 10 genes contributing to cell type determination
# * Not these are not really genes contributing to cell type determination as 
# this is not how cell type was determined (this is important for you when 
# writing the text), these are just the top 10 marker genes of each cluster.

Idents(seurat_data) <- seurat_data$RNA_combined_celltype

all_markers <- FindAllMarkers(seurat_data, assay = "RNA", 
                              logfc.threshold = 0.25, max.cells.per.ident = 2000,
                              only.pos = TRUE)

# Ran this to make sure the results didn't change with the downsampling,
# they did not.
# all_markers_short <- all_markers %>%
#   dplyr::filter(abs(avg_log2FC) > 1)

use_markers <- all_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 10, wt = avg_log2FC)

plot_heatmap(seurat_data, gene_list = use_markers$gene,
             colors = final_colors, meta_col = "RNA_combined_celltype",
             average_expression = TRUE, assay = "RNA")

plot_heatmap(seurat_data, gene_list = use_markers$gene,
             colors = final_colors, meta_col = "RNA_combined_celltype",
             average_expression = FALSE, assay = "RNA")

# 4.	UMAP - Antigen reactive cells cluster independently (not clustering by mRNA)
# INS/GAD/IA2/Multi-islet/TET/DNA
# * We decided thta this isn't a necessary plot
# 5.	Stacked bar chart: x axis = B cell subtype (updated), y axis = % of 
# total cells minus other_multi_reactive
# (INS/GAD/IA2/Multi-islet/TET/DNA/Negative)

no_other <- subset(seurat_data, subset = tet_hash_id != "Other_Multi_Reactive")

`%notin%` <- Negate(`%in%`)

tetramer_colors_use <- tetramer_colors[names(tetramer_colors) %notin%
                                         c("Other_Multi_Reactive", 
                                           "Islet_Reactive")]

stacked_barplots(no_other, meta_col = "tet_hash_id",
                 split_by = "RNA_combined_celltype",
                 color = tetramer_colors_use) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1))

# 6.	INSR mRNA expression violin plot separated by B cell subtype
# 7.	INSR mRNA expression violin plot separated by status
# 8.	IGHM mRNA expression violin plot separated by status
# 9.	IGHD mRNA expression violin plot separated by status
# 10.	IGHG mRNA expression violin plot separated by status --> IGHG1, 2, 3, 4
# 11.	IGHA mRNA expression violin plot separated by status --> IGHA1, 2
# *make y axis of plots 7-11 consistent to show difference in values*

# Also make IgM and IgD with protein

violin_1 <- featDistPlot(seurat_data, geneset = "INSR",
                         col_by = "RNA_combined_celltype",
                         sep_by = "RNA_combined_celltype",
                         color = final_colors, combine = FALSE)

all_violins <- featDistPlot(seurat_data, geneset = c("INSR", "IGHM", "IGHD", 
                                                     "IGHG1", "IGHA1"),
                            col_by = "Status", sep_by = "RNA_combined_celltype",
                            color = status_colors, combine = FALSE)

# seurat_data <- subset(seurat_data, subset = isotype != "IGHE")
# 
# all_violins <- featDistPlot(seurat_data, geneset = c("IGHG1", "IGHA1"),
#                             col_by = "isotype", sep_by = "isotype",
#                             color = status_colors, combine = FALSE)

# 8

violin_1 <- violin_1$INSR +
  ggplot2::ylim(0, 8)

all_violins <- lapply(all_violins, function(x){
  x <- x + 
    ggplot2::ylim(0, 8)
  
  return(x)
})


graphics.off()

pdf(file.path(image_dir, "gene_expression_violins_status.pdf"), height = 4,
    width = 8)
print(all_violins)

dev.off()


pdf(file.path(image_dir, "gene_expression_violins_celltype.pdf"), height = 5,
    width = 8)
print(violin_1)

dev.off()

## Differential expression ------------------------------------------------------
# 12.	New DEG comparisons heatmaps (first pass plots, will narrow down later)
# 13.	BCR signaling pathway (KEGG) genes heatmap

# Heatmap by sample/status + antigen binding

# Add in twin DE

## Isotype analysis -------------------------------------------------------------
# 14.	Stacked bar chart of all cells: y axis = isotype frequency, 
# x axis = disease status 

all_cells <- stacked_barplots(seurat_data, meta_col = "isotype",
                              split_by = "Status")

pdf(file.path(image_dir, "isotype_barplot_all.pdf"))

print(all_cells)

dev.off()

# 15.	Stacked bar chart of islet-reactive cells: y axis = isotype frequency, 
# x axis = disease status 
islet_reactive <- subset(seurat_data, 
                         subset = tet_hash_id %in% c("Islet_Multi_Reactive",
                                                     "IA2-tet", "GAD-tet",
                                                     "INS-tet"))


islet_cells <- stacked_barplots(islet_reactive, meta_col = "isotype",
                                split_by = "Status")

pdf(file.path(image_dir, "isotype_barplot_islet_reactive.pdf"))

print(islet_cells)

dev.off()

## SHM -------------------------------------------------------------------------
# 16.	SHM histogram - all cells
# 17.	SHM histogram - islet-reactive cells

## Diversity -------------------------------------------------------------------
# 18.	Shannon alpha diversity rarefaction curve of all cells by status
# 19.	Jaccard beta diversity rarefaction curve of islet-reactive cells by status


## Clonal analysis -------------------------------------------------------------
# 20.	Count number of islet-reactive clones found within each individual, 
# make a stacked bar chart 
# delineating whether the clones were public (black) or private (white), and 
# group the individuals on 
# the X axis according to disease status. Y axis is the number of clones present
# according to category. 
# 21.	Do the same as 20, but instead of raw # of clones, normalize to a percent
# of total clones found in 
# each individual, so if 4 clones total and 1 was public, 0.25 of the bar is 
# black and 0.75 of the bar is white.

## Antigen reactvitiy ----------------------------------------------------------
# 22.	Stacked bar chart of frequency of antigen reactivity type where 
# x axis = disease status and y 
# axis = INS/GAD/IA2/Multi-Islet TET/DNA
# 23.	Same as 22 but collapse all islet-reactive to one category so 
# y axis = All-Islet TET/DNA
# 24.	Same as 22 & 23 but include the negative fraction of cells
# 25.	Same as 22 & 23 but include the other-reactive fraction of cells
# 26.	Same as 22 & 23 but include the other-reactive and negative 
# fractions of cells

## Other -----------------------------------------------------------------------
# 27.	Replicate something similar to 
# [Figure 6 of this paper](https://journals.aai.org/jimmunol/article/198/4/1460/109668/Dysregulation-of-B-Cell-Repertoire-Formation-in) 
# for our samples with all cells (by status)
# 28.	Replicate something similar to Figure 6 of this paper for our samples with islet-reactive cells (by status)
# 29.	Replicate something similar to 
# [Figure 6 of this paper](https://journals.aai.org/jimmunol/article/198/4/1460/109668/Dysregulation-of-B-Cell-Repertoire-Formation-in) 
# for our samples with all clones (by status)
# 30.	Replicate something similar to 
# [Figure 6 of this paper](https://journals.aai.org/jimmunol/article/198/4/1460/109668/Dysregulation-of-B-Cell-Repertoire-Formation-in) 
# for our samples with public clones (by status)
# 31.	Replicate something similar to 
# [Figure 6 of this paper](https://journals.aai.org/jimmunol/article/198/4/1460/109668/Dysregulation-of-B-Cell-Repertoire-Formation-in) 
# for our samples with private clones (by status)