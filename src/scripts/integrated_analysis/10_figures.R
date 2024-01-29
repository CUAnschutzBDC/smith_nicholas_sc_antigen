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

seurat_data <- subset(seurat_data, subset = imcantation_isotype != "IGHE")

# Set up factors ---------------------------------------------------------------
seurat_data$Status <- factor(seurat_data$Status, levels = c("ND", "AAB", "T1D"))

plotting_levels <- c("Negative", "DNA-tet", "TET-tet", "Other_Multi_Reactive",
                     "Islet_Multi_Reactive", "IA2-tet", 
                     "GAD-tet", "INS-tet")

seurat_data$libra_tet_hash_id <- factor(seurat_data$libra_tet_hash_id, 
                                  levels = plotting_levels)

seurat_data$final_celltype <- factor(seurat_data$final_celltype,
                                     levels = c("Naive", "ABC", "Resting_Memory",
                                                "Memory", "Plasmablast"))

image_dir <- file.path(save_dir, "images", "final_figures")

ifelse(!dir.exists(image_dir), dir.create(image_dir), FALSE)

# Colors -----------------------------------------------------------------------
all_colors <- readRDS(file = file.path("files/all_colors.rds"))


final_colors <- all_colors$cell_type_colors

final_colors <- final_colors[names(final_colors) != "Translational_Intermediate"]

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

umap_one <- plotDimRed(seurat_data, col_by = "sample",
                       plot_type = "pca.umap", color = sample_colors,
                       ggrastr = TRUE)[[1]]

umap_two <- plotDimRed(seurat_data, col_by = "final_celltype",
                       plot_type = "pca.umap", color = final_colors,
                       ggrastr = TRUE)[[1]]


umap_three <- plotDimRed(seurat_data, col_by = "sample",
                         plot_type = "rna_mnn.umap", color = sample_colors,
                         ggrastr = TRUE)[[1]]

umap_four <- plotDimRed(seurat_data, col_by = "final_celltype",
                        plot_type = "rna_mnn.umap", color = final_colors,
                        ggrastr = TRUE)[[1]]

full_plot <- cowplot::plot_grid(umap_one, umap_two, umap_three, umap_four,
                                rel_widths = c(0.84, 1, 0.84, 1),
                                nrow = 2, ncol = 2)

pdf(file.path(image_dir, "quality_umaps.pdf"),
    width = 12, height = 10)

print(full_plot)

dev.off()

## Cell type information -------------------------------------------------------
# 2.	Updated UMAP of all cells clustered according to mRNA for cell type
# (Naïve, Transitional (intermediate), Resting Memory, Activated Memory, Plasmablast)

seurat_data$celltype_cluster <- paste(seurat_data$final_celltype, 
                                      seurat_data$RNA_cluster,
                                      sep = "_")



        

cluster_celltype_colors <- c("Naive_0" = "#58b44d",
                             "Naive_1" = "#9daf5d",
                             "Resting_Memory_2" = "#746ec8",
                             "Naive_3" = "#61732b",
                             "Memory_4" = "#9f7239",
                             "Naive_5" = "#b2b72e",
                             "Naive_6" = "#61a0d5",
                             "ABC_7" = "#4bae8b",
                             "Resting_Memory_8" = "#bd5abe",
                             "Resting_Memory_9" = "#c3618a",
                             "Plasmablast_10" = "#d89a40",
                             "Plasmablast_11" = "#bf673d")

seurat_data$celltype_cluster <- factor(seurat_data$celltype_cluster,
                                       levels = c("Naive_0", "Naive_1", 
                                                  "Naive_3", "Naive_5", 
                                                  "Naive_6", 
                                                  "Resting_Memory_2",
                                                  "Resting_Memory_8", 
                                                  "Resting_Memory_9",
                                                  "ABC_7",
                                                  "Memory_4",
                                                  "Plasmablast_10",
                                                  "Plasmablast_11"))

cluster_celltype_colors <- cluster_celltype_colors[order(match(names(cluster_celltype_colors),
                                                               levels(seurat_data$celltype_cluster)))]

pdf(file.path(image_dir, "rna_umap.pdf"), width = 8, height = 8)

print(plotDimRed(seurat_data, col_by = "celltype_cluster", 
                 plot_type = "rna_mnn.umap",
           color = cluster_celltype_colors, ggrastr = TRUE))

dev.off()

# 3.	List of top 10 genes contributing to cell type determination
# * Not these are not really genes contributing to cell type determination as 
# this is not how cell type was determined (this is important for you when 
# writing the text), these are just the top 10 marker genes of each cluster.

Idents(seurat_data) <- seurat_data$celltype_cluster

all_markers <- FindAllMarkers(seurat_data, assay = "RNA", 
                              logfc.threshold = 0.25, max.cells.per.ident = 2000,
                              only.pos = TRUE)

save_excel <- openxlsx::createWorkbook()

for(clust in unique(all_markers$cluster)){
  write_markers <- all_markers %>%
    dplyr::filter(cluster == clust, p_val_adj < 0.05)
  
  openxlsx::addWorksheet(wb = save_excel, sheetName = clust)
  
  openxlsx::writeData(wb = save_excel, sheet = clust, x = write_markers)
}

openxlsx::saveWorkbook(wb = save_excel, 
                       file = file.path(image_dir, "cluster_markers.xslx"),
                       overwrite = TRUE)

# Ran this to make sure the results didn't change with the downsampling,
# they did not.
# all_markers_short <- all_markers %>%
#   dplyr::filter(abs(avg_log2FC) > 1)

use_markers <- all_markers %>%
  dplyr::filter(!grepl("IGK|IGH|IGL", gene)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 10, wt = avg_log2FC) %>%
  dplyr::arrange(cluster, desc(avg_log2FC))

graphics.off()

pdf(file.path(image_dir, "marker_heatmap_average.pdf"),
    width = 8, height = 12)

plot_heatmap(seurat_data, gene_list = use_markers$gene,
             colors = cluster_celltype_colors, meta_col = "celltype_cluster",
             average_expression = TRUE, assay = "RNA")

dev.off()

pdf(file.path(image_dir, "marker_heatmap_all.pdf"),
    width = 8, height = 12)

plot_heatmap(seurat_data, gene_list = use_markers$gene,
             colors = cluster_celltype_colors, meta_col = "celltype_cluster",
             average_expression = FALSE, assay = "RNA")

dev.off()

graphics.off()

# Cell type by individual
pdf(file.path(image_dir, "sample_celltype_barplot.pdf"),
    width = 8, height = 8)
barplot <- stacked_barplots(seurat_data, meta_col = "celltype_cluster",
                            split_by = "sample",
                            color = cluster_celltype_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1))

print(barplot)
dev.off()

# Cell type by status
pdf(file.path(image_dir, "status_celltype_barplot.pdf"),
    width = 8, height = 8)
barplot <- stacked_barplots(seurat_data, meta_col = "celltype_cluster",
                            split_by = "Status",
                            color = cluster_celltype_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1))

print(barplot)

dev.off()

# 4.	UMAP - Antigen reactive cells cluster independently (not clustering by mRNA)
# INS/GAD/IA2/Multi-islet/TET/DNA
# * We decided thta this isn't a necessary plot
# 5.	Stacked bar chart: x axis = B cell subtype (updated), y axis = % of 
# total cells minus other_multi_reactive
# (INS/GAD/IA2/Multi-islet/TET/DNA/Negative)

no_other <- subset(seurat_data, subset = libra_tet_hash_id != "Other_Multi_Reactive")

`%notin%` <- Negate(`%in%`)

tetramer_colors_use <- tetramer_colors[names(tetramer_colors) %notin%
                                         c("Other_Multi_Reactive", 
                                           "Islet_Reactive")]

barplot <- stacked_barplots(no_other, meta_col = "libra_tet_hash_id",
                            split_by = "final_celltype",
                            color = tetramer_colors_use) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1))

pdf(file.path(image_dir, "tetramer_stacked_barplot.pdf"),
    height = 8, width = 8)
print(barplot)
dev.off()

# 6.	INSR mRNA expression violin plot separated by B cell subtype
# 7.	INSR mRNA expression violin plot separated by status
# 8.	IGHM mRNA expression violin plot separated by status
# 9.	IGHD mRNA expression violin plot separated by status
# 10.	IGHG mRNA expression violin plot separated by status --> IGHG1, 2, 3, 4
# 11.	IGHA mRNA expression violin plot separated by status --> IGHA1, 2
# *make y axis of plots 7-11 consistent to show difference in values*

# Also make IgM and IgD with protein
# Make all G and A
# Plot INSR by tet binding

violin_1 <- featDistPlot(seurat_data, geneset = "INSR",
                         col_by = "final_celltype",
                         sep_by = "final_celltype",
                         color = final_colors, combine = FALSE)

violin_2 <- featDistPlot(seurat_data, geneset = "INSR",
                         col_by = "libra_tet_hash_id",
                         sep_by = "libra_tet_hash_id",
                         color = tetramer_colors, combine = FALSE)

all_violins <- featDistPlot(seurat_data, geneset = c("INSR", "IGHM", "IGHD", 
                                                     "IGHG1", "IGHG2",
                                                     "IGHG3", "IGHG3", 
                                                     "IGHA1", "IGHA2"),
                            col_by = "Status", sep_by = "final_celltype",
                            color = status_colors, combine = FALSE)

# seurat_data <- subset(seurat_data, subset = isotype != "IGHE")
# 
# all_violins <- featDistPlot(seurat_data, geneset = c("IGHG1", "IGHA1"),
#                             col_by = "isotype", sep_by = "isotype",
#                             color = status_colors, combine = FALSE)

# 8

violin_1 <- violin_1$INSR +
  ggplot2::ylim(0, 8)

violin_2 <- violin_2$INSR +
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
print(violin_2)

dev.off()

## Differential expression ------------------------------------------------------
# 12.	New DEG comparisons heatmaps (first pass plots, will narrow down later)
# 13.	BCR signaling pathway (KEGG) genes heatmap

# Heatmap by sample/status + antigen binding

# Add in twin DE

## Isotype analysis -------------------------------------------------------------
# 14.	Stacked bar chart of all cells: y axis = isotype frequency, 
# x axis = disease status 

isotype <- subset(seurat_data, subset = imcantation_isotype != "")

colors <- MetBrewer::met.brewer(name = "Peru1", n = 6)
colors2 <- MetBrewer::met.brewer(name = "Peru2", n = 7)
isotype_colors <- c(colors[c(5, 6, 3, 4)], colors2[4:7])
names(isotype_colors) <- c("IGHM", "IGHD",  "IGHA1", "IGHA2",
                           "IGHG1", "IGHG2", "IGHG3", "IGHG4")
isotype$imcantation_isotype <- factor(isotype$imcantation_isotype,
                                      levels = names(isotype_colors))

all_cells <- stacked_barplots(isotype, meta_col = "imcantation_isotype",
                              split_by = "Status", color = isotype_colors)

pdf(file.path(image_dir, "isotype_barplot_all.pdf"), height = 8, width = 8)

print(all_cells)

dev.off()

# 15.	Stacked bar chart of islet-reactive cells: y axis = isotype frequency, 
# x axis = disease status 
islet_reactive <- subset(isotype, 
                         subset = libra_tet_hash_id %in% c("Islet_Multi_Reactive",
                                                     "IA2-tet", "GAD-tet",
                                                     "INS-tet"))


islet_cells <- stacked_barplots(islet_reactive, meta_col = "imcantation_isotype",
                                split_by = "Status", color = isotype_colors)

pdf(file.path(image_dir, "isotype_barplot_islet_reactive.pdf"),
    height = 8, width = 8)

print(islet_cells)

dev.off()

# TODO
# Sepaarate by cell type as well
# Plot next to insulin score

# isotype_violin <- featDistPlot(isotype, geneset = "INS-tet",
#                                assay = "SCAR_TET_LOG",
#                                sep_by = "Status", col_by = "imcantation_isotype",
#                                color = isotype_colors,
#                                combine = FALSE)


isotype_violin <- featDistPlot(isotype, geneset = c("INS-tet",
                                                    "GAD-tet",
                                                    "IA2-tet",
                                                    "TET-tet",
                                                    "DNA-tet"),
                               assay = "SCAR_TET_LOG",
                               col_by = "Status", sep_by = "imcantation_isotype",
                               color = status_colors,
                               combine = FALSE)


isotype_violin2 <- featDistPlot(islet_reactive, geneset = c("INS-tet",
                                                    "GAD-tet",
                                                    "IA2-tet",
                                                    "TET-tet",
                                                    "DNA-tet"),
                               assay = "SCAR_TET_LIBRA",
                               col_by = "Status", sep_by = "imcantation_isotype",
                               color = status_colors,
                               combine = FALSE)
pdf(file.path(image_dir, "isotype_status_tet_log.pdf"),
    height = 8, width = 8)

print(isotype_violin)

dev.off()

isotype_violin <- featDistPlot(isotype, geneset = c("INS-tet",
                                                    "GAD-tet",
                                                    "IA2-tet",
                                                    "TET-tet",
                                                    "DNA-tet"),
                               assay = "SCAR_TET_PROPORTIONS",
                               col_by = "Status", sep_by = "imcantation_isotype",
                               color = status_colors,
                               combine = FALSE)[[1]]+
  ggplot2::ylim(-1, 50)


pdf(file.path(image_dir, "isotype_status_tet_norm.pdf"),
    height = 8, width = 8)

print(isotype_violin)

dev.off()

# Make for each reactivity only


## SHM -------------------------------------------------------------------------
# 16.	SHM histogram - all cells
# 17.	SHM histogram - islet-reactive cells

sep_columns <- c("chains", "cdr3", "cdr3_length",
                 "cdr3_nt_length", "v_gene", "d_gene", "j_gene", "c_gene",
                 "reads", "umis", "productive", "full_length",
                 "v_ins", "v_del", "v_mis", "d_ins", "d_del",
                 "d_mis", "j_ins", "j_del", "j_mis", "c_ins", "c_del", "c_mis",
                 "all_ins", "all_del", "all_mis", "vd_ins", "vd_del", "dj_ins",
                 "dj_del", "v_mis_freq", "d_mis_freq", "j_mis_freq",
                 "c_mis_freq", "all_mis_freq")
keep_columns <- c("isotype", "final_celltype", "sample", "paired",
                  "clonotype_id", "Status", "libra_tet_hash_id", "full_libra_tet_hash_id",
                  "all_chains")

all_info <- seurat_data[[]] %>%
  dplyr::mutate(all_chains = chains) %>%
  dplyr::select(dplyr::all_of(c(sep_columns, keep_columns)), n_chains) %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::filter(!is.na(n_chains))


all_info_split <- cSplit(all_info, sep_columns, sep = ";", direction = "long") %>%
  dplyr::filter(!is.na(chains))


# CDR3 length

heavy_data <- all_info_split %>%
  dplyr::filter(chains == "IGH")


# ridge_cdr3_len <- ggplot2::ggplot(heavy_data, ggplot2::aes(x = cdr3_length,
#                                                            y = Status,
#                                                            fill = Status)) +
#   ggridges::geom_density_ridges() +
#   ggplot2::scale_fill_manual(values = status_colors)

violin_cdr3_len <- ggplot2::ggplot(heavy_data, ggplot2::aes(y = cdr3_length,
                                                           x = Status,
                                                           fill = Status)) +
  ggplot2::geom_violin(adjust = 1.5) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::stat_summary(fun = median, geom = "point", size = 2,
                        position = ggplot2::position_dodge(1))
  
pdf(file.path(image_dir, "heavy_cdr3_len_violin.pdf"),
    height = 6, width = 8)

print(violin_cdr3_len)

dev.off()

# all mis
# ridge_mis <- ggplot2::ggplot(heavy_data, ggplot2::aes(x = all_mis,
#                                                            y = Status,
#                                                            fill = Status)) +
#   ggridges::geom_density_ridges() +
#   ggplot2::scale_fill_manual(values = status_colors)

violin_mis <- ggplot2::ggplot(heavy_data, ggplot2::aes(y = all_mis,
                                                            x = Status,
                                                            fill = Status)) +
  ggplot2::geom_violin() +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::stat_summary(fun = median, geom = "point", size = 2,
                        position = ggplot2::position_dodge(1)) 

violin_mis_zoom <- ggplot2::ggplot(heavy_data, ggplot2::aes(y = all_mis,
                                                            x = Status,
                                                            fill = Status)) +
  ggplot2::geom_violin(adjust = 1.5) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::stat_summary(fun = median, geom = "point", size = 2,
                        position = ggplot2::position_dodge(1)) +
  ggplot2::ylim(-1, 50)


pdf(file.path(image_dir, "heavy_cdr3_smh_violin.pdf"),
    height = 6, width = 8)

print(violin_mis)
print(violin_mis_zoom)

dev.off()

# all miss freq

# ridge_mis_freq <- ggplot2::ggplot(heavy_data, ggplot2::aes(x = all_mis_freq,
#                                                            y = Status,
#                                                            fill = Status)) +
#   ggridges::geom_density_ridges() +
#   ggplot2::scale_fill_manual(values = status_colors) +
#   ggplot2::xlim(c(-0.01, 0.4))

violin_mis_freq <- ggplot2::ggplot(heavy_data, ggplot2::aes(y = all_mis_freq,
                                                            x = Status,
                                                            fill = Status)) +
  ggplot2::geom_violin(adjust = 2.5) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::stat_summary(fun = median, geom = "point", size = 2,
                        position = ggplot2::position_dodge(1)) 

violin_mis_freq_zoom <- ggplot2::ggplot(heavy_data, ggplot2::aes(y = all_mis_freq,
                                                            x = Status,
                                                            fill = Status)) +
  ggplot2::geom_violin(adjust = 2.5) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::stat_summary(fun = median, geom = "point", size = 2,
                        position = ggplot2::position_dodge(1)) +
  ggplot2::ylim(-0.01, 0.2)

pdf(file.path(image_dir, "heavy_cdr3_smh_freq_violin.pdf"),
    height = 6, width = 8)

print(violin_mis_freq)
print(violin_mis_freq_zoom)

dev.off()

# ridge_mis_freq_cell <- ggplot2::ggplot(heavy_data, ggplot2::aes(x = all_mis_freq,
#                                                                 y = final_celltype,
#                                                                 fill = final_celltype)) +
#   ggridges::geom_density_ridges() +
#   ggplot2::scale_fill_manual(values = final_colors) +
#   ggplot2::xlim(c(-0.01, 0.4))


violin_cell <- ggplot2::ggplot(heavy_data, ggplot2::aes(y = all_mis_freq,
                                                            x = final_celltype,
                                                            fill = Status)) +
  ggplot2::geom_violin(adjust = 2.5) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::stat_summary(fun = median, geom = "point", size = 2,
                        position = ggplot2::position_dodge(1)) +
  ggplot2::ylim(-0.01, 0.2)

violin_cell <- ggplot2::ggplot(heavy_data, ggplot2::aes(y = all_mis_freq,
                                                        x = final_celltype,
                                                        fill = final_celltype)) +
  ggplot2::geom_violin(adjust = 2.5) +
  ggplot2::scale_fill_manual(values = final_colors) +
  ggplot2::stat_summary(fun = median, geom = "point", size = 2,
                        position = ggplot2::position_dodge(1)) +
  ggplot2::ylim(-0.01, 0.2)


## Diversity -------------------------------------------------------------------
# 18.	Shannon alpha diversity rarefaction curve of all cells by status
# 19.	Jaccard beta diversity rarefaction curve of islet-reactive cells by status

seurat_data <- calc_diversity(seurat_data,
  data_col    = "final_clone",
  cluster_col = "Status",
  method      = abdiv::shannon
)

plot_diversity(seurat_data,
  data_col    = "final_clone",
  cluster_col = "Status",
  method      = abdiv::shannon,
  plot_colors = as.character(status_colors)
)

# Need to install iNEXT into docker image
plot_rarefaction(seurat_data, 
  data_col    = "final_clone",
  cluster_col = "Status",
  method      = c("shannon"),
  plot_colors = as.character(status_colors)
)
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