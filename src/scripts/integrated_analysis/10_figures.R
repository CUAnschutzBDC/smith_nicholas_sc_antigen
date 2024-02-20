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
library(KEGGREST)
library(org.Hs.eg.db)


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

seurat_data <- subset(seurat_data, subset = imcantation_isotype != "IGHE")

# Set up factors ---------------------------------------------------------------
seurat_data$Status <- factor(seurat_data$Status, levels = c("ND", "AAB", "T1D"))

plotting_levels <- c("Negative", "DNA.tet", "TET.tet", "Other_Multi_Reactive",
                     "Islet_Multi_Reactive", "IA2.tet", 
                     "GAD.tet", "INS.tet")

seurat_data$tet_name_cutoff <- factor(seurat_data$tet_name_cutoff, 
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

## Figure 1 --------------------------------------------------------------------
## Cell type information -------------------------------------------------------
# 2.	Updated UMAP of all cells clustered according to mRNA for cell type
# (Naïve, Transitional (intermediate), Resting Memory, Activated Memory, Plasmablast)

# Number of cells by status/tetramer
# With only singlets

# With multis included

seurat_data$celltype_cluster <- paste(seurat_data$final_celltype, 
                                      seurat_data$RNA_cluster,
                                      sep = "_")


new_clusters <- c("Naive_0" = "Naive_0",
                  "Naive_1" = "Naive_0",
                  "Resting_Memory_2" = "Resting_Memory_1",
                  "Naive_3" = "Naive_2",
                  "Memory_4" = "Memory_3",
                  "Naive_5" = "Naive_4",
                  "Naive_6" = "Naive_0",
                  "ABC_7" = "ABC_5",
                  "Resting_Memory_8" = "Resting_Memory_6",
                  "Resting_Memory_9" = "Resting_Memory_7",
                  "Plasmablast_10" = "Plasmablast_8",
                  "Plasmablast_11" = "Plasmablast_8")

seurat_data$celltype_cluster <- new_clusters[seurat_data$celltype_cluster]




cluster_celltype_colors <- c("Naive_0" = "#58b44d",
                             "Resting_Memory_1" = "#746ec8",
                             "Naive_2" = "#61732b",
                             "Memory_3" = "#9f7239",
                             "Naive_4" = "#b2b72e",
                             "ABC_5" = "#4bae8b",
                             "Resting_Memory_6" = "#bd5abe",
                             "Resting_Memory_7" = "#61a0d5",
                             "Plasmablast_8" = "#d89a40")

seurat_data$celltype_cluster <- factor(seurat_data$celltype_cluster,
                                       levels = c("Naive_0", "Naive_2", 
                                                  "Naive_4", 
                                                  "Resting_Memory_1",
                                                  "Resting_Memory_6", 
                                                  "Resting_Memory_7",
                                                  "ABC_5",
                                                  "Memory_3",
                                                  "Plasmablast_8"))

cluster_celltype_colors <- cluster_celltype_colors[order(match(names(cluster_celltype_colors),
                                                               levels(seurat_data$celltype_cluster)))]

pdf(file.path(image_dir, "1D_rna_umap.pdf"), width = 8, height = 8)

print(plotDimRed(seurat_data, col_by = "celltype_cluster", 
                 plot_type = "rna_mnn.umap",
                 color = cluster_celltype_colors, ggrastr = TRUE))

dev.off()

#	List of top 10 genes per cell type
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
                       file = file.path(image_dir, "cluster_markers.xlsx"),
                       overwrite = TRUE)

# Ran this to make sure the results didn't change with the downsampling,
# they did not.
# all_markers_short <- all_markers %>%
#   dplyr::filter(abs(avg_log2FC) > 1)

use_markers <- all_markers %>%
  dplyr::filter(!grepl("IGK|IGH|IGL|AL139020.1", gene)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 10, wt = avg_log2FC) %>%
  dplyr::arrange(cluster, desc(avg_log2FC))

graphics.off()

pdf(file.path(image_dir, "1E_marker_heatmap_average.pdf"),
    width = 8, height = 12)

plot_heatmap(seurat_data, gene_list = use_markers$gene,
             colors = cluster_celltype_colors, meta_col = "celltype_cluster",
             average_expression = TRUE, assay = "RNA")

dev.off()

# pdf(file.path(image_dir, "marker_heatmap_all.pdf"),
#     width = 8, height = 12)
# 
# plot_heatmap(seurat_data, gene_list = use_markers$gene,
#              colors = cluster_celltype_colors, meta_col = "celltype_cluster",
#              average_expression = FALSE, assay = "RNA")
# 
# dev.off()

graphics.off()

sample_order <- seurat_data[[]] %>%
  dplyr::select(sample, Status) %>%
  dplyr::distinct() %>%
  dplyr::arrange(Status)

seurat_data$sample <- factor(seurat_data$sample, levels = sample_order$sample)

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
pdf(file.path(image_dir, "1F_status_celltype_barplot.pdf"),
    width = 8, height = 8)
barplot <- stacked_barplots(seurat_data, meta_col = "celltype_cluster",
                            split_by = "Status",
                            color = cluster_celltype_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1))

print(barplot)

dev.off()

# AG+ group distribution by cell type
pdf(file.path(image_dir, "1G_ag_celltype_barplot.pdf"),
    width = 8, height = 8)
barplot <- stacked_barplots(seurat_data, meta_col = "celltype_cluster",
                            split_by = "Status",
                            color = cluster_celltype_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1))

print(barplot)

dev.off()

no_other <- subset(seurat_data, subset = tet_name_cutoff != "Other_Multi_Reactive")

`%notin%` <- Negate(`%in%`)

tetramer_colors_use <- tetramer_colors[names(tetramer_colors) %notin%
                                         c("Other_Multi_Reactive", 
                                           "Islet_Reactive")]

names(tetramer_colors_use) <- make.names(names(tetramer_colors_use))

barplot <- stacked_barplots(no_other, meta_col = "tet_name_cutoff",
                            split_by = "final_celltype",
                            color = tetramer_colors_use) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1))

pdf(file.path(image_dir, "1G_tetramer_stacked_barplot.pdf"),
    height = 8, width = 8)
print(barplot)
dev.off()


# OTHER
values <- seurat_data[[]] %>%
  dplyr::select(sample, Status, tet_name_cutoff) %>%
  dplyr::group_by(sample, Status) %>%
  dplyr::add_count(name = "sample_count") %>%
  dplyr::group_by(sample, Status, tet_name_cutoff) %>%
  dplyr::add_count(name = "tet_count") %>%
  dplyr::distinct() %>%
  dplyr::mutate(percent = tet_count / sample_count * 100)
  
write.csv(values, file.path(image_dir, "antigen_binding_percents.csv"))

p1 <- ggplot2::ggplot(values, ggplot2::aes(x = tet_name_cutoff, y = percent,
                                     fill = Status)) +
  ggplot2::geom_boxplot() +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.75)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1))

values <- seurat_data[[]] %>%
  dplyr::select(sample, Status, tet_name_cutoff) %>%
  dplyr::filter(tet_name_cutoff != "Negative") %>%
  dplyr::group_by(sample, Status) %>%
  dplyr::add_count(name = "sample_count") %>%
  dplyr::group_by(sample, Status, tet_name_cutoff) %>%
  dplyr::add_count(name = "tet_count") %>%
  dplyr::distinct() %>%
  dplyr::mutate(percent = tet_count / sample_count * 100)

p2 <- ggplot2::ggplot(values, ggplot2::aes(x = tet_name_cutoff, y = percent,
                                     fill = Status)) +
  ggplot2::geom_boxplot() +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.75)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1))


values <- seurat_data[[]] %>%
  dplyr::select(sample, Status, tet_name_cutoff) %>%
  dplyr::filter(tet_name_cutoff != "Negative", 
                tet_name_cutoff != "Other_Multi_Reactive") %>%
  dplyr::group_by(sample, Status) %>%
  dplyr::add_count(name = "sample_count") %>%
  dplyr::group_by(sample, Status, tet_name_cutoff) %>%
  dplyr::add_count(name = "tet_count") %>%
  dplyr::distinct() %>%
  dplyr::mutate(percent = tet_count / sample_count * 100)

p3 <- ggplot2::ggplot(values, ggplot2::aes(x = tet_name_cutoff, y = percent,
                                     fill = Status)) +
  ggplot2::geom_boxplot() +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.75)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1))

p4 <- ggplot2::ggplot(values, ggplot2::aes(x = tet_name_cutoff, y = percent,
                                           fill = Status, label = sample)) +
  ggplot2::geom_boxplot() +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.75)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1)) +
  ggrepel::geom_label_repel(size = 2, position = position_dodge(0.75))

pdf(file.path(image_dir, "tetramer_percent_by_status.pdf"))
print(p1)
print(p2)
print(p3)
print(p4)

dev.off()

facs_percents <- read.table(here("files/FACS_percent_PE_pos.csv"),
                            header = TRUE, sep = ",") %>%
  tidyr::pivot_longer(names_to = "old_status", values_to = "percent",
                      cols = colnames(.))
facs_number <- read.table(here("files/FACS_number_PE_pos.csv"),
                          header = TRUE, sep = ",") %>%
  tidyr::pivot_longer(names_to = "old_status", values_to = "count",
                      cols = colnames(.))

new_status <- c("FDR" = "ND",
                "Aab" = "AAB",
                "T1D" = "T1D")

facs_percents$Status <- new_status[facs_percents$old_status]
facs_number$Status <- new_status[facs_number$old_status]

facs_percents$type <- "PE_positive"
facs_number$type <- "PE_positive"


p1 <- ggplot2::ggplot(facs_percents, ggplot2::aes(x = type, y = percent,
                                           fill = Status)) +
  ggplot2::geom_boxplot() +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.75)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1))

p2 <- ggplot2::ggplot(facs_number, ggplot2::aes(x = type, y = count,
                                                  fill = Status)) +
  ggplot2::geom_boxplot() +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.75)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1))

pdf(file.path(image_dir, "1C_ag_positive_facs.pdf"),
    width = 6, height = 8)
print(p1)
dev.off()

pdf(file.path(image_dir, "ag_positive_count_facs.pdf"),
    width = 6, height = 8)
print(p2)
dev.off()

## Figure S1 -------------------------------------------------------------------

### QC plots -------------------------------------------------------------------
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

pdf(file.path(image_dir, "Supp1_quality_plots.pdf"))
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

pdf(file.path(image_dir, "Supp2_quality_umaps.pdf"),
    width = 12, height = 10)

print(full_plot)

dev.off()

## Figure 2 -------------------------------------------------------------------
# Read in de genes
markers_sig <- read.csv(file.path(save_dir, "files", "mast_de.csv"))

# Make heatmap of all DE
de_genes <- markers_sig[markers_sig$cluster !=
                          "T1D_Islet_Reactive_AAB_Islet_Reactive",]$gene

new_sample_order <- c("110", "116", "108", "107", "113", "114", "118",
                      "106", "117", "115", "105", "111", "102", "112",
                      "119", "109")

seurat_data$sample <- factor(seurat_data$sample, levels = new_sample_order)

seurat_islet <- subset(seurat_data, subset = tet_name_cutoff %in%
                         c("Islet_Multi_Reactive", "IA2.tet",
                           "GAD.tet", "INS.tet"))

heatmap_data <- plot_heatmap(seurat_islet, gene_list = de_genes,
                             colors = sample_colors, meta_col = "sample",
                             average_expression = TRUE, assay = "RNA",
                             plot_rownames = FALSE, cluster_rows = TRUE, 
                             return_data = TRUE)
blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", 
                "#2DB7A3", "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")

sample_info <- sample_order %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("sample")

rownames(sample_info) <- make.names(rownames(sample_info))

annotation_colors <- list("Status" = status_colors)

heatmap <- pheatmap::pheatmap(heatmap_data$z_score, cluster_rows = TRUE, 
                              cluster_cols = FALSE, show_rownames = FALSE, 
                              show_colnames = TRUE, annotation_col = sample_info, 
                              annotation_colors = annotation_colors, 
                              color = blueYellow, border_color = NA, 
                              clustering_method = "complete", silent = TRUE)

graphics.off()

pdf(file.path(image_dir, "2A_de_heatmap_average.pdf"),
    width = 8, height = 12)

print(heatmap)

dev.off()

make_heatmap <- function(kegg_list = NULL, gene_list = NULL){

  if(!is.null(kegg_list)){
    # Pull out genes from kegg list
    path <- keggLink("pathway", "hsa")
    
    all_genes <- path[grepl(kegg_list, path)]
    
    all_gene_id <- keggConv( "ncbi-geneid", names(all_genes))
    
    all_gene_id <- gsub("ncbi-geneid:", "", all_gene_id)
    
    gene_symbols <- mapIds(org.Hs.eg.db, keys = all_gene_id, 
                           keytype = "ENTREZID", column = "SYMBOL")
    
    gene_symbols <- gene_symbols[!grepl("LOC[0-9]+", gene_symbols)]    
  } else if (!is.null(gene_list)){
    gene_symbols <- gene_list
  } else {
    stop("Must provide gene list or kegg id")
  }

  
  
  heatmap_data <- plot_heatmap(seurat_islet, gene_list = gene_symbols,
                               colors = sample_colors, meta_col = "sample",
                               average_expression = TRUE, assay = "RNA",
                               plot_rownames = FALSE, cluster_rows = TRUE, 
                               return_data = TRUE)
  blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", 
                  "#2DB7A3", "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")
  
  sample_info <- sample_order %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("sample")
  
  rownames(sample_info) <- make.names(rownames(sample_info))
  
  annotation_colors <- list("Status" = status_colors)
  
  heatmap <- pheatmap::pheatmap(heatmap_data$z_score, cluster_rows = TRUE, 
                                cluster_cols = FALSE, show_rownames = TRUE, 
                                show_colnames = TRUE, annotation_col = sample_info, 
                                annotation_colors = annotation_colors, 
                                color = blueYellow, border_color = NA, 
                                clustering_method = "complete", silent = TRUE)
  return(heatmap)
  
}

# Heatmaps of gsea gene lists
paths_heatmap <- c("bcr" = "hsa04662", "antigen_processing" = "hsa04612",
                   "EBV" = "hsa05169")

all_heatmaps <- lapply(paths_heatmap, function(x){
  make_heatmap(kegg_list = x)
})

pro_anti_inflammatory <- c("IFNG", "IL1B", "CCL2", "IL23A", "CXCL13", 
                           "IL12A", "EBI3", "TGFB1", "IL10", "IL11", 
                           "IL6", "IL13", "IL1RA", "TNFA")

cytokine_heatmap <- make_heatmap(gene_list = pro_anti_inflammatory)

graphics.off()

pdf(file.path(image_dir, "2B_de_heatmap_bcr.pdf"),
    width = 8, height = 15)

print(all_heatmaps$bcr)

dev.off()
graphics.off()


pdf(file.path(image_dir, "2C_de_heatmap_antigen_processing.pdf"),
    width = 8, height = 12)

print(all_heatmaps$antigen_processing)

dev.off()
graphics.off()


pdf(file.path(image_dir, "2D_de_heatmap_cytokines.pdf"),
    width = 8, height = 3)

print(cytokine_heatmap)

dev.off()
graphics.off()


pdf(file.path(image_dir, "2E_de_heatmap_ebv.pdf"),
    width = 8, height = 30)

print(all_heatmaps$EBV)

dev.off()
graphics.off()


# 4.	UMAP - Antigen reactive cells cluster independently (not clustering by mRNA)
# INS/GAD/IA2/Multi-islet/TET/DNA
# * We decided thta this isn't a necessary plot
# 5.	Stacked bar chart: x axis = B cell subtype (updated), y axis = % of 
# total cells minus other_multi_reactive
# (INS/GAD/IA2/Multi-islet/TET/DNA/Negative)



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
                         col_by = "tet_name_cutoff",
                         sep_by = "tet_name_cutoff",
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
                         subset = tet_name_cutoff %in% c("Islet_Multi_Reactive",
                                                     "IA2-tet", "GAD-tet",
                                                     "INS-tet"))


islet_cells <- stacked_barplots(isotype, meta_col = "imcantation_isotype",
                                split_by = "Status", color = isotype_colors)

pdf(file.path(image_dir, "isotype_barplot_islet_reactive.pdf"),
    height = 8, width = 8)

print(islet_cells)

dev.off()

isotype$celltype_cluster <- factor(isotype$celltype_cluster, levels =
                                     c("Naive_0", "Naive_1", "Naive_3", 
                                       "Naive_5", "Naive_6", "ABC_7",
                                       "Resting_Memory_2", "Resting_Memory_8",
                                       "Resting_Memory_9",
                                       "Memory_4", "Plasmablast_10",
                                       "Plasmablast_11"))


barplot <- stacked_barplots(isotype, meta_col = "imcantation_isotype",
                            split_by = "celltype_cluster",
                            color = isotype_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1))

pdf(file.path(image_dir, "isotype_barplot_celltype.pdf"),
    height = 8, width = 8)

print(barplot)

dev.off()

seurat_data$celltype_cluster <- factor(seurat_data$celltype_cluster, levels =
                                     c("Naive_0", "Naive_1", "Naive_3", 
                                       "Naive_5", "Naive_6", "ABC_7",
                                       "Resting_Memory_2", "Resting_Memory_8",
                                       "Resting_Memory_9",
                                       "Memory_4", "Plasmablast_10",
                                       "Plasmablast_11"))

p1 <- featDistPlot(seurat_data, geneset = c("IgD", "IgM", "CD27.1"),
             assay = "ADT", sep_by = "celltype_cluster",
             color = cluster_celltype_colors, combine = FALSE)

pdf(file.path(image_dir, "adts_by_cluster.pdf"),
    height = 8, width = 8)

print(p1)

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
                  "clonotype_id", "Status", 
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