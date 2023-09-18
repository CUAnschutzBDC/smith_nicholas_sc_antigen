# Check CLR vs DSB normalization
library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)
library(viridis)


plot_heatmap2 <- function(seurat_object, gene_list, meta_col,
                          colors = NULL, meta_df = NULL, color_list = NULL,
                          max_val = 2.5, min_val = -2.5, cluster_rows = FALSE,
                          cluster_cols = FALSE, average_expression = FALSE,
                          plot_meta_col = TRUE, plot_rownames = TRUE,
                          cell_order = NULL, return_data = FALSE, 
                          assay = "RNA", ...){
  
  if(! assay %in% Assays(seurat_object)){
    stop(paste0("Your chosen assay: ", assay, " is not in the provided object"))
  }
  
  
  if(average_expression){
    # Find average expression of genes in clusters
    Idents(seurat_object) <- meta_col
    heatmap_df <- AverageExpression(seurat_object, seurat = FALSE,
                                    group.by = meta_col)
    heatmap_df <- heatmap_df[[assay]]
    # Test if the colnames look like integers
    character_vals <- 
      suppressWarnings(all(!is.na(as.numeric(as.character(colnames(heatmap_df))))))
    if(is.null(meta_df)){
      sample_info <- seurat_object[[meta_col]]
      # Add levels
      if(is.null(levels(sample_info[[meta_col]]))){
        sample_info[[meta_col]] <- factor(sample_info[[meta_col]])
      }
      meta_df <- data.frame(levels(sample_info[[meta_col]]))
      colnames(meta_df) <- meta_col
      if(character_vals){
        rownames(meta_df) <- paste0("X", meta_df[[meta_col]])
      } else {
        rownames(meta_df) <- meta_df[[meta_col]]
      }
      if(is.null(colors)){
        colors <- brewer.pal(length(levels(sample_info[[meta_col]])), "Set1")
        names(colors) <- levels(sample_info[[meta_col]])
      } 
      # make a list for the column labeing
      color_list <- list(colors)
      names(color_list) <- meta_col
    }
  } else {
    # Pull out data and subset to genes of interest
    heatmap_df <- GetAssayData(seurat_object, slot = "data", assay = assay)
  }
  heatmap_df <- heatmap_df[rownames(heatmap_df) %in% gene_list, ]
  heatmap_df <- data.frame(heatmap_df)
  
  heatmap_df <- heatmap_df[order(match(rownames(heatmap_df), gene_list)), ]
  
  # remove any zero values
  heatmap_df <- heatmap_df[rowSums(heatmap_df) > 0,]
  
  if(is.null(meta_df)){
    # Make a df for the column labeling
    sample_info <- seurat_object[[meta_col]]
    # Add levels
    if(is.null(levels(sample_info[[meta_col]]))){
      sample_info[[meta_col]] <- factor(sample_info[[meta_col]])
    }
    if(is.null(colors)){
      colorcount <- length(levels(sample_info[[meta_col]]))
      colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(colorcount)
      names(colors) <- levels(sample_info[[meta_col]])
    } 
    # make a list for the column labeing
    coloring <- list(colors)
    names(coloring) <- meta_col
  } else {
    # The sample info and color list must be provided
    sample_info <- meta_df
    coloring <- color_list
  }
  
  # Set cluster order
  
  cluster_order <- levels(sample_info[[meta_col]])
  heatmap_scale <- t(scale(t(heatmap_df), scale = TRUE))
  # Colors for heatmap (from the ArchR package)
  blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
                  "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")
  
  if(!is.null(cell_order)){
    sample_info <- sample_info[order(match(rownames(sample_info),
                                           cell_order)), , drop = FALSE]
    rownames(sample_info) <- make.names(rownames(sample_info))
    if (!identical(colnames(heatmap_scale), rownames(sample_info))) {
      heatmap_scale <- heatmap_scale[, rownames(sample_info)]
    }
  } else if (!cluster_cols) {
    sample_info <- sample_info[order(match(sample_info[[meta_col]], 
                                           cluster_order)), , drop = FALSE]
    rownames(sample_info) <- make.names(rownames(sample_info))
    if (!identical(colnames(heatmap_scale), rownames(sample_info))) {
      heatmap_scale <- heatmap_scale[, rownames(sample_info)]
    }
  } else {
    rownames(sample_info) <- make.names(rownames(sample_info))
  }
  
  if(!plot_meta_col){
    sample_info[[meta_col]] <- NULL
  }
  
  # This makes the values more even
  heatmap_scale <- ifelse(heatmap_scale > max_val, max_val, heatmap_scale)
  heatmap_scale <- ifelse(heatmap_scale < min_val, min_val, heatmap_scale)
  
  heatmap <- pheatmap::pheatmap(heatmap_scale, cluster_rows = cluster_rows,
                                cluster_cols = cluster_cols,
                                show_rownames = plot_rownames,
                                show_colnames = FALSE, annotation_col = sample_info,
                                annotation_colors = coloring, color = blueYellow,
                                border_color = NA, clustering_method = "complete",
                                silent = TRUE, ...)
  if(return_data){
    return(list("heatmap" = heatmap,
                "z_score" = heatmap_scale,
                "counts" = heatmap_df))
  } else {
    return(heatmap)
  }
}

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log
# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "JH310-12_AP"

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- args[[4]]
#sample_info <- here("files/sample_info.tsv")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident

ADT_pca <- sample_info$adt_pca
run_adt_umap <- sample_info$adt_umap

RNA_pcs <- sample_info$PCs
resolution <- sample_info$resolution


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


dsb_norm <- GetAssayData(seurat_data, slot = "data", assay = "DSB_ADT")
dsb_norm <- t(dsb_norm) %>% data.frame()
# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

raw_vals <- GetAssayData(seurat_data, slot = "counts", assay = "CLR_ADT")

raw_vals <- log1p(raw_vals) %>% as.matrix

raw_vals <- data.frame(t(raw_vals))

clr_vals <- GetAssayData(seurat_data, slot = "data", assay = "CLR_ADT")

clr_vals <- data.frame(t(clr_vals))

scar_vals <- GetAssayData(seurat_data, slot = "data", assay = "SCAR_ADT_LOG")

scar_vals <- data.frame(t(as.matrix(scar_vals)))

scar_clr <- GetAssayData(seurat_data, slot = "data", assay = "SCAR_ADT")

scar_clr <- data.frame(t(scar_clr))


# Plots ---------------------------------------------
p1 <- ggplot2::ggplot(dsb_norm, ggplot2::aes(x = IgM, y = IgD)) +
  ggplot2::geom_point() +
  ggplot2::ggtitle("dsb norm")

p2 <- ggplot2::ggplot(raw_vals, ggplot2::aes(x = IgM, y = IgD)) +
  ggplot2::geom_point() +
  ggplot2::ggtitle("raw vals")

p3 <- ggplot2::ggplot(clr_vals, ggplot2::aes(x = IgM, y = IgD)) +
  ggplot2::geom_point() +
  ggplot2::ggtitle("clr norm")

p4 <- ggplot2::ggplot(scar_vals, ggplot2::aes(x = IgM, y = IgD)) +
  ggplot2::geom_point() +
  ggplot2::ggtitle("scar log norm")

p5 <- ggplot2::ggplot(scar_clr, ggplot2::aes(x = IgM, y = IgD)) +
  ggplot2::geom_point() +
  ggplot2::ggtitle("scar clr norm")


p6 <- ggplot2::ggplot(dsb_norm, ggplot2::aes(x = IA2.tet, y = GAD.tet)) +
  ggplot2::geom_point() +
  ggplot2::ggtitle("dsb norm")

p7 <- ggplot2::ggplot(raw_vals, ggplot2::aes(x = IA2.tet, y = GAD.tet)) +
  ggplot2::geom_point() +
  ggplot2::ggtitle("raw vals")

p8 <- ggplot2::ggplot(clr_vals, ggplot2::aes(x = IA2.tet, y = GAD.tet)) +
  ggplot2::geom_point() +
  ggplot2::ggtitle("clr norm")

p9 <- ggplot2::ggplot(scar_vals, ggplot2::aes(x = IA2.tet, y = GAD.tet)) +
  ggplot2::geom_point() +
  ggplot2::ggtitle("scar log norm")

p10 <- ggplot2::ggplot(scar_clr, ggplot2::aes(x = IA2.tet, y = GAD.tet)) +
  ggplot2::geom_point() +
  ggplot2::ggtitle("scar clr norm")

# Heatmaps ---------------------------------------------------------------------
pdf(file.path(save_dir, "images", "clr_vs_dsb_normalization.pdf"),
    height = 8, width = 8.5)

DefaultAssay(seurat_data) <- "TET"
genes_plot <- rownames(seurat_data)
DefaultAssay(seurat_data) <- "RNA"

sample_info <- seurat_data[[]] %>%
  dplyr::select(tet_hash_id, scar_hash_id, dplyr::matches("tet.*_clusters")) %>%
  dplyr::select(!tetdsb_nd_clusters)

set1_palette <- grDevices::colorRampPalette(
  colors = RColorBrewer::brewer.pal(n = 9, name = "Set1")
  )

tet_clr <- set1_palette(length(unique(sample_info$tetclr_clusters)))
names(tet_clr) <- unique(sample_info$tetclr_clusters)

tet_dsb <- set1_palette(length(unique(sample_info$tetdsb_clusters)))
names(tet_dsb) <- unique(sample_info$tetdsb_clusters)

tet_scar <- set1_palette(length(unique(sample_info$tetscar_clusters)))
names(tet_scar) <- unique(sample_info$tetscar_clusters)

tet_scar_log <- set1_palette(length(unique(sample_info$tetscarlog_clusters)))
names(tet_scar_log) <- unique(sample_info$tetscarlog_clusters)

hash_tet <- set1_palette(length(unique(sample_info$tet_hash_id)))
names(hash_tet) <- unique(sample_info$tet_hash_id)

hash_scar <- set1_palette(length(unique(sample_info$scar_hash_id)))
names(hash_scar) <- unique(sample_info$scar_hash_id)

sample_info <- rev(sample_info)

all_colors <- list("tet_hash_id" = hash_tet,
                   "scar_hash_id" = hash_tet,
                   "tetdsb_clusters" = tet_dsb,
                   "tetclr_clusters" = tet_clr,
                   "tetscar_clusters" = tet_scar,
                   "tetscarlog_clusters" = tet_scar_log)

plot_heatmap2(seurat_object = seurat_data, gene_list = genes_plot, 
              meta_col = "tet_hash_id", colors = NULL, 
              meta_df = sample_info, color_list = all_colors, max_val = 2.5, 
              min_val = -2.5,  cluster_rows = TRUE, 
              cluster_cols = TRUE, average_expression = FALSE, 
              plot_meta_col = TRUE, plot_rownames = TRUE,
              cell_order = NULL, return_data = FALSE, 
              assay = "CLR_TET", main = "CLR")

grid::grid.newpage()

plot_heatmap2(seurat_object = seurat_data, gene_list = genes_plot, 
              meta_col = "tet_hash_id", colors = NULL, 
              meta_df = sample_info, color_list = all_colors, max_val = 2.5, 
              min_val = -2.5,  cluster_rows = TRUE, 
              cluster_cols = TRUE, average_expression = FALSE, 
              plot_meta_col = TRUE, plot_rownames = TRUE,
              cell_order = NULL, return_data = FALSE, 
              assay = "DSB_TET", main = "DSB")

grid::grid.newpage()

plot_heatmap2(seurat_object = seurat_data, gene_list = genes_plot, 
              meta_col = "tet_hash_id", colors = NULL, 
              meta_df = sample_info, color_list = all_colors, max_val = 2.5, 
              min_val = -2.5,  cluster_rows = TRUE, 
              cluster_cols = TRUE, average_expression = FALSE, 
              plot_meta_col = TRUE, plot_rownames = TRUE,
              cell_order = NULL, return_data = FALSE, 
              assay = "SCAR_TET", main = "SCAR CLR")

grid::grid.newpage()

plot_heatmap2(seurat_object = seurat_data, gene_list = genes_plot, 
              meta_col = "tet_hash_id", colors = NULL, 
              meta_df = sample_info, color_list = all_colors, max_val = 2.5, 
              min_val = -2.5,  cluster_rows = TRUE, 
              cluster_cols = TRUE, average_expression = FALSE, 
              plot_meta_col = TRUE, plot_rownames = TRUE,
              cell_order = NULL, return_data = FALSE, 
              assay = "SCAR_TET_LOG", main = "SCAR LOG")

## Heatmap with nd -------------------------------------------------------------
seurat_no_doub <- subset(seurat_data, subset = tet_hash_id != "Doublet")

sample_info <- seurat_no_doub[[]] %>%
  dplyr::select(tet_hash_id, dplyr::matches("tetdsb_nd_clusters")) 

set1_palette <- grDevices::colorRampPalette(
  colors = RColorBrewer::brewer.pal(n = 9, name = "Set1")
)

tet_dsb <- set1_palette(length(unique(sample_info$tetdsb_nd_clusters)))
names(tet_dsb) <- unique(sample_info$tetdsb_nd_clusters)

hash <- set1_palette(length(unique(sample_info$tet_hash_id)))
names(hash) <- unique(sample_info$tet_hash_id)

all_colors <- list("tetdsb_nd_clusters" = tet_dsb,
                   "tet_hash_id" = hash)

grid::grid.newpage()

plot_heatmap2(seurat_object = seurat_no_doub, gene_list = genes_plot, 
              meta_col = "tet_hash_id", colors = NULL, 
              meta_df = sample_info, color_list = all_colors, max_val = 2.5, 
              min_val = -2.5,  cluster_rows = TRUE, 
              cluster_cols = TRUE, average_expression = FALSE, 
              plot_meta_col = TRUE, plot_rownames = TRUE,
              cell_order = NULL, return_data = FALSE, 
              assay = "DSB_TET", main = "DSB")

# Correlations -----------------------------------------------------------------
corr_dsb <- cor(x = dsb_norm)

pheatmap::pheatmap(corr_dsb, main = "dsb")

#grid::grid.newpage()
corr_raw <- cor(x = raw_vals)

pheatmap::pheatmap(corr_raw, main = "raw")

#grid::grid.newpage()

corr_clr <- cor(x = clr_vals)

pheatmap::pheatmap(corr_clr, main = "clr")

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
print(p9)
print(p10)

dev.off()