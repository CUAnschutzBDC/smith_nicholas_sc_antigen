library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(clustree)
library(harmony)
library(singlecellmethods)
library(batchelor)


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

remove_ambience <- TRUE

normalization_method <- "log" # can be SCT or log

batch_columns <- c("sample", "date.processed.for.scSeq")

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


pca_list <- list("rna" = c(assay = seurat_assay, reduction_name = "pca"))

if(remove_ambience){
  pca_list <- c(pca_list,
                list("ambience" = c(assay = "AMBRNA", 
                                    reduction_name = "ambpca")))
  
} 

# Set directories
save_dir <- file.path(results_dir, "R_analysis", sample)

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

# Colors -----------------------------------------------------------------------
all_colors <- readRDS(file = file.path("files/all_colors.rds"))


celltype_colors <- all_colors$cell_type_colors
# UMAP -------------------------------------------------------------------------
for(assay_type in names(pca_list)){
  seurat_assay <- pca_list[[assay_type]][["assay"]]
  reduction_name <- pca_list[[assay_type]][["reduction_name"]]
  #seurat_data$RNA_cluster <- NULL
  
  # Remove previous clustering
  remove_cols <- colnames(seurat_data[[]])[grepl("res\\.[0-9]",
                                                 colnames(seurat_data[[]]))]
  
  for (i in remove_cols){
    seurat_data[[i]] <- NULL
  }
  
  set.seed(0)
  
  # UMAP of gene expression
  set.seed(0)
  umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                           resolution = 0.8, assay = seurat_assay, HTO = HTO,
                           reduction = reduction_name)
  
  seurat_data <- umap_data$object
  
  gene_plots <- umap_data$plots
  
  rm(umap_data)
  rm(gene_plots)
  
  seurat_data[[paste0(assay_type, "_uncorrected_cluster")]] <- 
    seurat_data[[paste0(seurat_assay, "_cluster")]][[1]]
  
  # Batch correction -------------------------------------------------------------
  
  # NOTE: Decide if batch correction is needed
  # NOTE: Test Harmony and MNN to decide what is needed.
  
  ## Harmony ---------------------------------------------------------------------
  pdf(file.path(save_dir, "images", paste0(assay_type, "_harmony_convergence.pdf")))
  seurat_data <- RunHarmony(seurat_data, c("sample"),
                            plot_convergence = TRUE,
                            theta = 5, reduction = reduction_name,
                            reduction.save = paste0(assay_type, "_harmony"),
                            assay.use = seurat_assay)
  
  harmony_plot <- plotDimRed(seurat_data, col_by = c(batch_columns,
                                                     "RNA_celltype"),
                             plot_type = paste0(assay_type, "_harmony"),
                             ggrastr = TRUE)
  
  print(harmony_plot)
  dev.off()
  
  ## MNN -------------------------------------------------------------------------
  DefaultAssay(seurat_data) <- seurat_assay
  sce_data <- as.SingleCellExperiment(seurat_data)
  
  set.seed(0)
  # Run fastMNN on batch
  corrected_data <- fastMNN(sce_data, batch = sce_data$sample,
                            subset.row = VariableFeatures(seurat_data))
  
  # Check reduced dims
  if(!identical(rownames(SingleCellExperiment::reducedDim(x = corrected_data)), 
                colnames(seurat_data))){
    stop("names are not the same between the object and fastmnn output!!!")
  }
  
  # Change names
  colnames(SingleCellExperiment::reducedDim(x = corrected_data)) <-
    paste0(assay_type, "_mnn_", 
           1:ncol(SingleCellExperiment::reducedDim(x = corrected_data)))
  
  # Add to seurat object
  seurat_data[[paste0(assay_type, "_mnn")]] <- CreateDimReducObject(
    embeddings = SingleCellExperiment::reducedDim(x = corrected_data),
    loadings = as.matrix(SingleCellExperiment::rowData(x = corrected_data)),
    assay = DefaultAssay(object = seurat_data),
    key = paste0(assay_type, "_mnn_")
  )
  
  rm(corrected_data)
  
  # UMAP -------------------------------------------------------------------------
  
  
  # Remove previous clustering
  remove_cols <- colnames(seurat_data[[]])[grepl("res\\.[0-9]",
                                                 colnames(seurat_data[[]]))]
  
  for (i in remove_cols){
    seurat_data[[i]] <- NULL
  }
  
  ## Harmony ---------------------------------------------------------------------
  set.seed(0)
  # UMAP of gene expression
  umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                           resolution = 0.6, assay = seurat_assay, HTO = HTO,
                           reduction = paste0(assay_type, "_harmony"))
  
  seurat_data <- umap_data$object
  
  gene_plots <- umap_data$plots
  
  rm(umap_data)
  rm(gene_plots)
  
  seurat_data[[paste0(assay_type, "_harmony_clust")]] <- 
    seurat_data[[paste0(seurat_assay, "_cluster")]][[1]]
  
  cm <- confusionMatrix(seurat_data[[paste0(assay_type, "_harmony_clust")]][[1]],
                        seurat_data$RNA_celltype)
  
  cm <- cm / rowSums(cm)
  
  harmony_heatmap <- pheatmap::pheatmap(cm, silent = TRUE)
  
  ## MNN -------------------------------------------------------------------------
  
  set.seed(0)
  # UMAP of gene expression
  umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                           resolution = 0.6, assay = seurat_assay, HTO = HTO,
                           reduction = paste0(assay_type, "_mnn"), size = 1)
  
  seurat_data <- umap_data$object
  
  gene_plots <- umap_data$plots
  
  rm(umap_data)
  rm(gene_plots)
  
  seurat_data[[paste0(assay_type, "_mnn_clust")]] <-
    seurat_data[[paste0(seurat_assay, "_cluster")]][[1]]
  
  cm <- confusionMatrix(seurat_data[[paste0(assay_type, "_mnn_clust")]][[1]],
                        seurat_data$RNA_celltype)  
  cm <- cm / rowSums(cm)
  
  mnn_heatmap <- pheatmap::pheatmap(cm, silent = TRUE)
  
  ## Compare --------------------------------------------------------------------
  
  cm <- confusionMatrix(seurat_data[[paste0(assay_type, "_harmony_clust")]][[1]],
                        seurat_data[[paste0(assay_type, "_mnn_clust")]][[1]])  
  cm <- cm / rowSums(cm)
  
  compare_heatmap <- pheatmap::pheatmap(cm, silent = TRUE)
  
  pdf(file.path(save_dir, "images",
                paste0(assay_type, "_testing_batch_correction.pdf")))
  
  print(plotDimRed(seurat_data, col_by = "RNA_celltype",
                   plot_type = paste0(reduction_name, ".umap"),
                   ggrastr = TRUE, color = celltype_colors))
  print(plotDimRed(seurat_data, col_by = batch_columns, 
                   plot_type = paste0(reduction_name, ".umap"),
                   ggrastr = TRUE))
  
  
  grid::grid.newpage()
  print(harmony_heatmap)
  print(plotDimRed(seurat_data, col_by = "RNA_celltype",
                   plot_type = paste0(assay_type, "_harmony.umap"),
                   ggrastr = TRUE))
  print(plotDimRed(seurat_data, col_by = batch_columns, 
                   plot_type = paste0(assay_type, "_harmony.umap"),
                   ggrastr = TRUE))
  
  grid::grid.newpage()
  print(mnn_heatmap)
  
  print(plotDimRed(seurat_data, col_by = "RNA_celltype", 
                   plot_type = paste0(assay_type, "_mnn.umap"),
                   ggrastr = TRUE, color = celltype_colors))
  print(plotDimRed(seurat_data, col_by = batch_columns,
                   plot_type = paste0(assay_type, "_mnn.umap"),
                   ggrastr = TRUE))
  
  grid::grid.newpage()
  print(compare_heatmap)
  
  dev.off()
  
  pdf(file.path(save_dir, "images", 
                paste0(assay_type, 
                       "_testing_batch_correction_separated_plot.pdf")),
      width = 15, height = 15)
  all_plots <- lapply(unique(seurat_data$sample), function(x){
    plotDimRed(seurat_data, col_by = "RNA_celltype",
               plot_type = paste0(reduction_name, ".umap"),
               highlight_group = TRUE, meta_data_col = "sample",
               group = x, ggrastr = TRUE, color = celltype_colors)[[1]]
  })
  
  print(cowplot::plot_grid(plotlist = all_plots,
                           nrow = length(unique(seurat_data$sample))/4))
  
  
  all_plots <- lapply(unique(seurat_data$sample), function(x){
    plotDimRed(seurat_data, col_by = "RNA_celltype",
               plot_type = paste0(assay_type, "_harmony.umap"),
               highlight_group = TRUE, meta_data_col = "sample",
               group = x, ggrastr = TRUE, color = celltype_colors)[[1]]
  })
  
  print(cowplot::plot_grid(plotlist = all_plots, 
                           nrow = length(unique(seurat_data$sample))/4))
  
  
  all_plots <- lapply(unique(seurat_data$sample), function(x){
    plotDimRed(seurat_data, col_by = "RNA_celltype", 
               plot_type = paste0(assay_type, "_mnn.umap"),
               highlight_group = TRUE, meta_data_col = "sample",
               group = x, ggrastr = TRUE, color = celltype_colors)[[1]]
  })
  
  print(cowplot::plot_grid(plotlist = all_plots, 
                           nrow = length(unique(seurat_data$sample))/4))
  
  dev.off()
  
}
saveRDS(seurat_data, file.path(save_dir, "rda_obj/seurat_processed.rds"))
