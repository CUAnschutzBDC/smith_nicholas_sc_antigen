library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(clustree)
library(harmony)
library(singlecellmethods)
library(batchelor)
library(clustifyr)


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

remove_ambience <- TRUE
batch_columns <- c("sample", "date.processed.for.scSeq")
cor_cutoff <- 0.5

all_ref_dir <-
  "/beevol/home/wellskri/Analysis/references/single_cell_references"
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
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))

celltype_levels <- c("Naive_1", "Naive_3",
                     "BND2", "B.intermediate",
                     "Memory_IgA", "Resting_memory", 
                     "Activated_memory",
                     "Plasmablast",
                     "CD14.Mono", "CD8.TEM")

seurat_data <- subset(seurat_data, 
                      subset = RNA_combined_celltype %in% celltype_levels[1:8])

pca_list <- list("rna" = c(assay = seurat_assay, reduction_name = "pca"))

if(remove_ambience){
  pca_list <- c(pca_list,
                list("ambience" = c(assay = "AMBRNA", 
                                    reduction_name = "ambpca")))
} 
# Colors -----------------------------------------------------------------------

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

all_celltypes <- unique(seurat_data$RNA_celltype)
other_celltypes <- all_celltypes[!all_celltypes %in% names(new_colors)]

set_one_brewer <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9,
                                                                       name = "Set1"))
other_colors <- set_one_brewer(length(other_celltypes))

names(other_colors) <- other_celltypes

celltype_colors <- readRDS(file = file.path(save_dir, "color_palette.rds"))


# PCA --------------------------------------------------------------------------
for(assay_type in names(pca_list)){
  seurat_assay <- pca_list[[assay_type]][["assay"]]
  reduction_name <- pca_list[[assay_type]][["reduction_name"]]
  DefaultAssay(seurat_data) <- seurat_assay
  
  seurat_data <- FindVariableFeatures(seurat_data,
                                      assay = seurat_assay)
  
  seurat_data <- ScaleData(seurat_data, assay = seurat_assay)
  
  VariableFeatures(seurat_data) <- VariableFeatures(seurat_data)[!grepl("IG[H|L|K]",
                                                                        VariableFeatures(seurat_data))]
  
  # PCA of gene expression, weighted by cell number per sample
  seurat_data <- singlecellmethods::RunBalancedPCA(obj = seurat_data, 
                                                   weight.by = "orig.ident",
                                                   npcs = 50,
                                                   assay.use = seurat_assay,
                                                   reduction.name = reduction_name)
  
  seurat_data$date.processed.for.scSeq <- factor(seurat_data$date.processed.for.scSeq)
}


if(ADT & run_adt_umap){
  # PCA of surface protein
  seurat_data <- PCA_dimRed(seurat_data, assay = "ADT")
  
  ADT_plots <- plot_PCA(HTO = HTO, assay = "ADT", sample_object = seurat_data)
  
  pdf(file.path(save_dir, "images/ADT_pca.pdf"))
  plot(ADT_plots)
  dev.off()
}

# Batch Correction -------------------------------------------------------------
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
  pdf(file.path(save_dir, "images", paste0(assay_type, "_harmony_convergence_repeat.pdf")))
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
                paste0(assay_type, "_testing_batch_correction_repeat.pdf")))
  
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
                       "_testing_batch_correction_separated_plot_repeat.pdf")),
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

# UMAP -------------------------------------------------------------------------
for(assay_type in names(pca_list)){
  seurat_assay <- pca_list[[assay_type]][["assay"]]
  reduction_name <- pca_list[[assay_type]][["reduction_name"]]
  
  reduction_save <- paste(assay_type, batch_correction, sep = "_")
  
  # Remove previous clustering
  remove_cols <- colnames(seurat_data[[]])[grepl("res\\.[0-9]",
                                                 colnames(seurat_data[[]]))]
  
  for (i in remove_cols){
    seurat_data[[i]] <- NULL
  }
  
  ## Cluster --------------------------------------------------------------------
  
  # UMAP of gene expression
  set.seed(0)
  umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                           resolution = resolution, assay = seurat_assay,
                           HTO = HTO,
                           reduction = reduction_save)
  seurat_data <- umap_data$object
  
  plot_type <- paste0(reduction_save, ".umap")
  
  seurat_data[[paste0(assay_type, "_corrected_cluster")]] <- 
    seurat_data[[paste0(seurat_assay, "_cluster")]][[1]]
}

if(ADT & run_adt_umap){
  # UMAP of surface protein
  umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = ADT_pcs,
                           resolution = 0.6, assay = "ADT", HTO = TRUE)
  
  seurat_data <- umap_data$object
  
  adt_plots <- umap_data$plots
  
  
  # UMAP of combined
  # Identify multimodal neighbors. These will be stored in the neighbors slot, 
  # and can be accessed using bm[['weighted.nn']]
  # The WNN graph can be accessed at bm[["wknn"]], 
  # and the SNN graph used for clustering at bm[["wsnn"]]
  # Cell-specific modality weights can be accessed at bm$RNA.weight
  if(SCT){
    pca_slot <- "sctpca"
    weight_name <- "SCT.weight"
  } else{
    pca_slot <- "pca"
    weight_name <- "RNA.weight"
  }
  seurat_data <- FindMultiModalNeighbors(
    seurat_data, reduction.list = list(pca_slot, "apca"), 
    dims.list = list(1:RNA_pcs, 1:ADT_pcs),
    modality.weight.name = c(weight_name, "ADT.weight")
  )
  
  
  seurat_data <- RunUMAP(seurat_data, nn.name = "weighted.nn",
                         reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  seurat_data <- FindClusters(seurat_data, graph.name = "wsnn",
                              algorithm = 3, resolution = 2, verbose = FALSE)
  
  seurat_data[["combined_cluster"]] <- Idents(seurat_data)
  col_by_list <- c("combined_cluster", "orig.ident")
  if(HTO){
    col_by_list <- c(col_by_list, "HTO_classification")
  }
  save_plot <- file.path(save_dir, "images/combined_umap.pdf")
  plot_list <- plotDimRed(sample_object = seurat_data,
                          save_plot = NULL,
                          col_by = col_by_list, return_plot = TRUE,
                          plot_type = "wnn.umap")
}

# Name clusters ----------------------------------------------------------------
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
  seurat_data$RNA_combined_celltype_old <- seurat_data$RNA_combined_celltype
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


saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_reprocessed.rds"))
