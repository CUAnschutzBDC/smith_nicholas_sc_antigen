library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)
library(viridis)
library(clustifyr)
library(singlecellmethods)
library(batchelor)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

remove_ambience <- TRUE

cor_cutoff <- 0.5

all_ref_dir <-
  "/beevol/home/wellskri/Analysis/references/single_cell_references"

normalization_method <- "log" # can be SCT or log
# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

#args <- c("merged", here("results"), "", here("files/sample_info.tsv"))


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
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

pca_list <- list("rna" = c(assay = seurat_assay, reduction_name = "pca"))

if(remove_ambience){
  pca_list <- c(pca_list,
                list("ambience" = c(assay = "AMBRNA", 
                                    reduction_name = "ambpca")))
}


# Remove B cells with more than n = 2 chains -----------------------------------
# The chains we want to remove are:
# IGH;IGH;IGK
# IGH;IGH;IGK;IGK
# IGH;IGH;IGK;IGL
# IGH;IGH;IGL
# IGH;IGH;IGL;IGL
# IGH;IGK;IGK
# IGH;IGL;IGL
cm <- confusionMatrix(seurat_data$chains, seurat_data$Doublet_finder)
cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

cm <- confusionMatrix(seurat_data$Doublet_finder, seurat_data$chains)
cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

# So remove all the chains with multiple heavy, remove "doublets" from
# all non-plasmablasts

remove_chains <- c("IGH;IGH;IGK", "IGH;IGH;IGK;IGK", "IGH;IGH;IGK;IGL",
                   "IGH;IGH;IGL", "IGH;IGH;IGL;IGL")

remove_cells_chain <- seurat_data[[]] %>%
  dplyr::filter(chains %in% remove_chains)

remove_cells_chain <- rownames(remove_cells_chain)

remove_cells_doublet <- seurat_data[[]] %>%
  dplyr::filter(Doublet_finder == "Doublet" & 
                  RNA_combined_celltype != "Plasmablast")

remove_cells_doublet <- rownames(remove_cells_doublet)

remove_cells <- unique(c(remove_cells_chain, remove_cells_doublet))

keep_cells <- colnames(seurat_data)[!colnames(seurat_data) %in% remove_cells]

seurat_data <- subset(seurat_data, cells = keep_cells)

# Remove Other Multi Reactive
seurat_data$tet_name_cutoff <- ifelse(seurat_data$tet_name_cutoff == "Other_Multi_Reactive" & 
                                        grepl("INS|GAD|IA2", seurat_data$full_tet_name_cutoff),
                                      "Islet_Multi_Reactive", seurat_data$tet_name_cutoff)

seurat_data <- subset(seurat_data, 
                      subset = tet_name_cutoff != "Other_Multi_Reactive")


#-------------------------------------------------------------------------------
for(assay_type in names(pca_list)){
  seurat_assay <- pca_list[[assay_type]][["assay"]]
  reduction_name <- pca_list[[assay_type]][["reduction_name"]]

  plot_type <- paste0(assay_type, "_", batch_correction, ".umap")

  # Reprocess data ---------------------------------------------------------------
  DefaultAssay(seurat_data) <- seurat_assay
  seurat_data <- seurat_data %>%
    FindVariableFeatures() %>%
    ScaleData(vars.to.regress = vars.to.regress)

  # Update variable features
  VariableFeatures(seurat_data) <- VariableFeatures(seurat_data)[!grepl("IG[H|L|K]",
                                                                        VariableFeatures(seurat_data))]

  ## PCA -------------------------------------------------------------------------

  # PCA of gene expression, weighted by cell number per sample
  seurat_data <- singlecellmethods::RunBalancedPCA(obj = seurat_data, 
                                                   weight.by = "Sample.Name",
                                                   npcs = 50,
                                                   assay.use = seurat_assay,
                                                   reduction.name = reduction_name)
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

  ## UMAP ------------------------------------------------------------------------

  # Remove previous clustering
  remove_cols <- colnames(seurat_data[[]])[grepl("res\\.[0-9]",
                                                 colnames(seurat_data[[]]))]

  for (i in remove_cols){
    seurat_data[[i]] <- NULL
  }

  # UMAP of gene expression
  set.seed(0)
  # UMAP of gene expression
  umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                           resolution = 0.6, assay = seurat_assay, HTO = HTO,
                           reduction = paste0(assay_type, "_mnn"), size = 1)
  
  seurat_data <- umap_data$object

  seurat_data[[paste0(assay_type, "_corrected_cluster")]] <-
    seurat_data[[paste0(seurat_assay, "_cluster")]][[1]]

  ## Name cell types -------------------------------------------------------------

  ###############
  # Seurat pbmc #
  ###############

  # Information for cell mapping
  ref_dir <- file.path(all_ref_dir, "pbmc/seurat")

  ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference_l2.csv"),
                      header = TRUE, row.names = 1)

  pdf(file.path(save_dir, "images", 
                paste0(assay_type, "no_doublet_celltype_mapping.pdf")))

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
                      paste0(assay_type, "_no_doublet_celltype_mapping_seurat.csv")))



  ####################
  # Published object #
  ####################

  ref_dir <- here("files/210825_object")


  ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference.csv"),
                      header = TRUE, row.names = 1)

  DefaultAssay(seurat_data) <- "RNA"

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
                      paste0(assay_type, "_no_doublet_celltype_mapping_bnd.csv")))


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

cell_type_mapping <- c("BND2" = "ABC",
                       "Memory_IgA" = "Memory",
                       "Memory_IgE_IgG" = "Memory",
                       "Naive_1" = "Naive",
                       "Naive_3" = "Naive",
                       "Plasmablast" = "Plasmablast",
                       "Resting_memory" = "Resting_Memory")


seurat_data$final_celltype <- cell_type_mapping[seurat_data$RNA_combined_celltype]

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))