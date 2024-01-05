# ADT and tet analysis
# 1. separate out tet
# 2. UMAP on ADT_CLR, TET_CLR, ADT_DSB, TET_DSB
# 3. Name the clusters from each with the sample name
# 4. HTO demux on TET_DSB, TET_CLR
library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)
library(viridis)

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

image_dir <- file.path(save_dir, "images", "tetramer_demux")

ifelse(!dir.exists(image_dir), dir.create(image_dir, recursive = TRUE), FALSE)

# Function to run UMAP on the specified assay
# Have the function also return plots of all members of the assay and the
# clustering, these should be violin + umap plots returned as a cowplot
adt_dim_red <- function(seurat_object, assay, reduction_name,
                        sample_name, resolution = 0.8, hash_id_col = "tet_hash_id"){
  DefaultAssay(seurat_object) <- assay
  # hack seurat to use normalized protein values as a dimensionality reduction object.
  VariableFeatures(seurat_object) <- rownames(seurat_object)
  
  nPCs <- length(VariableFeatures(seurat_object))
  
  # run true pca to initialize dr pca slot for WNN 
  seurat_object <- ScaleData(seurat_object, verbose = FALSE) %>%
    RunPCA(reduction.name = reduction_name,
           features = VariableFeatures(seurat_object),
           verbose = FALSE)
  
  # make matrix of norm values to add as dr embeddings
  pseudo <- t(as.matrix(GetAssayData(seurat_object, slot = "data")))
  pseudo_colnames <- paste('pseudo', 1:ncol(pseudo), sep = "_")
  colnames(pseudo) <- pseudo_colnames
  # add to object 
  seurat_object@reductions[[reduction_name]]@cell.embeddings = pseudo
  
  # Clustering and UMAP
  
  seurat_object <- FindNeighbors(seurat_object, 
                                 features = rownames(seurat_object), 
                                 dims = 1:nPCs, reduction = reduction_name)
  
  seurat_object <- FindClusters(seurat_object, resolution = resolution, 
                                graph.name = paste0(assay, "_snn"))
  
  plot_type <- paste0(reduction_name, ".umap")
  
  seurat_object <- RunUMAP(seurat_object, metric = "correlation", 
                           dims = 1:nPCs, reduction = reduction_name,
                           assay = assay, 
                           reduction.key = paste0(reduction_name, "UMAP_"), 
                           reduction.name = plot_type)
  
  cluster_col <- paste0(reduction_name, "_clusters")
  
  seurat_object[[cluster_col]] <- paste0(sample_name, "_",
                                         Idents(seurat_object))
  
  full_plot <- plotDimRed(seurat_object, col_by = cluster_col,
                          plot_type = plot_type)
  
  all_plots <- lapply(VariableFeatures(seurat_object), function(x){
    umap <- plotDimRed(seurat_object, col_by = x, plot_type = plot_type)
    violin <- featDistPlot(seurat_object = seurat_object, 
                           geneset = x, combine = FALSE, sep_by = cluster_col)
    
    final_plot <- cowplot::plot_grid(full_plot[[1]], NULL, umap[[1]],
                                     violin[[1]], 
                                     align = "hv", axis = "lr", 
                                     nrow = 2, ncol = 2)
    
    return(final_plot)
  })
  
  remove_clusters <- colnames(seurat_object[[]])[grepl(paste0(assay, "_snn_res"), 
                                                       colnames(seurat_object[[]]))]
  
  for(i in remove_clusters){
    seurat_object[[i]] <- NULL
  }
  
  # Test against the hto demux
  cM <- confusionMatrix(seurat_object[[cluster_col]][[1]],
                        seurat_object[[hash_id_col]][[1]])
  
  cM <- cM /rowSums(cM)
  
  heatmap1 <- pheatmap::pheatmap(cM, silent = TRUE)
  
  cM <- confusionMatrix(seurat_object[[hash_id_col]][[1]],
                        seurat_object[[cluster_col]][[1]])
  
  cM <- cM /rowSums(cM)
  
  heatmap2 <- pheatmap::pheatmap(cM, silent = TRUE)
  
  
  return(list("plots" = all_plots, "object" = seurat_object,
              "heatmap1" = heatmap1, "heatmap2" = heatmap2))
  
}

# Starting assays
# ADT (DSB)
# CLR ADT
# CLR TET (named TET)
# Need to make TET DSB
# Rename TET assay
seurat_data[["CLR_TET"]] <- seurat_data[["TET"]]

# Create new TET assay with dsb normalization
dsb_assay <- GetAssayData(seurat_data, slot = "data", assay = "DSB_ADT")

all_tetramers <- dsb_assay[grepl("tet", rownames(dsb_assay)),]

seurat_data[["DSB_TET"]] <- CreateAssayObject(data = all_tetramers)

# Cluster demux ----------------------------------------------------------------

# Now run through the UMAP for each subset, save the plots
umap_res <- adt_dim_red(seurat_object = seurat_data, assay = "DSB_ADT",
                        reduction_name = "adtdsb", sample_name = sample,
                        resolution = 0.8)

seurat_data <- umap_res$object

pdf(file.path(image_dir, "adt_dsb_plots.pdf"))
print(umap_res$plots)
grid::grid.newpage()
print(umap_res$heatmap1)
grid::grid.newpage()
print(umap_res$heatmap2)
dev.off()


umap_res <- adt_dim_red(seurat_object = seurat_data, assay = "CLR_ADT",
                        reduction_name = "adtclr", sample_name = sample,
                        resolution = 0.8)

seurat_data <- umap_res$object

pdf(file.path(image_dir, "adt_clr_plots.pdf"))
print(umap_res$plots)
grid::grid.newpage()
print(umap_res$heatmap1)
grid::grid.newpage()
print(umap_res$heatmap2)
dev.off()

umap_res <- adt_dim_red(seurat_object = seurat_data, assay = "SCAR_ADT",
                        reduction_name = "adtscar", sample_name = sample,
                        resolution = 0.8, hash_id_col = "scar_hash_id")

seurat_data <- umap_res$object

pdf(file.path(image_dir, "adt_scar_plots.pdf"))
print(umap_res$plots)
grid::grid.newpage()
print(umap_res$heatmap1)
grid::grid.newpage()
print(umap_res$heatmap2)
dev.off()

umap_res <- adt_dim_red(seurat_object = seurat_data, assay = "DSB_TET",
                        reduction_name = "tetdsb", sample_name = sample,
                        resolution = 0.8)

seurat_data <- umap_res$object

pdf(file.path(image_dir, "tet_dsb_plots.pdf"))
print(umap_res$plots)
grid::grid.newpage()
print(umap_res$heatmap1)
grid::grid.newpage()
print(umap_res$heatmap2)
dev.off()

umap_res <- adt_dim_red(seurat_object = seurat_data, assay = "CLR_TET",
                        reduction_name = "tetclr", sample_name = sample,
                        resolution = 0.8)

seurat_data <- umap_res$object

pdf(file.path(image_dir, "tet_clr_plots.pdf"))
print(umap_res$plots)
grid::grid.newpage()
print(umap_res$heatmap1)
grid::grid.newpage()
print(umap_res$heatmap2)
dev.off()

umap_res <- adt_dim_red(seurat_object = seurat_data, assay = "SCAR_TET",
                        reduction_name = "tetscar", sample_name = sample,
                        resolution = 0.8, hash_id_col = "scar_hash_id")

seurat_data <- umap_res$object

pdf(file.path(image_dir, "tet_scar_plots.pdf"))
print(umap_res$plots)
grid::grid.newpage()
print(umap_res$heatmap1)
grid::grid.newpage()
print(umap_res$heatmap2)
dev.off()

umap_res <- adt_dim_red(seurat_object = seurat_data, assay = "SCAR_TET_LOG",
                        reduction_name = "tetscarlog", sample_name = sample,
                        resolution = 0.8, hash_id_col = "scar_hash_id")

seurat_data <- umap_res$object

pdf(file.path(image_dir, "tet_scar_log_plots.pdf"))
print(umap_res$plots)
grid::grid.newpage()
print(umap_res$heatmap1)
grid::grid.newpage()
print(umap_res$heatmap2)
dev.off()

# Try again without doublets, only on the DSB
seurat_no_doub <- subset(seurat_data, subset = tet_hash_id != "Doublet")

umap_res <- adt_dim_red(seurat_object = seurat_no_doub, assay = "DSB_TET",
                        reduction_name = "tetdsb_nd", sample_name = sample,
                        resolution = 0.8)

seurat_no_doub <- umap_res$object

new_meta <- seurat_no_doub[[]] %>%
  dplyr::select(tetdsb_nd_clusters)

seurat_data <- AddMetaData(seurat_data, metadata = new_meta)

pdf(file.path(image_dir, "tet_dsb_nd_plots.pdf"))
print(umap_res$plots)
grid::grid.newpage()
print(umap_res$heatmap1)
grid::grid.newpage()
print(umap_res$heatmap2)
dev.off()

# Plot hto demux as well -------------------------------------------------------
violins <- featDistPlot(seurat_object = seurat_data, 
                       geneset = rownames(all_tetramers),
                       combine = FALSE, sep_by = "tet_hash_id",
                       assay = "DSB_TET")

pdf(file.path(image_dir, "hto_demux_dsb.pdf"))
print(violins)
dev.off()

violins <- featDistPlot(seurat_object = seurat_data, 
                        geneset = rownames(all_tetramers),
                        combine = FALSE, sep_by = "tet_hash_id",
                        assay = "CLR_TET")

pdf(file.path(image_dir, "hto_demux_clr.pdf"))
print(violins)
dev.off()

violins <- featDistPlot(seurat_object = seurat_data, 
                        geneset = rownames(all_tetramers),
                        combine = FALSE, sep_by = "tet_hash_id",
                        assay = "SCAR_TET")

pdf(file.path(image_dir, "hto_demux_scar.pdf"))
print(violins)
dev.off()

violins <- featDistPlot(seurat_object = seurat_data, 
                        geneset = rownames(all_tetramers),
                        combine = FALSE, sep_by = "tet_hash_id",
                        assay = "SCAR_TET_LOG")

pdf(file.path(image_dir, "hto_demux_scar_log.pdf"))
print(violins)
dev.off()

violins <- featDistPlot(seurat_object = seurat_data, 
                        geneset = rownames(all_tetramers),
                        combine = FALSE, sep_by = "scar_hash_id",
                        assay = "DSB_TET")

pdf(file.path(image_dir, "scar_hto_demux_dsb.pdf"))
print(violins)
dev.off()

violins <- featDistPlot(seurat_object = seurat_data, 
                        geneset = rownames(all_tetramers),
                        combine = FALSE, sep_by = "scar_hash_id",
                        assay = "CLR_TET")

pdf(file.path(image_dir, "scar_hto_demux_clr.pdf"))
print(violins)
dev.off()

violins <- featDistPlot(seurat_object = seurat_data, 
                        geneset = rownames(all_tetramers),
                        combine = FALSE, sep_by = "scar_hash_id",
                        assay = "SCAR_TET")

pdf(file.path(image_dir, "scar_hto_demux_scar.pdf"))
print(violins)
dev.off()

violins <- featDistPlot(seurat_object = seurat_data, 
                        geneset = rownames(all_tetramers),
                        combine = FALSE, sep_by = "scar_hash_id",
                        assay = "SCAR_TET_LOG")

pdf(file.path(image_dir, "scar_hto_demux_scar_log.pdf"))
print(violins)
dev.off()


saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))