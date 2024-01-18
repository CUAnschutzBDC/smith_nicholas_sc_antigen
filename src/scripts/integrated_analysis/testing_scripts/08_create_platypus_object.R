library(here)
library(scAnalysisR)
library(pheatmap)
library(tidyverse)
library(Platypus)
library(splitstackshape)

normalization_method <- "log" # can be SCT or log
# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

#sample <- args[[1]]
#sample <- gsub("__.*", "", sample)
sample <- "merged"

#results_dir <- args[[2]]
results_dir <- here("results")

#sample_info <- args[[4]]
sample_info <- here("files/sample_info.tsv")

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

ifelse(!dir.exists(file.path(save_dir, "images", "vdj_analysis")),
       dir.create(file.path(save_dir, "images", "vdj_analysis")),
       FALSE)

# Colors -----------------------------------------------------------------------
all_colors <- readRDS(file = file.path("files/all_colors.rds"))


final_colors <- all_colors$cell_type_colors

tetramer_colors <- all_colors$tetramer_colors

sample_colors <- all_colors$sample_colors


status_colors <- all_colors$status_colors

celltypes_keep <- c("Naive_1", "Naive_3",
                    "BND2", "B.intermediate",
                    "Memory_IgA", "Resting_memory", 
                    "Activated_memory",
                    "Plasmablast")

seurat_data <- subset(seurat_data,
                      subset = RNA_combined_celltype %in% celltypes_keep)

# Read in data with platypus ---------------------------------------------------
all_samples <- unique(seurat_data$sample)

vdj_out_directory <- lapply(all_samples, function(x){
  return(here("results", x, "outs", "per_sample_outs",
              x, "vdj_b"))
})

gex_out_directory <- lapply(all_samples, function(x){
  return(here("results", x, "outs", "per_sample_outs",
              x, "count"))
})

names(gex_out_directory) <- all_samples
names(vdj_out_directory) <- all_samples

# This takes a very long time to run
platypus_all <- VDJ_GEX_matrix(VDJ.out.directory.list = vdj_out_directory,
                               GEX.out.directory.list = gex_out_directory,
                               GEX.integrate = T,
                               VDJ.combine = T,
                               integrate.GEX.to.VDJ = T,
                               integrate.VDJ.to.GEX = T)


saveRDS(platypus_all,
        file.path(save_dir, "rda_obj", "platypus_obj.rds"))

vdj_data <- platypus_all$VDJ

# Map sample name back to my sample name
all_samples <- unique(seurat_data$sample)
names(all_samples) <- paste0("s", rep(1:length(all_samples)))
vdj_data$sample <- all_samples[vdj_data$sample_id]

vdj_data$sample_barcode <- gsub("s[0-9]+_", "", vdj_data$barcode)
vdj_data$sample_barcode <- paste(vdj_data$sample, vdj_data$sample_barcode,
                                 sep = "_")

meta_data <- seurat_data[[]] %>%
  dplyr::select(sample, Status, tet_hash_id, cdr3, n_chains, isotype, v_gene,
                c_gene, j_gene, d_gene) %>%
  tibble::rownames_to_column("barcode") %>% 
  dplyr::mutate(sample_barcode = gsub("-.*", "", barcode)) %>%
  dplyr::mutate(sample_barcode = paste(sample, sample_barcode, sep = "_")) %>%
  dplyr::select(-barcode, -sample)

colnames(meta_data) <- c("Status", "tet_hash_id",
                         "tenx_cdr3", "tenx_nchains", "tenx_isotype",
                         "tenx_vgene", "tenx_cgene",
                         "tenx_jgene", "tenx_dgene", "sample_barcode")

# Subset to only cells that are in the seurat object
use_vdj <- vdj_data %>%
  dplyr::filter(sample_barcode %in% meta_data$sample_barcode)

# Needed to install stringdist
new_clones <- VDJ_clonotype(VDJ = use_vdj,
                            VDJ.VJ.1chain = F, #Keeping cells with aberrant chain numbers
                            clone.strategy = "VDJcdr3.homology",
                            homology.threshold = 0.2, # corresponds to 80% sequence similarity
                            global.clonotype = FALSE) # Clonotype can't occur across samples  

new_clones_sample <- VDJ_clonotype(VDJ = use_vdj,
                                   VDJ.VJ.1chain = F, #Keeping cells with aberrant chain numbers
                                   clone.strategy = "VDJcdr3.homology",
                                   homology.threshold = 0.2, # corresponds to 80% sequence similarity
                                   global.clonotype = TRUE) # Clonotype can occur across samples  

save_list <- list("within_samples" = new_clones,
                  "across_samples" = new_clones_sample)


saveRDS(save_list,
        file.path(save_dir, "rda_obj", "platypus_clones.rds"))
