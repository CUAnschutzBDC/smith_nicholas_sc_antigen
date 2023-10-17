# Use this to figure out how immcantation is working

library(here)
library(scAnalysisR)
library(pheatmap)
library(tidyverse)
library(alakazam)
library(shazam)

normalization_method <- "log" # can be SCT or log
# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

#sample <- args[[1]]
#sample <- gsub("__.*", "", sample)
sample <- "JH310-12_AP"

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

immcantation_path <- file.path(here("results", sample, "outs", 
                                    "immcantation", "filtered_contig_igblast_db-pass.tsv"))
  
vdj_data <- alakazam::readChangeoDb(immcantation_path)

tenx_data <- read.csv(file.path("results", sample, "outs",
                                "per_sample_outs", sample, "vdj_b",
                                "filtered_contig_annotations.csv"))
  

tenx_data$sequence_id <- tenx_data$contig_id

merged_data <- merge(tenx_data, vdj_data, by = "sequence_id",
                     all.x = TRUE, all.y = TRUE)



# Lots of NAs at this step
seurat_meta <- seurat_data[[]] %>%
  dplyr::select(RNA_celltype, chains, n_chains,
                cdr3, v_gene, d_gene, j_gene, isotype)

colnames(seurat_meta) <- paste0("seurat_", colnames(seurat_meta))

seurat_meta <- seurat_meta %>%
  tibble::rownames_to_column("barcode")

tenx_small <- tenx_data %>%
  dplyr::select(v_gene, d_gene, j_gene, c_gene, cdr3, barcode)

merge_with_seurat <- merge(tenx_small, seurat_meta, by = "barcode",
                           all.x = FALSE, all.y = TRUE)


# Now 100% of these are NA for the 10x data, so anything that isn't
# included in the output from 10x is given the NA values. So places
# where there is gene expression but no VDJ
na_vals <- merge_with_seurat[is.na(merge_with_seurat$seurat_v_gene),]


vdj_small <- vdj_data %>%
  dplyr::select(v_call, d_call, j_call, cell_id, cdr3) %>%
  dplyr::rename(barcode = cell_id)


merge_immcantation_seurat <- merge(vdj_small, seurat_meta, by = "barcode",
                                   all.x = FALSE, all.y = TRUE)


# So immcantation also just doesn't have any output when there is no v gene
# detected.
na_vals <- merge_immcantation_seurat[is.na(merge_immcantation_seurat$seurat_v_gene),]
