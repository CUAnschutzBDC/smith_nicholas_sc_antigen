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

# The clones com from Chang-o define clones
# based only on heavy or light
# Based on the distance found in 07_run_immcantation using shazam
# Note, these appear to be found within samples individually
clone_info <- read.table(file.path(save_dir, "define_clones",
                                   "immcantation_combined_clone-pass.tsv"),
                         sep = "\t", header = TRUE)

sample_clone_info <- read.table(file.path(save_dir, "define_clones",
                                   "immcantation_combined_clone-pass_sample.tsv"),
                         sep = "\t", header = TRUE)

sample_clone_info <- sample_clone_info %>%
  dplyr::select(sequence_id, clone_id, locus, sample) %>%
  dplyr::rename(sample_clone_id = clone_id,
                sample_locus = locus) %>%
  dplyr::mutate(sequence_id_sample = paste(sequence_id,
                                             sample, sep = "_")) %>%
  dplyr::select(-sequence_id, -sample)

clone_info <- clone_info %>%
  dplyr::mutate(sequence_id_sample = paste(sequence_id,
                                             sample, sep = "_")) %>%
  merge(sample_clone_info, by = "sequence_id_sample", all.x = TRUE)

# Find expanded clones
clone_info_expanded <- clone_info %>%
  dplyr::group_by(clone_id) %>%
  dplyr::add_count(name = "clone_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sample_clone_id) %>%
  dplyr::add_count(name = "sample_clone_count") %>%
  dplyr::ungroup()

# This just makes sure that heavy and light are always in different clones
# They are so it doesn't need to be run
# heavy_light <- clone_info %>%
#   dplyr::select(clone_id, locus) %>%
#   dplyr::ungroup() %>%
#   dplyr::distinct() %>%
#   dplyr::group_by(clone_id) %>%
#   dplyr::add_count(name = "locus_count") 

# See how the calling works
expanded <- clone_info_expanded %>%
  dplyr::filter(locus == "IGH")

# Seems like the V, J and junciton are all similar/the same
expanded <- expanded %>%
  dplyr::filter(clone_count == 107)

expanded_2 <- clone_info_expanded %>%
  dplyr::filter(locus == "IGH", clone_count > 1)

# Add heavy and light clones separately to the seurat object
# First make a barcode_sample column, do the same with seurat
# metadata, merge by this, and rename the rownames as the barcodes
# from seurat. Then add to so.

# Find cases of shared clones across samples
sample_clones <- clone_info %>%
  dplyr::filter(locus == "IGH") %>%
  dplyr::group_by(sample_clone_id) %>%
  dplyr::add_count(name = "sample_clone_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(clone_id) %>%
  dplyr::add_count(name= "clone_count") %>%
  dplyr::ungroup() %>%
  dplyr::select(sample, clone_id, sample_clone_count, clone_count) %>%
  dplyr::distinct() %>%
  dplyr::group_by(clone_id) %>%
  dplyr::add_count(name = "sample_count") %>%
  dplyr::filter(sample_count > 1)
