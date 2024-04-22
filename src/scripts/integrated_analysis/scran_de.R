# DE scran
library(here)
library(scAnalysisR)
library(pheatmap)
library(tidyverse)
library(scran)
library(Seurat)

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

DefaultAssay(seurat_data) <- "RNA"

seurat_data <- DietSeurat(object = seurat_data, assays = "RNA")

test_mapping <- c("DNA.tet" = "Other",
                  "GAD.tet" = "Islet_Reactive",
                  "IA2.tet" = "Islet_Reactive",
                  "INS.tet" = "Islet_Reactive",
                  "Islet_Multi_Reactive" = "Islet_Reactive",
                  "Negative" = "Other",
                  "Other_Multi_Reactive" = "Other",
                  "TET.tet" = "Other")

seurat_data$test_id <- test_mapping[seurat_data$tet_name_cutoff]

sce <- as.SingleCellExperiment(seurat_data)

sce$full_test <- paste(sce$Status, sce$test_id, sep = "_")

pairing_data <- data.frame("A" = c("T1D_Islet_Reactive",
                                   "AAB_Islet_Reactive"),
                           "B" = c("ND_Islet_Reactive",
                                   "ND_Islet_Reactive"))

all_res <- scoreMarkers(sce, groups = sce$full_test, 
                        block = sce$date.processed.for.scSeq,
                        pairings = pairing_data)


t1d <- data.frame(all_res$T1D_Islet_Reactive)

ggplot2::ggplot(t1d, ggplot2::aes(median.AUC)) +
  ggplot2::geom_histogram(binwidth = 0.005) +
  ggplot2::xlim(0.4, 0.6)

ggplot2::ggplot(t1d, ggplot2::aes(median.logFC.cohen)) +
  ggplot2::geom_histogram(binwidth = 0.005) +
  ggplot2::xlim(-0.5, 0.5)

ggplot2::ggplot(t1d, ggplot2::aes(median.logFC.detected)) +
  ggplot2::geom_histogram(binwidth = 0.05) +
  ggplot2::xlim(-2, 2)

aab <- data.frame(all_res$AAB_Islet_Reactive)

ggplot2::ggplot(aab, ggplot2::aes(median.AUC)) +
  ggplot2::geom_histogram(binwidth = 0.005) +
  ggplot2::xlim(0.4, 0.6)

ggplot2::ggplot(aab, ggplot2::aes(median.logFC.cohen)) +
  ggplot2::geom_histogram(binwidth = 0.005) +
  ggplot2::xlim(-0.5, 0.5)


ggplot2::ggplot(aab, ggplot2::aes(median.logFC.detected)) +
  ggplot2::geom_histogram(binwidth = 0.05) +
  ggplot2::xlim(-2, 2)



cutoff_percentile <- quantile(t1d$median.AUC, 0.95)

t1d_sig <- t1d[t1d$median.AUC > cutoff_percentile,]

t1d_sig <- t1d_sig[t1d_sig$median.logFC.detected > 0.25,]

cutoff_percentile <- quantile(aab$median.AUC, 0.95)

aab_sig <- aab[aab$median.AUC > cutoff_percentile,]

aab_sig <- aab_sig[aab_sig$median.logFC.detected > 0.25,]

