library(here)
library(scAnalysisR)
library(Seurat)
library(tidyverse)

# Naive - resting memory/activate naive/intermediate - ABC/memory - Plasmablast

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

# Set directories
save_dir <- file.path(results_dir, "R_analysis", sample)

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj",
                                 "seurat_processed_no_doublet.rds"))


# Make geo directory
ifelse(!dir.exists(here("GEO_submission")), dir.create(here("GEO_submission")),
       FALSE)

# Make file to move and rename fastq -------------------------------------------

all_meta <- seurat_data[[]] %>%
  dplyr::select(ID, Sample.Name) %>%
  dplyr::distinct()

rownames(all_meta) <- NULL

write.csv(all_meta, here("files", "meta_mapping.csv"), quote = FALSE)

rm(all_meta)

# Will need (python):
# All fastq files renamed by the sample id
# All features, barcodes, and matrix files
# filtered_contig_annotations.csv.gz, concat_ref.bam, airr_rearrangement.tsv.gz for vdj
# Tar all of these

# Make all meta data file ------------------------------------------------------

select_cols <- c("nCount_RNA", "nFeature_RNA", "nCount_ADT", "nFeature_ADT",
                 "percent.mt", "Sex", "Collection.Date",
                 "Age.at.Collection..years.", "Age.at.Collection..Months.",
                 "Status", "Date.of.Diagnosis", "Days.post.onset",
                 "millions.of.cells.frozen", "HLA.type", "Autoantibodies",
                 "date.processed.for.scSeq", "Notes..FDR.relationship.",
                 "sample", "Phase", "chains", "n_chains", "cdr3", "cdr3_nt",
                 "cdr3_length", "cdr3_nt_length", "v_gene", "d_gene", "j_gene",
                 "c_gene", "isotype", "reads", "umis", "productive",
                 "full_length", "paired", "all_mis_freq", "libra_tet_hash_id",
                 "libra_full_hash_id", "scar_libra_tet_hash_id",
                 "scar_libra_full_hash_id", "tet_name_cutoff",
                 "full_tet_name_cutoff", "RNA_cluster", "cluster_celltype",
                 "final_clone", "imcantation_isotype")
                 
pub_meta <- seurat_data[[]] %>%
  dplyr::select(dplyr::all_of(select_cols))

pub_meta$Autoantibodies <- NULL
pub_meta$HLA.type <- NULL

# Update autoantibodies, HLA, ethnicity
manuscript_meta <- openxlsx::read.xlsx(here("files/manuscript_meta.xlsx"))

manuscript_meta$sample <- manuscript_meta$ID
manuscript_meta <- manuscript_meta[, c("sample", "HLA", "Autoantibodies", 
                                       "Ethnicity")]

pub_meta <- merge(pub_meta, manuscript_meta, by = "sample")


write.csv(pub_meta, here("GEO_submission", "metadata.csv"))

rm(pub_meta)

# Make all raw and nomralized counts files -------------------------------------

# Counts we need
# RNA raw and normalized for merged
raw_RNA <- GetAssayData(seurat_data, assay = "RNA", slot = "counts")
saveRDS(raw_RNA, here("GEO_submission", "raw_RNA_counts.rda"))
rm(raw_RNA)

normalized_RNA <- GetAssayData(seurat_data, assay = "RNA", slot = "data")
saveRDS(normalized_RNA, here("GEO_submission", "normalized_RNA_counts.rda"))
rm(normalized_RNA)

# ADT raw, normalized, libra, and quantile score for merged
raw_ADT <- GetAssayData(seurat_data, assay = "ADT", slot = "counts")
scar_ADT <- GetAssayData(seurat_data, assay = "SCAR_ADT", slot = "counts")

raw_ADT <- raw_ADT[!rownames(raw_ADT) %in% c("HTO1", "HTO3"),]
scar_ADT <- scar_ADT[!rownames(scar_ADT) %in% c("HTO1", "HTO3"),]

saveRDS(raw_ADT, here("GEO_submission", "raw_ADT.rda"))
saveRDS(scar_ADT, here("GEO_submission", "scar_ADT.rda"))

rm(raw_ADT)
rm(scar_ADT)

# TET_LIBRA
libra_score <- GetAssayData(seurat_data, assay = "TET_LIBRA", slot = "data")

# NEW_TET_PROPORTIONS
quantile_score <- GetAssayData(seurat_data, assay = "NEW_TET_PROPORTIONS",
                               slot = "data")


saveRDS(libra_score, here("GEO_submission", "libra_score.rda"))
saveRDS(quantile_score, here("GEO_submission", "quantile_score.rda"))

# PCA
pca <- Embeddings(seurat_data, reduction = "pca")

# MNN
mnn <- Embeddings(seurat_data, reduction = "rna_mnn")

# Corrected UMAP
umap <- Embeddings(seurat_data, reduction = "rna_mnn.umap")

# write.csv(pca, here("GEO_submission", "pca.csv"))
# write.csv(mnn, here("GEO_submission", "mnn.csv"))
write.csv(umap, here("GEO_submission", "umap.csv"))
