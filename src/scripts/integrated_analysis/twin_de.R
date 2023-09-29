library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log


cell_types <- "RNA_celltype"
clusters <- "RNA_cluster"
pval <- 0.05
logfc <- 0.5

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
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))

# Colors -----------------------------------------------------------------------
final_colors <- c("Resting_memory" = "#924bdb", # Resting memory
                  "Naive_1" = "#69ba3d", # Naive 1
                  "Naive_2" = "#9a43a4", # Naive 2
                  "Memory_IgE_IgG" = "#bf9b31", # Memory IgE/IgG1
                  "Naive_3" = "#6477ce", # Naive 3
                  "Memory_IgA" = "#d15131", # Memory IA
                  "Early_memory" = "#4c9e8e", # Early Memory
                  "BND2" = "#cc4570", #Bnd2
                  "DN2" = "#648d4f", # DN2
                  "Activated_memory" = "#985978", # Activated memory
                  "Activated_naive" = "#a06846", # Activated naive
                  "B.intermediate" = "#00008b",
                  "CD14.Mono" = "#e0205a",
                  "pDC" = "#ffb6d3",
                  "Plasmablast" = "#ffac14",
                  "CD8.TEM" = "#000000") 

tetramer_colors <- MetBrewer::met.brewer(name = "Juarez", n = 7,
                                         type = "continuous")
names(tetramer_colors) <- c("Negative", "Doublet", "INS-tet", "TET-tet",
                            "IA2-tet", "GAD-tet", "DNA-tet")

sample_colors <- MetBrewer::met.brewer(name = "Archambault", n = 16,
                                       type = "continuous")


all_samples <- unique(seurat_data$sample)

names(sample_colors) <- all_samples

antigen_colors <- MetBrewer::met.brewer(name = "Demuth", n = 4,
                                        type = "discrete")

names(antigen_colors) <- c("Negative", "Tet_antigen", "other",
                           "diabetes_antigen")


status_colors <- MetBrewer::met.brewer(name = "Egypt", n = 4,
                                       type = "discrete")

seurat_data$Status <- gsub(" ", "_", seurat_data$Status)

names(status_colors) <- c("nd", "no", "aab_stage_1", "aab_stage_2")

# DE ---------------------------------------------------------------------------
memory_cells <- c("Activated_memory", "Memory_IgA", 
                  "Plasmablast", "Resting_memory")

naive_cells <- c("Naive_1", "Naive_2", "Naive_3")

all_b <- c(memory_cells, naive_cells, "B.intermediate", "BND2", "Plasmablast")

seurat_data_mem <- subset(seurat_data, 
                          subset = RNA_combined_celltype %in% memory_cells)

seurat_data_naive <- subset(seurat_data,
                            subset = RNA_combined_celltype %in% naive_cells)

seurat_data_b <- subset(seurat_data,
                        subset = RNA_combined_celltype %in% all_b)

Idents(seurat_data_mem) <- "sample"
Idents(seurat_data_naive) <- "sample"
Idents(seurat_data_b) <- "sample"

DE_genes_mem <- FindMarkers(seurat_data_mem, ident.1 = "JH310-12_FG",
                            ident.2 = "JH310-12_NG")

DE_genes_niave <- FindMarkers(seurat_data_naive, ident.1 = "JH310-12_FG",
                               ident.2 = "JH310-12_NG")

DE_genes_b <- FindMarkers(seurat_data_b, ident.1 = "JH310-12_FG",
                            ident.2 = "JH310-12_NG")

# Check against others

processed_date <- read.csv(file.path(save_dir, "files",
                                     "pseudobulk_all_bcell_de_processed_date.csv"))

collected_date_bin <- read.csv(file.path(save_dir, "files",
                                         "pseudobulk_all_bcell_de.csv"))

capture_date <- read.csv(file.path(save_dir, "files",
                                   "pseudobulk_all_bcell_de_capture_correct.csv"))


previous_res <- read.csv("/beevol/home/wellskri/Analysis/Mia_Smith/Catherine_Nicolas/20230224_BND_HC_NPOD_pLN_spl_3_NO_T1D_1_FDR/results/Catherine_combined/R_analysis/files/pseudobulk_all_bcell_de.csv")

processed_date_nd_no <- processed_date %>%
  dplyr::filter(contrast == "NO - ND")


collection_date_nd_no <- collected_date_bin %>%
  dplyr::filter(contrast == "NO - ND")

previous_res_nd_no <- previous_res %>%
  dplyr::filter(contrast == "NO - ND")

capture_res_nd_no <- capture_date %>%
  dplyr::filter(contrast == "NO - ND")

DE_genes_b_sig <- DE_genes_b %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  tibble::rownames_to_column("gene")

dim(processed_date_nd_no)
dim(collection_date_nd_no)
dim(previous_res_nd_no)
dim(capture_res_nd_no)
dim(DE_genes_b_sig)
length(intersect(processed_date_nd_no$gene, DE_genes_b_sig$gene))
length(intersect(collection_date_nd_no$gene, DE_genes_b_sig$gene))
length(intersect(capture_res_nd_no$gene, DE_genes_b_sig$gene))
length(intersect(previous_res_nd_no$gene, DE_genes_b_sig$gene))



# Heatmap of DE genes
# Change sample order
sample_order <- c("JH310-12_AP",
                  "JH310-12_TS",
                  "JH313-15_OF",
                  "JH310-12_MV",
                  "JH310-12_NG",
                  "JH313-15_DB",
                  "JH313-15_MH",
                  "JH310-12_BH",
                  "JH313-15_AG",
                  "JH313-15_NF",
                  "JH313-15_PH",
                  "JH310-12_FG",
                  "JH310-12_EB",
                  "JH310-12_CP",
                  "JH313-15_FD",
                  "JH313-15_JB")

antigen_order <- c("diabetes_antigen", "other", "Negative", "Tet_antigen")

# Any diabetes antigen and no binding
name_mapping <- c("INS-tet" = "diabetes_antigen",
                  "GAD-tet" = "diabetes_antigen",
                  "IA2-tet" = "diabetes_antigen",
                  "TET-tet" = "Tet_antigen",
                  "Negative" = "Negative",
                  "Doublet" = "other",
                  "DNA-tet" = "other")

seurat_data$test_id <- name_mapping[as.character(seurat_data$hash.ID)]


seurat_data$meta_col <- paste(seurat_data$sample, seurat_data$test_id,
                                sep = "_")
meta_df <- seurat_data[[c("sample", "test_id", "Status", "meta_col")]]
meta_ave <- meta_df
rownames(meta_ave) <- NULL
meta_ave <- distinct(meta_ave)
rownames(meta_ave) <- meta_ave$meta_col

meta_ave$sample <- factor(meta_ave$sample, levels = sample_order)
meta_ave$test_id <- factor(meta_ave$test_id, levels = antigen_order)

meta_ave <- meta_ave %>%
  dplyr::arrange(sample, test_id)

meta_ave$Status <- gsub(" ", "_", meta_ave$Status)

color_list <- list("sample" = sample_colors[names(sample_colors) %in% 
                                              unique(meta_df$sample)],
                   "test_id" = antigen_colors,
                   "Status" = status_colors)


print(scAnalysisR::plot_heatmap(seurat_object = seurat_data, meta_df = meta_ave,
                                gene_list = DE_genes_b_sig$gene,
                                average_expression = TRUE,
                                meta_col = "meta_col", plot_meta_col = FALSE,
                                max_val = 2.5, min_val = -2.5,
                                cluster_rows = TRUE, cluster_cols = FALSE,
                                color_list = color_list, plot_rownames = FALSE))

