library(here)
library(scAnalysisR)
library(pheatmap)
library(tidyverse)
library(splitstackshape)
library(circlize)
library(viridis)
library(plotly)
library(ggalluvial)

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


# Pull out v_gene information
sep_columns <- c("chains", "cdr3", "cdr3_length",
                 "cdr3_nt_length", "v_gene", "d_gene", "j_gene", "c_gene",
                 "reads", "umis", "productive", "full_length",
                 "v_ins", "v_del", "v_mis", "d_ins", "d_del",
                 "d_mis", "j_ins", "j_del", "j_mis", "c_ins", "c_del", "c_mis",
                 "all_ins", "all_del", "all_mis", "vd_ins", "vd_del", "dj_ins",
                 "dj_del", "v_mis_freq", "d_mis_freq", "j_mis_freq",
                 "c_mis_freq", "all_mis_freq")
keep_columns <- c("isotype", "RNA_combined_celltype", "sample", "paired",
                  "clonotype_id", "Status", "tet_hash_id", "all_chains")

all_info <- seurat_data[[]] %>%
  dplyr::mutate(all_chains = chains) %>%
  dplyr::select(dplyr::all_of(c(sep_columns, keep_columns)), n_chains) %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::filter(!is.na(n_chains))


all_info_split <- cSplit(all_info, sep_columns, sep = ";", direction = "long") %>%
  dplyr::filter(!is.na(chains))


# Things to pull out
# What is the most common J and light chain associated with v4-4
# Most common J and light chain associated with any odds ratio v gene
# Subset to just class switched with 3+ SHM
  # Find percent of all cells within a clone that bind insulin/gad/ia2
  # Find number SHM
  # Order by both percent ins/gad/ia2 and number SHM

# Odds ratio v genes
# Read in odds results
v_counting_dir <- file.path(save_dir, "files", "v_gene_counting")

ifelse(!dir.exists(v_counting_dir), dir.create(v_counting_dir), FALSE)

tests_use <- c("no_nd", "aab1_nd", "aab2_nd")
all_v_genes <- lapply(tests_use, function(x){
  
  # Keep only genes that were used significantly more highly in 
  # the condition vs control
  odds_data <- openxlsx::read.xlsx(file.path(save_dir, "files", "vdj_files",
                                             "stats", paste0(x, ".xlsx")),
                                   sheet = "odds_ratio") %>%
    dplyr::filter(p_adj < 0.05, odds_ratio > 1) %>%
    dplyr::select(v_gene, test, odds_ratio)
  
  return(odds_data)
  
})

all_v_genes <- do.call(rbind, all_v_genes)

wide_v_genes <- all_v_genes %>%
  dplyr::mutate(test = paste0("odds_ratio_", test)) %>%
  tidyr::pivot_wider(names_from = "test", values_from = "odds_ratio")

save_wb <- openxlsx::createWorkbook()

count_v_j_gene <- lapply(unique(all_v_genes$v_gene), function(x){
  short_info <- all_info_split %>%
    dplyr::filter(v_gene == x) %>%
    dplyr::group_by(j_gene) %>%
    dplyr::add_count(name = "j_frequency") %>%
    dplyr::group_by(j_gene, isotype) %>%
    dplyr::add_count(name = "j_gene_isotype_frequency") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(j_gene, RNA_combined_celltype) %>%
    dplyr::add_count(name = "j_celltype_frequency") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(j_gene, tet_hash_id) %>%
    dplyr::add_count(name = "j_tetramer_frequency") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(j_gene, Status) %>%
    dplyr::add_count(name = "j_status_frequency") %>%
    dplyr::select(v_gene, j_gene, isotype, RNA_combined_celltype, tet_hash_id,
                  Status, dplyr::contains("frequency")) %>%
    dplyr::distinct() %>%
    dplyr::arrange(desc(j_frequency), desc(j_tetramer_frequency),
                   desc(j_status_frequency),
                   desc(j_celltype_frequency),
                   desc(j_gene_isotype_frequency))
  
  short_odds <- wide_v_genes %>%
    dplyr::filter(v_gene == x) %>%
    dplyr::select(-v_gene)
  
  short_info <- cbind(short_info, short_odds)
  
  openxlsx::addWorksheet(wb = save_wb, sheetName = x)
  openxlsx::writeData(wb = save_wb, sheet = x, x = short_info)
  
  return(short_info)
})

openxlsx::saveWorkbook(save_wb, file.path(v_counting_dir,
                                          "v_j_frequency_high_odds_ratio.xlsx"), 
                       overwrite = TRUE)

keep_chains <- c("IGH;IGK", "IGH;IGL")

all_info_split_hl <- all_info_split %>%
  dplyr::filter(all_chains %in% keep_chains) %>%
  dplyr::select(isotype, RNA_combined_celltype, tet_hash_id, Status,
                chains, barcode, v_gene) %>%
  tidyr::pivot_wider(names_from = chains, values_from = v_gene) %>%
  dplyr::mutate(IGK_IGL = ifelse(!is.na(IGK), IGK, IGL)) %>%
  dplyr::select(IGH, IGK_IGL, isotype, RNA_combined_celltype, tet_hash_id,
                Status)

save_wb <- openxlsx::createWorkbook()

count_v_hl_gene <- lapply(unique(all_v_genes$v_gene), function(x){
  short_info <- all_info_split_hl %>%
    dplyr::filter(IGH == x) %>%
    dplyr::group_by(IGK_IGL) %>%
    dplyr::add_count(name = "light_frequency") %>%
    dplyr::group_by(IGK_IGL, isotype) %>%
    dplyr::add_count(name = "light_gene_isotype_frequency") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(IGK_IGL, RNA_combined_celltype) %>%
    dplyr::add_count(name = "light_celltype_frequency") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(IGK_IGL, tet_hash_id) %>%
    dplyr::add_count(name = "light_tetramer_frequency") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(IGK_IGL, Status) %>%
    dplyr::add_count(name = "light_status_frequency") %>%
    dplyr::select(IGH, IGK_IGL, isotype, RNA_combined_celltype, tet_hash_id,
                  Status,dplyr::contains("frequency")) %>%
    dplyr::distinct() %>%
    dplyr::arrange(desc(light_frequency), desc(light_tetramer_frequency),
                   desc(light_status_frequency),
                   desc(light_celltype_frequency), 
                   desc(light_gene_isotype_frequency))
  
  short_odds <- wide_v_genes %>%
    dplyr::filter(v_gene == x) %>%
    dplyr::select(-v_gene)
  
  short_info <- cbind(short_info, short_odds)
  
  openxlsx::addWorksheet(wb = save_wb, sheetName = x)
  openxlsx::writeData(wb = save_wb, sheet = x, x = short_info)
  
  return(short_info)
})

openxlsx::saveWorkbook(save_wb, file.path(v_counting_dir,
                                          "heavy_light_frequency_high_odds_ratio.xlsx"), 
                       overwrite = TRUE)
# Plots
# CDR3 length as density plot
# CDR3 length subset to just Class switched with 3+ SHM
# Percent SHM

# CDR3 length as density plot

