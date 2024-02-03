library(here)
library(scAnalysisR)
library(pheatmap)
library(tidyverse)
library(splitstackshape)
library(circlize)
library(viridis)
library(plotly)
library(ggalluvial)
library(Seurat)

normalization_method <- "log" # can be SCT or log
# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

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

# Get mean and sd of insulin, gad and ia2

assay_data1 <- GetAssayData(seurat_data, slot = "data", assay = "NEW_TET_PROPORTIONS") %>%
  as.matrix() %>%
  t() %>%
  data.frame

colnames(assay_data1) <- paste0("bimodal_cutoff_", colnames(assay_data1))

seurat_data <- AddMetaData(seurat_data, metadata = assay_data1)


assay_data2 <- GetAssayData(seurat_data, slot = "data", assay = "SCAR_TET_LIBRA") %>%
  as.matrix() %>%
  t() %>%
  data.frame

colnames(assay_data2) <- paste0("libra_cutoff_", colnames(assay_data2))


seurat_data <- AddMetaData(seurat_data, metadata = assay_data2)


# Pull out v_gene information
sep_columns <- c("chains", "cdr3", "cdr3_length",
                 "cdr3_nt_length", "v_gene", "d_gene", "j_gene", "c_gene",
                 "reads", "umis", "productive", "full_length",
                 "v_ins", "v_del", "v_mis", "d_ins", "d_del",
                 "d_mis", "j_ins", "j_del", "j_mis", "c_ins", "c_del", "c_mis",
                 "all_ins", "all_del", "all_mis", "vd_ins", "vd_del", "dj_ins",
                 "dj_del", "v_mis_freq", "d_mis_freq", "j_mis_freq",
                 "c_mis_freq", "all_mis_freq")
keep_columns <- c("isotype", "final_celltype", "sample", "paired",
                  "clonotype_id", "Status", "tet_name_cutoff", "full_tet_name_cutoff",
                  "scar_libra_tet_hash_id", "scar_libra_full_hash_id",
                  "all_chains", colnames(assay_data1), colnames(assay_data2))

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

tests_use <- c("all_t1d_nd", "all_aab_nd")
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
    dplyr::mutate("INS_mean" = mean(bimodal_cutoff_INS.tet),
                  "INS_sd" = sd(bimodal_cutoff_INS.tet),
                  "TET_mean" = mean(bimodal_cutoff_TET.tet),
                  "TET_sd" = sd(bimodal_cutoff_TET.tet),
                  "GAD_mean" = mean(bimodal_cutoff_GAD.tet),
                  "GAD_sd" = sd(bimodal_cutoff_GAD.tet),
                  "IA2_mean" = mean(bimodal_cutoff_IA2.tet),
                  "IA2_sd" = sd(bimodal_cutoff_IA2.tet),
                  "DNA_mean" = mean(bimodal_cutoff_DNA.tet),
                  "DNA_sd" = sd(bimodal_cutoff_DNA.tet)) %>%
    dplyr::group_by(j_gene, isotype) %>%
    dplyr::add_count(name = "j_gene_isotype_frequency") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(j_gene, final_celltype) %>%
    dplyr::add_count(name = "j_celltype_frequency") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(j_gene, tet_name_cutoff) %>%
    dplyr::add_count(name = "j_tetramer_frequency") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(j_gene, Status) %>%
    dplyr::add_count(name = "j_status_frequency") %>%
    dplyr::select(v_gene, j_gene, isotype, final_celltype, tet_name_cutoff,
                  full_tet_name_cutoff, Status, dplyr::contains("frequency"), 
                  dplyr::contains("mean"),
                  dplyr::contains("sd")) %>%
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
  dplyr::select(isotype, final_celltype, tet_name_cutoff,
                full_tet_name_cutoff, Status,
                chains, barcode, v_gene, dplyr::contains("tet")) %>%
  tidyr::pivot_wider(names_from = chains, values_from = v_gene) %>%
  dplyr::mutate(IGK_IGL = ifelse(!is.na(IGK), IGK, IGL)) %>%
  dplyr::select(IGH, IGK_IGL, isotype, final_celltype, tet_name_cutoff,
                full_tet_name_cutoff, Status, dplyr::contains("tet"))

save_wb <- openxlsx::createWorkbook()

count_v_hl_gene <- lapply(unique(all_v_genes$v_gene), function(x){
  short_info <- all_info_split_hl %>%
    dplyr::filter(IGH == x) %>%
    dplyr::group_by(IGK_IGL) %>%
    dplyr::add_count(name = "light_frequency") %>%
    dplyr::mutate("INS_mean" = mean(bimodal_cutoff_INS.tet),
                  "INS_sd" = sd(bimodal_cutoff_INS.tet),
                  "TET_mean" = mean(bimodal_cutoff_TET.tet),
                  "TET_sd" = sd(bimodal_cutoff_TET.tet),
                  "GAD_mean" = mean(bimodal_cutoff_GAD.tet),
                  "GAD_sd" = sd(bimodal_cutoff_GAD.tet),
                  "IA2_mean" = mean(bimodal_cutoff_IA2.tet),
                  "IA2_sd" = sd(bimodal_cutoff_IA2.tet),
                  "DNA_mean" = mean(bimodal_cutoff_DNA.tet),
                  "DNA_sd" = sd(bimodal_cutoff_DNA.tet)) %>%
    dplyr::group_by(IGK_IGL, isotype) %>%
    dplyr::add_count(name = "light_gene_isotype_frequency") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(IGK_IGL, final_celltype) %>%
    dplyr::add_count(name = "light_celltype_frequency") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(IGK_IGL, tet_name_cutoff) %>%
    dplyr::add_count(name = "light_tetramer_frequency") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(IGK_IGL, Status) %>%
    dplyr::add_count(name = "light_status_frequency") %>%
    dplyr::select(IGH, IGK_IGL, isotype, final_celltype, tet_name_cutoff,
                  full_tet_name_cutoff, Status,
                  dplyr::contains("frequency"), dplyr::contains("mean"),
                  dplyr::contains("sd")) %>%
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

# Repeat with clones -----------------------------------------------------------
# The clones com from Chang-o define clones
# based only on heavy or light
# Based on the distance found in 07_run_immcantation using shazam
# Note, these appear to be found on the whole dataset
clone_info <- read.table(file.path(save_dir, "define_clones",
                                   "immcantation_combined_clone-pass.tsv"),
                         sep = "\t", header = TRUE)

# add c_call to the so

# These are found within samples individually
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
  merge(sample_clone_info, by = "sequence_id_sample", all.x = TRUE) %>%
  dplyr::mutate(cell_sample = paste(sample, cell_id, sep = "_")) %>%
  dplyr::select(-sample)

meta_data <- seurat_data[[]] %>% 
  dplyr::select(sample, final_celltype, Status, tet_name_cutoff, v_gene,
                j_gene, chains, dplyr::all_of(c(colnames(assay_data1),
                colnames(assay_data2))),
                isotype, full_tet_name_cutoff, scar_libra_tet_hash_id, 
                scar_libra_full_hash_id, Sample.Name) %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::mutate(barcode = gsub("_[0-9]+", "", barcode)) %>%
  dplyr::mutate(cell_sample = barcode) %>%
  dplyr::select(-Sample.Name)

clone_info <- merge(clone_info, meta_data, by = "cell_sample", all.x = FALSE,
                    all.y = TRUE)


# Check number of non-b cells with v gene call
# We see exactly the number we expect 
table(is.na(seurat_data$v_gene), seurat_data$final_celltype)
clone_info_check <- clone_info %>%
  dplyr::select(cell_sample, final_celltype) %>%
  dplyr::distinct()
table(clone_info_check$final_celltype)

# clone_id, sample_clone_id
# clone_id is called across all samples, sample_clone_id is called only
# within samples

count_clones <- clone_info %>%
  dplyr::group_by(cell_sample) %>%
  dplyr::mutate(sequences = paste(sequence, collapse = ";"),
                   alignment_sequences = paste(sequence_alignment, collapse = ";")) %>%
  dplyr::filter(locus == "IGH") %>%
  dplyr::group_by(clone_id, v_gene, j_gene) %>%
  dplyr::add_count(name = "clone_count") %>%
  dplyr::filter(clone_count >= 3) %>%
  dplyr::select(v_gene, j_gene, cdr3, productive,
                cdr3, clone_id, sample, final_celltype,
                Status, tet_name_cutoff, chains, isotype,
                dplyr::all_of(c(colnames(assay_data1),
                              colnames(assay_data2))),
                clone_count, full_tet_name_cutoff,
                sequences, alignment_sequences, scar_libra_tet_hash_id,
                scar_libra_full_hash_id)


# Add in counts for individuals --> how many individuals are seen?
count_clones <- count_clones %>%
  dplyr::group_by(clone_id, v_gene, j_gene) %>%
  dplyr::mutate(number_of_samples = length(unique(sample)))

# Add in counts for isotype, cell type and binding
count_clones <- count_clones %>% 
  dplyr::group_by(clone_id, v_gene, j_gene, isotype) %>%
  dplyr::add_count(name = "isotype_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(clone_id, v_gene, j_gene, final_celltype) %>%
  dplyr::add_count(name = "cell_type_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(clone_id, v_gene, j_gene, tet_name_cutoff) %>%
  dplyr::add_count(name = "tet_binding_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(clone_id, v_gene, j_gene, Status) %>%
  dplyr::add_count(name = "status_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(clone_id, v_gene, j_gene, sample) %>%
  dplyr::add_count(name = "sample_count") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(isotype_percent = isotype_count / clone_count * 100,
                cell_type_percent = cell_type_count / clone_count * 100,
                tet_binding_percent = tet_binding_count / clone_count * 100,
                status_percent = status_count / clone_count * 100) %>%
  dplyr::group_by(clone_id) %>%
  dplyr::arrange(v_gene, j_gene) %>%
  dplyr::mutate(final_clone = paste(clone_id, cumsum(!duplicated(v_gene, j_gene)),
                                    sep = "_"))

count_clones %>%
  dplyr::filter(clone_id == "19880") %>%
  data.frame()

# Break up into clones that are shared between individuals and those that are now
column_order <- c("clone_id", "final_clone", "v_gene", "j_gene", "cdr3",
                  "chains", "isotype", "productive", "sample", 
                  "final_celltype",  "Status", "tet_name_cutoff",
                  "full_tet_name_cutoff", "scar_libra_tet_hash_id",
                  "scar_libra_full_hash_id", colnames(assay_data1),
                  colnames(assay_data2),"clone_count", "sample_count",
                  "number_of_samples", "isotype_count", "isotype_percent", 
                  "cell_type_count", "cell_type_percent", "tet_binding_count", 
                  "tet_binding_percent", "status_count", "status_percent",
                  "sequences", "alignment_sequences")

shared_clones <- count_clones %>%
  dplyr::filter(number_of_samples > 1) %>%
  dplyr::arrange(desc(clone_count)) %>%
  dplyr::select(dplyr::all_of(column_order))

expanded_clones <- count_clones %>%
  dplyr::filter(sample_count > 1) %>%
  dplyr::arrange(desc(clone_count)) %>%
  dplyr::select(dplyr::all_of(column_order))


new_wb <- openxlsx::createWorkbook() 

openxlsx::addWorksheet(wb = new_wb, sheetName = "public_clones")
openxlsx::writeData(wb = new_wb, sheet = "public_clones", x = shared_clones)

openxlsx::addWorksheet(wb = new_wb, sheetName = "expanded_clones")
openxlsx::writeData(wb = new_wb, sheet = "expanded_clones", x = expanded_clones)

openxlsx::saveWorkbook(new_wb, file.path(v_counting_dir,
                                          "clone_expansion.xlsx"), 
                       overwrite = TRUE)

test <- count_clones %>%
  dplyr::filter(clone_id == "110770") %>%
  data.frame


# Pull out any clones that are 75% or higher for INS, GAD or IA2
high_binding <- count_clones %>%
  dplyr::filter(tet_name_cutoff %in% c("INS-tet", "GAD-tet", "IA2-tet") &
                  tet_binding_percent >= 0.75) %>%
  dplyr::select(v_gene, j_gene, clone_id, tet_name_cutoff, tet_binding_percent) %>%
  dplyr::distinct()


# Add clone info to unmodified seurat object
seurat_data <- readRDS(file.path(save_dir, "rda_obj", 
                                 "seurat_processed_no_doublet.rds"))


add_info <- clone_info %>%
  dplyr::filter(locus == "IGH") %>%
  dplyr::filter(cell_sample %in% colnames(seurat_data)) %>%
  dplyr::group_by(clone_id) %>%
  dplyr::arrange(v_gene, j_gene) %>%
  dplyr::mutate(final_clone = paste(clone_id, cumsum(!duplicated(v_gene, j_gene)),
                                    sep = "_")) %>%
  dplyr::ungroup() %>%
  dplyr::select(final_clone, cell_sample) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("cell_sample")


seurat_data <- AddMetaData(seurat_data, add_info)

# Add in the c gene call
c_gene <- clone_info %>%
  dplyr::filter(locus == "IGH") %>%
  dplyr::select(c_call, cell_sample) %>%
  dplyr::filter(cell_sample %in% colnames(seurat_data)) %>%
  dplyr::distinct() %>%
  dplyr::rename(imcantation_isotype = c_call) %>%
  tibble::column_to_rownames("cell_sample")

seurat_data <- AddMetaData(seurat_data, c_gene)


saveRDS(seurat_data, file.path(save_dir, "rda_obj", 
                                 "seurat_processed_no_doublet.rds"))


# test_data <- seurat_data[[]] %>%
#   dplyr::select(v_gene, j_gene, final_clone) %>%
#   dplyr::mutate(clone_id = gsub("_.*", "", final_clone))
# 
# test_data %>%
#   dplyr::filter(clone_id == 64635) %>%
#   head()
