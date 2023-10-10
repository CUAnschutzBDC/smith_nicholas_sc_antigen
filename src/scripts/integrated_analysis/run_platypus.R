library(Platypus)
library(here)
library(scAnalysisR)
library(tidyverse)
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

vdj_dir <- file.path(save_dir, "images", "vdj_plots")

ifelse(!dir.exists(vdj_dir), dir.create(vdj_dir), FALSE)

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

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

# Analysis ---------------------------------------------------------------------

# The clones com from Chang-o define clones
# based only on heavy or light
# Based on the distance found in 07_run_immcantation using shazam
# Note, these appear to be found within samples individually
# clone_info <- read.table(file.path(save_dir, "define_clones",
#                                    "immcantation_combined_clone-pass.tsv"),
#                          sep = "\t", header = TRUE)
# 
# sample_clone_info <- read.table(file.path(save_dir, "define_clones",
#                                           "immcantation_combined_clone-pass_sample.tsv"),
#                                 sep = "\t", header = TRUE)
# 
# sample_clone_info <- sample_clone_info %>%
#   dplyr::select(sequence_id, clone_id, locus, sample) %>%
#   dplyr::rename(sample_clone_id = clone_id,
#                 sample_locus = locus) %>%
#   dplyr::mutate(sequence_id_sample = paste(sequence_id,
#                                            sample, sep = "_")) %>%
#   dplyr::select(-sequence_id, -sample)
# 
# clone_info <- clone_info %>%
#   dplyr::mutate(sequence_id_sample = paste(sequence_id,
#                                            sample, sep = "_")) %>%
#   merge(sample_clone_info, by = "sequence_id_sample", all.x = TRUE)
# 

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

keep_chains <- c("IGH", "IGH;IGK", "IGH;IGK;IGK",
                 "IGH;IGK;IGL", "IGH;IGL", 
                 "IGH;IGL;IGL")
chain_use <- c("IGH")

all_info_split <- all_info_split %>%
  dplyr::filter(all_chains %in% keep_chains,
                chains %in% chain_use)

all_plots <- VDJ_abundances(VDJ = all_info_split,
                            feature.columns = 'cdr3', 
                            grouping.column = 'tet_hash_id', 
                            sample.column = 'sample', 
                            output.format = 'density.ridges')



names(all_plots) <- unique(all_info_split$sample)

# Try to make dim reduced plots ------------------------------------------------
sample_info <- all_info_split %>%
  dplyr::select(Status, sample) %>%
  dplyr::distinct()


pca_plot <- VDJ_ordination(all_info_split, 
                           feature.columns = 'cdr3', 
                           grouping.column = 'sample',
                           method = 'pca',
                           VDJ.VJ.1chain = FALSE,
                           reduction.level = 'groups')

pca_plot_data <- pca_plot$data %>%
  dplyr::rename(sample = group) %>%
  merge(sample_info)

status_pca <- ggplot2::ggplot(pca_plot_data, ggplot2::aes(x = DIM1,
                                                           y = DIM2,
                                                           color = Status)) +
  ggplot2::geom_point(size = 5) +
  ggplot2::scale_color_manual(values = status_colors)

# install.packages("umap")
umap_plot <- VDJ_ordination(all_info_split, 
                            feature.columns = 'cdr3', 
                            grouping.column = 'sample',
                            method = 'umap',
                            VDJ.VJ.1chain = FALSE,
                            reduction.level = 'groups')


umap_plot_data <- umap_plot$data %>%
  dplyr::rename(sample = group) %>%
  merge(sample_info)

status_umap <- ggplot2::ggplot(umap_plot_data, ggplot2::aes(x = DIM1,
                                                            y = DIM2,
                                                            color = Status)) +
  ggplot2::geom_point(size = 5) +
  ggplot2::scale_color_manual(values = status_colors)

umap_vj_plot <- VDJ_ordination(all_info_split, 
                            feature.columns = c('v_gene', "j_gene"), 
                            grouping.column = 'sample',
                            method = 'umap',
                            VDJ.VJ.1chain = FALSE,
                            reduction.level = 'groups')


umap_plot_data <- umap_vj_plot$data %>%
  dplyr::rename(sample = group) %>%
  merge(sample_info)

status_umap <- ggplot2::ggplot(umap_plot_data, ggplot2::aes(x = DIM1,
                                                            y = DIM2,
                                                            color = Status)) +
  ggplot2::geom_point(size = 5) +
  ggplot2::scale_color_manual(values = status_colors)



pca_vj_plot <- VDJ_ordination(all_info_split, 
                               feature.columns = c('v_gene', "j_gene"), 
                               grouping.column = 'sample',
                               method = 'pca',
                               VDJ.VJ.1chain = FALSE,
                               reduction.level = 'groups')


pca_plot_data <- pca_vj_plot$data %>%
  dplyr::rename(sample = group) %>%
  merge(sample_info)

status_pca <- ggplot2::ggplot(pca_plot_data, ggplot2::aes(x = DIM1,
                                                            y = DIM2,
                                                            color = Status)) +
  ggplot2::geom_point(size = 5) +
  ggplot2::scale_color_manual(values = status_colors)

# Figure out how to run this PCA on your own so you can pull out genes

#Chao1 index for CDRH3s
diversity_plot <- VDJ_diversity(all_info_split, 
                                feature.columns = 'cdr3', 
                                grouping.column = 'sample', 
                                metric = 'shannon',
                                VDJ.VJ.1chain = FALSE,
                                subsample.to.same.n = T)
