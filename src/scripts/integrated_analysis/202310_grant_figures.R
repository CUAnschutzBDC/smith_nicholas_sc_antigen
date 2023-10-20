# Plots for Mia's grant
# updated UMAP of all B cell subsets identified (maybe removing CD14,
# CD8, and any other non-B cells)
# Barplot of % of each B cell subset found within ND, Stage 1, Stage 2, and 
# T1D. (x-axis is ND, Stage 1, Stage 2, etc) and y axis is percent of
# total B cells 
# Barplot of % of antigen bound cells found in ND, Stage 1, Stage 2, and 
# T1D—remove doublets but can keep negative bound cells (x-axis is ND, Stage 1,
# Stage 2, etc) and y axis is percent of total B cells 
# Barplot of clonally expanded antigen-binding B cells in Stage 1
# Barplot of clonally expanded antigen-binding B cells in Stage 2
# Barplot of clonally expanded antigen-binding B cells in T1D


library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)

source(here("src/scripts/muscat_plotting_functions.R"))
source(here("src/scripts/immcantation_plotting_functions.R"))

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "merged"

all_ref_dir <-
  "/beevol/home/wellskri/Analysis/references/single_cell_references"

HTO <- FALSE
ADT <- FALSE

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
base_dir <- here()

save_dir <- file.path(base_dir, "results", "R_analysis", sample)
fig_dir <- file.path(save_dir, "images", "202310_grant_figures")

ifelse(!dir.exists(fig_dir), dir.create(fig_dir), FALSE)

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_reprocessed.rds"))

plot_tetramers <- c("INS-tet", "TET-tet", "IA2-tet", "GAD-tet", "DNA-tet")

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
# Set levels -------------------------------------------------------------------
seurat_data$tet_hash_id <- factor(seurat_data$tet_hash_id,
                                  levels = names(tetramer_colors))

celltype_levels <- c("Naive_1", "Naive_3",
                     "BND2", "B.intermediate",
                     "Memory_IgA", "Resting_memory", 
                     "Activated_memory",
                     "Plasmablast")

seurat_data$RNA_combined_celltype <- factor(seurat_data$RNA_combined_celltype,
                                            levels = celltype_levels)

# Make seurat subsets ----------------------------------------------------------
seurat_b_no_dub <- subset(seurat_data, subset = tet_hash_id != "Doublet")

# Celltype plots ---------------------------------------------------------------
# updated UMAP of all B cell subsets identified (maybe removing CD14,
# CD8, and any other non-B cells)
# UMAP of celltypes
celltype_plot <- plotDimRed(seurat_data, col_by = "RNA_combined_celltype",
                            color = final_colors, plot_type = "rna_mnn.umap",
                            ggrastr = TRUE, size = 0.25)[[1]]

pdf(file.path(fig_dir, "celltype_umap.pdf"), height = 6, width = 8)

print(celltype_plot)

dev.off()

# Barplot of % of each B cell subset found within ND, Stage 1, Stage 2, and 
# T1D. (x-axis is ND, Stage 1, Stage 2, etc) and y axis is percent of
# total B cells 
sample_barplot <- scAnalysisR::stacked_barplots(seurat_data,
                                                meta_col = "RNA_combined_celltype",
                                                split_by = "Status", 
                                                color = final_colors) +
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(fig_dir, "celltype_barplot.pdf"), height = 8, width = 6)

print(sample_barplot)

dev.off()


# Antigen plots ----------------------------------------------------------------
# Barplot of % of antigen bound cells found in ND, Stage 1, Stage 2, and 
# T1D—remove doublets but can keep negative bound cells (x-axis is ND, Stage 1,
# Stage 2, etc) and y axis is percent of total B cells 
diabetes_barplot <- scAnalysisR::stacked_barplots(seurat_data,
                                                  meta_col = "tet_hash_id",
                                                  split_by = "Status", 
                                                  color = tetramer_colors) +
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(fig_dir, "tetramer_barplot_with_doublet.pdf"), height = 8, width = 6)

print(diabetes_barplot)

dev.off()

tetramer_no_dub_colors <- tetramer_colors[names(tetramer_colors) != "Doublet"]

diabetes_barplot <- scAnalysisR::stacked_barplots(seurat_b_no_dub,
                                                  meta_col = "tet_hash_id",
                                                  split_by = "Status", 
                                                  color = tetramer_no_dub_colors) +
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(fig_dir, "tetramer_barplot.pdf"), height = 8, width = 6)

print(diabetes_barplot)

dev.off()


# Clonal expansion -------------------------------------------------------------
clone_info <- read.table(file.path(save_dir, "define_clones",
                                   "immcantation_combined_clone-pass.tsv"),
                         sep = "\t", header = TRUE)

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
  dplyr::mutate(cell_sample = paste(cell_id, sample, sep = "_"))

meta_data <- seurat_data[[]] %>% 
  dplyr::select(sample, RNA_combined_celltype, Status, tet_hash_id, v_gene) %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::mutate(barcode = gsub("_[0-9]+", "", barcode)) %>%
  dplyr::mutate(cell_sample = paste(barcode, sample, sep = "_")) %>%
  dplyr::select(-sample)

clone_info <- merge(clone_info, meta_data, by = "cell_sample", all.x = TRUE,
                    all.y = FALSE)

## Barplot of clonally expanded antigen-binding B cells in Stage 1--------------
# Expanded clones only aab stage 1
clone_info_expanded <- clone_info %>%
  dplyr::filter(Status == "aab_stage_1") %>%
  dplyr::group_by(clone_id) %>%
  dplyr::add_count(name = "clone_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sample_clone_id) %>%
  dplyr::add_count(name = "sample_clone_count") %>%
  dplyr::ungroup()

return_plots <- make_plots(clone_info_expanded = clone_info_expanded,
                           name = "autoantibody stage 1",
                           percent_keep = 25)
pdf(file.path(fig_dir, "aab_stage_one_clones_with_doublet.pdf"),
    width = 10, height = 8)

print(return_plots)

dev.off()

# Without doublets
clone_info_expanded <- clone_info %>%
  dplyr::filter(Status == "aab_stage_1", tet_hash_id != "Doublet") %>%
  dplyr::group_by(clone_id) %>%
  dplyr::add_count(name = "clone_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sample_clone_id) %>%
  dplyr::add_count(name = "sample_clone_count") %>%
  dplyr::ungroup()

return_plots <- make_plots(clone_info_expanded = clone_info_expanded,
                           name = "autoantibody stage 1",
                           percent_keep = 25)
pdf(file.path(fig_dir, "aab_stage_one_clones.pdf"),
    width = 10, height = 8)

print(return_plots)

dev.off()


## Barplot of clonally expanded antigen-binding B cells in Stage 2 -------------
# Expanded clones only aab stage 2
clone_info_expanded <- clone_info %>%
  dplyr::filter(Status == "aab_stage_2") %>%
  dplyr::group_by(clone_id) %>%
  dplyr::add_count(name = "clone_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sample_clone_id) %>%
  dplyr::add_count(name = "sample_clone_count") %>%
  dplyr::ungroup()

return_plots <- make_plots(clone_info_expanded = clone_info_expanded,
                           name = "autoantibody stage 2",
                           percent_keep = 25)
pdf(file.path(fig_dir, "aab_stage_two_clones_with_doublet.pdf"),
    width = 10, height = 8)

print(return_plots)

dev.off()

# Without doublets
clone_info_expanded <- clone_info %>%
  dplyr::filter(Status == "aab_stage_2", tet_hash_id != "Doublet") %>%
  dplyr::group_by(clone_id) %>%
  dplyr::add_count(name = "clone_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sample_clone_id) %>%
  dplyr::add_count(name = "sample_clone_count") %>%
  dplyr::ungroup()

return_plots <- make_plots(clone_info_expanded = clone_info_expanded,
                           name = "autoantibody stage 2",
                           percent_keep = 25)
pdf(file.path(fig_dir, "aab_stage_two_clones.pdf"),
    width = 10, height = 8)

print(return_plots)

dev.off()

## Barplot of clonally expanded antigen-binding B cells in T1D -----------------
# Expanded clones only aab stage 1
clone_info_expanded <- clone_info %>%
  dplyr::filter(Status == "no") %>%
  dplyr::group_by(clone_id) %>%
  dplyr::add_count(name = "clone_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sample_clone_id) %>%
  dplyr::add_count(name = "sample_clone_count") %>%
  dplyr::ungroup()

return_plots <- make_plots(clone_info_expanded = clone_info_expanded,
                           name = "new onset",
                           percent_keep = 25)
pdf(file.path(fig_dir, "no_clones_with_doublet.pdf"),
    width = 10, height = 8)

print(return_plots)

dev.off()

# Without doublets
clone_info_expanded <- clone_info %>%
  dplyr::filter(Status == "no", tet_hash_id != "Doublet") %>%
  dplyr::group_by(clone_id) %>%
  dplyr::add_count(name = "clone_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sample_clone_id) %>%
  dplyr::add_count(name = "sample_clone_count") %>%
  dplyr::ungroup()

return_plots <- make_plots(clone_info_expanded = clone_info_expanded,
                           name = "new onset",
                           percent_keep = 25)
pdf(file.path(fig_dir, "no_clones.pdf"),
    width = 10, height = 8)

print(return_plots)

dev.off()
