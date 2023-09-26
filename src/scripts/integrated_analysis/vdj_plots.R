library(here)
library(scAnalysisR)
library(pheatmap)
library(tidyverse)
library(splitstackshape)
library(circlize)
library(viridis)

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



# Circos plots to make
# 1. Heavy V+J pairs
# 2. Light V+J pairs
# 3. Heavy + Light V+J pairs
# 4. Heavy VJ pairs between samples
# 5. Light VJ pairs between samples
# 6. Heavy + Light VJ pairs between samples
# For each, also separate by chain...

# Functions

make_circos_plot <- function(save_name, circos_df, color = NULL){
  pdf(save_name, height = 17, width = 17)
  
  chordDiagram(circos_df, annotationTrack = "grid", 
               preAllocateTracks = 1, col = color)
  
  # we go back to the first track and customize sector labels
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(-0.1, 0.5))
  }, bg.border = NA) # here set bg.border to NA is important
  
  dev.off()
  
  circos.clear()
}

count_genes <- function(starting_df, group_by = "all",
                        chain = "IGH",
                        keep_chains = c("IGH", "IGH;IGK",
                                        "IGH;IGK;IGK", "IGH;IGK;IGL",
                                        "IGH;IGL", "IGH;IGL;IGL"),
                        color_list = NULL){
  if(group_by == "all"){
    select_cols <- c("v_gene", "j_gene")
  } else if (group_by == "sample") {
    select_cols <- c("v_gene", "j_gene", "sample")
  } else if (group_by == "Status") {
    select_cols <- c("v_gene", "j_gene", "Status")
  } else {
    stop("group_by must be 'all', 'sample', or 'Status'")
  }
  
  return_df <- starting_df %>%
    dplyr::filter(all_chains %in% keep_chains) %>%
    dplyr::filter(chains == chain) %>%
    dplyr::select(dplyr::all_of(select_cols))
  
  return_df <- return_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(select_cols))) %>% 
    dplyr::add_count(name = "value") %>%
    dplyr::distinct() %>%
    dplyr::ungroup()
  
  if(group_by != "all"){
    return_df$color <- color_list[return_df[[group_by]]]
    color_values <- return_df$color
    return_df <- return_df %>%
      dplyr::select(-dplyr::all_of(group_by), -color)
  } else{
    color_values <- NULL
  }
  
  return(list("df" = return_df, "color_list" = color_values))
}

# Heavy V+J pairs --------------------------------------------------------------
# I'll start with 10x
all_data <- count_genes(starting_df = all_info_split,
                        group_by = "all")
  
make_circos_plot(save_name = file.path(vdj_dir,
                                       "heavy_v_j_circos.pdf"),
                 circos_df = all_data$df,
                 color = all_data$color_list)


## Repeat by sample ------------------------------------------------------------
sample_data <- count_genes(starting_df = all_info_split,
                           group_by = "sample", color_list = sample_colors)

make_circos_plot(save_name = file.path(vdj_dir,
                                       "heavy_v_j_circos_sample_color.pdf"),
                 circos_df = sample_data$df,
                 color = sample_data$color_list)

## Repeat by status ------------------------------------------------------------
status_data <- count_genes(starting_df = all_info_split,
                           group_by = "Status", color_list = status_colors)

make_circos_plot(save_name = file.path(vdj_dir,
                                       "heavy_v_j_circos_status_color.pdf"),
                 circos_df = status_data$df,
                 color = status_data$color_list)

# Separate by Status -----------------------------------------------------------
invisible(lapply(unique(all_info_split$Status), function(x){
  test_data <- all_info_split %>%
    dplyr::filter(Status == x)
  
  ## All samples
  all_data <- count_genes(starting_df = test_data,
                          group_by = "all")
  
  make_circos_plot(save_name = file.path(vdj_dir,
                                         paste0(x, "_heavy_v_j_circos.pdf")),
                   circos_df = all_data$df,
                   color = all_data$color_list)
  
  
  # Repeat by sample
  sample_data <- count_genes(starting_df = test_data,
                             group_by = "sample", color_list = sample_colors)
  
  make_circos_plot(save_name = file.path(vdj_dir,
                                         paste0(x, "_heavy_v_j_circos_sample_color.pdf")),
                   circos_df = sample_data$df,
                   color = sample_data$color_list)
  
  # Repeat by status
  sample_data <- count_genes(starting_df = test_data,
                             group_by = "Status", color_list = status_colors)
  
  make_circos_plot(save_name = file.path(vdj_dir,
                                         paste0(x, "_heavy_v_j_circos_status_color.pdf")),
                   circos_df = sample_data$df,
                   color = sample_data$color_list)
}))



# Isotype switched status ------------------------------------------------------
isotypes_use <- c("IGHA", "IGHD", "IGHG", "IGHM")

invisible(lapply(isotypes_use, function(x){
  test_data <- all_info_split %>%
    dplyr::filter(isotype == x)
  
  ## All samples
  all_data <- count_genes(starting_df = test_data,
                          group_by = "all")
  
  make_circos_plot(save_name = file.path(vdj_dir,
                                         paste0(x, "_heavy_v_j_circos.pdf")),
                   circos_df = all_data$df,
                   color = all_data$color_list)
  
  
  # Repeat by sample
  sample_data <- count_genes(starting_df = test_data,
                             group_by = "sample", color_list = sample_colors)
  
  make_circos_plot(save_name = file.path(vdj_dir,
                                         paste0(x, "_heavy_v_j_circos_sample_color.pdf")),
                   circos_df = sample_data$df,
                   color = sample_data$color_list)
  
  # Repeat by status
  sample_data <- count_genes(starting_df = test_data,
                             group_by = "Status", color_list = status_colors)
  
  make_circos_plot(save_name = file.path(vdj_dir,
                                         paste0(x, "_heavy_v_j_circos_status_color.pdf")),
                   circos_df = sample_data$df,
                   color = sample_data$color_list)
}))

# Tetramer ---------------------------------------------------------------------
tetramers_use <- c("DNA-tet", "Doublet", "GAD-tet", "IA2-tet",
                   "INS-tet", "Negative", "TET-tet")
invisible(lapply(isotypes_use, function(x){
  test_data <- all_info_split %>%
    dplyr::filter(tet_hash_id == x)
  
  ## All samples
  all_data <- count_genes(starting_df = test_data,
                          group_by = "all")
  
  make_circos_plot(save_name = file.path(vdj_dir,
                                         paste0(x, "_heavy_v_j_circos.pdf")),
                   circos_df = all_data$df,
                   color = all_data$color_list)
  
  
  # Repeat by sample
  sample_data <- count_genes(starting_df = test_data,
                             group_by = "sample", color_list = sample_colors)
  
  make_circos_plot(save_name = file.path(vdj_dir,
                                         paste0(x, "_heavy_v_j_circos_sample_color.pdf")),
                   circos_df = sample_data$df,
                   color = sample_data$color_list)
  
  # Repeat by status
  sample_data <- count_genes(starting_df = test_data,
                             group_by = "Status", color_list = status_colors)
  
  make_circos_plot(save_name = file.path(vdj_dir,
                                         paste0(x, "_heavy_v_j_circos_status_color.pdf")),
                   circos_df = sample_data$df,
                   color = sample_data$color_list)
}))


# Heatmap ----------------------------------------------------------------------
all_info_split_heatmap <- all_info_split %>%
  dplyr::filter(all_chains %in% c("IGH", "IGH;IGK", "IGH;IGK;IGK",
                                  "IGH;IGK;IGL", "IGH;IGL", "IGH;IGL;IGL")) %>%
  dplyr::filter(chains == "IGH") %>%
  dplyr::select(v_gene, j_gene)

all_info_split_heatmap <- all_info_split_heatmap %>%
  dplyr::group_by(v_gene, j_gene) %>% 
  dplyr::add_count(name = "value") %>%
  dplyr::distinct()

new_info <- all_info_split_heatmap %>%
  tidyr::pivot_wider(names_from = j_gene, values_from = value) %>%
  tibble::column_to_rownames("v_gene")

new_info[is.na(new_info)] <- 0

pdf(file.path(vdj_dir, "heavy_v_j_heatmap.pdf"))
pheatmap(new_info, color = magma(n = 10))

grid::grid.newpage()
new_info <- t(scale(t(new_info)))

pheatmap(new_info, color = magma(n = 10))

dev.off()
## Heatmap by sample -----------------------------------------------------------
all_info_split_heatmap <- all_info_split %>%
  dplyr::filter(all_chains %in% c("IGH", "IGH;IGK", "IGH;IGK;IGK",
                                  "IGH;IGK;IGL", "IGH;IGL", "IGH;IGL;IGL")) %>%
  dplyr::filter(chains == "IGH") %>%
  dplyr::select(v_gene, j_gene, sample, Status)

all_info_split_heatmap <- all_info_split_heatmap %>%
  dplyr::group_by(v_gene, j_gene, sample) %>% 
  dplyr::add_count(name = "value") %>%
  dplyr::distinct() %>%
  dplyr::mutate(j_sample = paste(j_gene, sample, sep = "_")) %>%
  dplyr::ungroup()

new_info <- all_info_split_heatmap %>%
  dplyr::select(v_gene, j_sample, value) %>%
  tidyr::pivot_wider(names_from = j_sample, values_from = value) %>%
  tibble::column_to_rownames("v_gene")

new_info[is.na(new_info)] <- 0

sample_info <- all_info_split_heatmap %>%
  dplyr::select(j_sample, j_gene, sample, Status) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("j_sample")

j_colors <- RColorBrewer::brewer.pal(n = length(unique(sample_info$j_gene)),
                                     name = "Set1")
names(j_colors) <- unique(sample_info$j_gene)

coloring <- list(sample = sample_colors,
                 j_gene = j_colors,
                 Status = status_colors)

if(!identical(rownames(sample_info), colnames(new_info))){
  new_info <- new_info[ , rownames(sample_info)]
}


pdf(file.path(vdj_dir, "heavy_v_j_heatmap_sample.pdf"))
pheatmap(new_info, color = magma(n = 10), annotation_col = sample_info,
         annotation_colors = coloring)

grid::grid.newpage()
new_info <- t(scale(t(new_info)))

pheatmap(new_info, color = magma(n = 10), annotation_col = sample_info,
         annotation_colors = coloring)

dev.off()

## Heatmap by status -----------------------------------------------------------
all_info_split_heatmap <- all_info_split %>%
  dplyr::filter(all_chains %in% c("IGH", "IGH;IGK", "IGH;IGK;IGK",
                                  "IGH;IGK;IGL", "IGH;IGL", "IGH;IGL;IGL")) %>%
  dplyr::filter(chains == "IGH") %>%
  dplyr::select(v_gene, j_gene, Status)

all_info_split_heatmap <- all_info_split_heatmap %>%
  dplyr::group_by(v_gene, j_gene, Status) %>% 
  dplyr::add_count(name = "value") %>%
  dplyr::distinct() %>%
  dplyr::mutate(j_status = paste(j_gene, Status, sep = "_")) %>%
  dplyr::ungroup()

new_info <- all_info_split_heatmap %>%
  dplyr::select(v_gene, j_status, value) %>%
  tidyr::pivot_wider(names_from = j_status, values_from = value) %>%
  tibble::column_to_rownames("v_gene")

new_info[is.na(new_info)] <- 0

sample_info <- all_info_split_heatmap %>%
  dplyr::select(j_status, j_gene, Status) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("j_status")

j_colors <- RColorBrewer::brewer.pal(n = length(unique(sample_info$j_gene)),
                                     name = "Set1")
names(j_colors) <- unique(sample_info$j_gene)

coloring <- list(j_gene = j_colors,
                 Status = status_colors)

if(!identical(rownames(sample_info), colnames(new_info))){
  new_info <- new_info[ , rownames(sample_info)]
}


pdf(file.path(vdj_dir, "heavy_v_j_heatmap_status.pdf"))
pheatmap(new_info, color = magma(n = 10), annotation_col = sample_info,
         annotation_colors = coloring)

grid::grid.newpage()
new_info <- t(scale(t(new_info)))

pheatmap(new_info, color = magma(n = 10), annotation_col = sample_info,
         annotation_colors = coloring)

dev.off()
