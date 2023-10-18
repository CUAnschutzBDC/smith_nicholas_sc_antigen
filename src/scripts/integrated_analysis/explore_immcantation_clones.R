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
# Note, these appear to be found on the whole dataset
clone_info <- read.table(file.path(save_dir, "define_clones",
                                   "immcantation_combined_clone-pass.tsv"),
                         sep = "\t", header = TRUE)

# These are found within samples individually
sample_clone_info <- read.table(file.path(save_dir, "define_clones",
                                   "immcantation_combined_clone-pass_sample.tsv"),
                         sep = "\t", header = TRUE)

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


# Check number of non-b cells with v gene call
# We see exactly the number we expect 
table(is.na(seurat_data$v_gene), seurat_data$RNA_combined_celltype)
clone_info_check <- clone_info %>%
  dplyr::select(cell_sample, RNA_combined_celltype) %>%
  dplyr::distinct()
table(clone_info_check$RNA_combined_celltype)

# Subset to only B cells
b_cells <- c("Activated_memory", "B.intermediate", "BND2", "Memory_IgA",
             "Naive_1", "Naive_3", "Plasmablast", "Resting_memory")

clone_info <- clone_info %>%
  dplyr::filter(RNA_combined_celltype %in% b_cells)

# This just makes sure that heavy and light are always in different clones
# They are so it doesn't need to be run
# heavy_light <- clone_info %>%
#   dplyr::select(clone_id, locus) %>%
#   dplyr::ungroup() %>%
#   dplyr::distinct() %>%
#   dplyr::group_by(clone_id) %>%
#   dplyr::add_count(name = "locus_count") 

ifelse(!dir.exists(file.path(save_dir, "images", "immcantation_clones")),
       dir.create(file.path(save_dir, "images", "immcantation_clones")),
       FALSE)

make_plots <- function(clone_info_expanded, name, percent_keep = 25,
                       locus_use = "IGH", clone_min = 3){
  expanded <- clone_info_expanded %>%
    dplyr::filter(locus == locus_use, clone_count > clone_min) 
  
  # Number of clones across samples
  tetramer_group_across <- expanded %>%
    dplyr::group_by(clone_id, tet_hash_id) %>%
    dplyr::add_count(name = "tetramer_count") %>%
    dplyr::mutate(tetramer_fraction = tetramer_count / clone_count * 100) %>%
    dplyr::mutate(clone_name = paste0("clone_", clone_id)) %>%
    dplyr::select(clone_id, clone_name, tet_hash_id, tetramer_count,
                  tetramer_fraction, clone_count) %>%
    dplyr::distinct()

  plot_one <- single_plot(tetramer_group = tetramer_group_across,
                          name = name, type = "across",
                          percent_keep = percent_keep)
  
  # Number of clones within samples
  expanded <- clone_info_expanded %>%
    dplyr::filter(locus == locus_use, sample_clone_count > clone_min) 
  
  tetramer_group_within <- expanded %>%
    dplyr::group_by(sample_clone_id, tet_hash_id) %>%
    dplyr::add_count(name = "tetramer_count") %>%
    dplyr::mutate(tetramer_fraction = tetramer_count / sample_clone_count * 100) %>%
    dplyr::mutate(sample_clone_name = paste0("clone_", sample_clone_id)) %>%
    dplyr::select(sample_clone_id, sample_clone_name, tet_hash_id, tetramer_count,
                  tetramer_fraction, sample_clone_count, sample) %>%
    dplyr::distinct()
  
  colnames(tetramer_group_within) <- gsub("sample_", "", 
                                          colnames(tetramer_group_within))
  
  plot_two <- single_plot(tetramer_group = tetramer_group_within,
                          name = name, type = "within",
                          percent_keep = percent_keep)
  
  # Clones shared between samples
  expanded <- clone_info_expanded %>%
    dplyr::filter(locus == locus_use, clone_count > clone_min) 
  
  tetramer_group_sample <- expanded %>%
    dplyr::group_by(clone_id, tet_hash_id) %>%
    dplyr::add_count(name = "tetramer_count") %>%
    dplyr::mutate(tetramer_fraction = tetramer_count / clone_count * 100) %>%
    dplyr::mutate(clone_name = paste0("clone_", clone_id)) %>%
    dplyr::select(clone_id, clone_name, tet_hash_id, tetramer_count,
                  tetramer_fraction, clone_count, sample)
  
  keep_clones <- tetramer_group_sample %>%
    dplyr::ungroup() %>%
    dplyr::select(clone_name, sample) %>%
    dplyr::distinct() %>%
    dplyr::group_by(clone_name) %>%
    dplyr::add_count(name = "number_samples") %>%
    dplyr::filter(number_samples > 1)
  
  tetramer_group_sample <- tetramer_group_sample %>%
    dplyr::filter(clone_name %in% unique(keep_clones$clone_name)) %>%
    dplyr::select(clone_id, clone_name, tet_hash_id, tetramer_count,
                  tetramer_fraction, clone_count) %>%
    dplyr::distinct()
  
  plot_three <- single_plot(tetramer_group = tetramer_group_sample,
                            name = name, type = "between",
                            percent_keep = percent_keep)
  
  return(list("across" = plot_one,
              "within" = plot_two,
              "between" = plot_three))
  
 
}

single_plot <- function(tetramer_group, name,
                        type, percent_keep = 25){
  
  if(type == "across"){
    save_name <- paste0("All clones across ", name, " samples")
  } else if(type == "within"){
    save_name <- paste0("All clones within ", name, " samples")
  } else if(type == "between"){
    save_name <- paste0("All clones shared between ", name, " samples")
  } else {
    stop("Type can only be across, within, or between")
  }
  
  interesting_clones <- tetramer_group %>%
    dplyr::filter(tet_hash_id %in% c("GAD-tet", "IA2-tet",
                                     "INS-tet")) %>%
    dplyr::group_by(clone_name) %>%
    dplyr::mutate(percent_diabetes = sum(tetramer_fraction)) %>%
    dplyr::filter(percent_diabetes > percent_keep)
  
  plot_data <- tetramer_group %>%
    dplyr::filter(clone_name %in% unique(interesting_clones$clone_name))
  
  count_data <- plot_data %>%
    dplyr::select(clone_name, clone_count) %>%
    dplyr::mutate(clone_count = paste0(clone_count, " cells")) %>%
    dplyr::distinct()
  
  return_plot <-ggplot2::ggplot(plot_data, ggplot2::aes(x = clone_name,
                                                        y = tetramer_fraction,
                                                        fill = tet_hash_id)) +
    ggplot2::geom_bar(position = "stack", stat = "identity") +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggplot2::scale_fill_manual(values = tetramer_colors) +
    ggplot2::geom_text(
      data = count_data, 
      ggplot2::aes(clone_name, 50, label = clone_count),
      show.legend = F,
      inherit.aes = F,
      color = "white", 
      angle = 90
    ) +
    ggplot2::ggtitle(save_name)
  
  if(nrow(plot_data > 0)){
    return(return_plot)
  } else {
    return(NULL)
  }
}

pdf(file.path(save_dir, "images", "immcantation_clones", "expansion_plots.pdf"),
    width = 12, height = 8)

# Find expanded clones all samples
# Find expanded clones
clone_info_expanded <- clone_info %>%
  dplyr::group_by(clone_id) %>%
  dplyr::add_count(name = "clone_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sample_clone_id) %>%
  dplyr::add_count(name = "sample_clone_count") %>%
  dplyr::ungroup()

return_plots <- make_plots(clone_info_expanded = clone_info_expanded, 
                           name = "all",
                           percent_keep = 25)


print(return_plots)
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

print(return_plots)

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

print(return_plots)

# Expanded clones only no 
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

print(return_plots)

# Expanded clones only nd
clone_info_expanded <- clone_info %>%
  dplyr::filter(Status == "nd") %>%
  dplyr::group_by(clone_id) %>%
  dplyr::add_count(name = "clone_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sample_clone_id) %>%
  dplyr::add_count(name = "sample_clone_count") %>%
  dplyr::ungroup()

return_plots <- make_plots(clone_info_expanded = clone_info_expanded,
                           name = "non diabetic",
                           percent_keep = 25)

print(return_plots)

dev.off()