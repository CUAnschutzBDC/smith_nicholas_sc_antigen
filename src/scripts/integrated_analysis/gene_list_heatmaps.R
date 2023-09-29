library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)
library(djvdj)
library(UpSetR)
library(muscat)
library(SingleCellExperiment)

source(here("src/scripts/muscat_plotting_functions.R"))

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "merged"

# Set directories
base_dir <- here()


save_dir <- file.path(base_dir, "results", "R_analysis", sample)
fig_dir <- file.path(save_dir, "images", "figures")

ifelse(!dir.exists(fig_dir), dir.create(fig_dir), FALSE)

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))
gene_lists <- readRDS(here("results/mia_gene_lists.rds"))
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
seurat_data <- subset(seurat_data, subset = RNA_combined_celltype != "undetermined")
seurat_data$tet_hash_id <- factor(seurat_data$tet_hash_id,
                                  levels = names(tetramer_colors))

celltype_levels <- c("Naive_1", "Naive_3",
                     "BND2", "B.intermediate",
                     "Memory_IgA", "Resting_memory", 
                     "Activated_memory",
                     "Plasmablast",
                     "CD14.Mono", "CD8.TEM")

seurat_data$RNA_combined_celltype <- factor(seurat_data$RNA_combined_celltype,
                                            levels = celltype_levels)


# gene set plots ---------------------------------------------------------------
ifelse(!dir.exists(file.path(fig_dir, "gene_list_heatmaps")),
       dir.create(file.path(fig_dir, "gene_list_heatmaps")), FALSE)
# First remove all non-b cells
celltypes_keep <- c("Naive_1", "Naive_3",
                    "BND2", "B.intermediate",
                    "Memory_IgA", "Resting_memory", 
                    "Activated_memory",
                    "Plasmablast")
seurat_data <- subset(seurat_data, subset = RNA_combined_celltype %in% celltypes_keep)


# Any diabetes antigen and no binding
name_mapping <- c("INS-tet" = "diabetes_antigen",
                  "GAD-tet" = "diabetes_antigen",
                  "IA2-tet" = "diabetes_antigen",
                  "TET-tet" = "Tet_antigen",
                  "Negative" = "Negative",
                  "Doublet" = "other",
                  "DNA-tet" = "other")

seurat_data$test_id <- name_mapping[as.character(seurat_data$tet_hash_id)]

make_de_plots <- function(genes_plot, save_name, vsd, seurat_object,
                          fig_dir){
  
  ifelse(!dir.exists(file.path(fig_dir)),
         dir.create(file.path(fig_dir),
                    recursive = TRUE), FALSE)
  fig_width <- 8
  genes_plot_full <- genes_plot[genes_plot %in% rownames(seurat_object)]
  
  fig_height <- length(genes_plot_full) / 2
  
  # Heatmap of genes
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
  
  seurat_object$meta_col <- paste(seurat_object$sample, seurat_object$test_id,
                                  sep = "_")
  meta_df <- seurat_object[[c("sample", "test_id", "Status", "meta_col")]]
  meta_ave <- meta_df
  rownames(meta_ave) <- NULL
  meta_ave <- distinct(meta_ave)
  rownames(meta_ave) <- meta_ave$meta_col
  
  meta_ave$sample <- factor(meta_ave$sample, levels = sample_order)
  meta_ave$test_id <- factor(meta_ave$test_id, levels = antigen_order)
  
  meta_ave <- meta_ave %>%
    dplyr::arrange(sample, test_id)
  
  color_list <- list("sample" = sample_colors[names(sample_colors) %in% 
                                                unique(meta_df$sample)],
                     "test_id" = antigen_colors,
                     "Status" = status_colors)
  
  
  pdf(file.path(fig_dir,
                paste0(save_name, "_no_nd_all.pdf")),
      width = fig_width, height = fig_height)
  print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave,
                               gene_list = genes_plot_full,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = TRUE))
  dev.off()
  
  # plot only the diabetes antigen
  subset_meta <- meta_ave %>%
    dplyr::filter(test_id == "diabetes_antigen") 
  
  subset_meta$sample <- factor(subset_meta$sample, levels = sample_order)
  
  subset_meta <- subset_meta %>%
    dplyr::arrange(sample)
  
  pdf(file.path(fig_dir, paste0(save_name, "_no_nd_diabetes_antigen.pdf")),
      width = fig_width, height = fig_height)
  
  print(make_corrected_heatamp(vsd = vsd, meta_ave = subset_meta,
                               gene_list = genes_plot_full,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = TRUE))
  
  dev.off()
  
  # plot only the diabetes antigen genes with tet
  subset_meta <- meta_ave %>%
    dplyr::filter(test_id %in% c("Tet_antigen", "diabetes_antigen"))
  
  subset_meta$sample <- factor(subset_meta$sample, levels = sample_order)
  subset_meta$test_id <- factor(subset_meta$test_id, levels = antigen_order)
  
  subset_meta <- subset_meta %>%
    dplyr::arrange(test_id, sample)
  
  subset_meta2 <- subset_meta %>%
    dplyr::arrange(sample, test_id)
  
  pdf(file.path(fig_dir,
                paste0(save_name, "_no_nd_diabetes_antigen_with_tet.pdf")),
      width = fig_width, height = fig_height)
  
  print(make_corrected_heatamp(vsd = vsd, meta_ave = subset_meta,
                               gene_list = genes_plot_full,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = TRUE))
  
  grid::grid.newpage()
  
  print(make_corrected_heatamp(vsd = vsd, meta_ave = subset_meta2,
                               gene_list = genes_plot_full,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = TRUE))
  
  dev.off()
  
}

# De genes made in 07_de_analysis.R
vsd_bcells <- readRDS(file.path(save_dir, "rda_obj", 
                                "psuedobulk_capture_batch_corrected_object.rds"))
all_res_bcells <- read.csv(file.path(save_dir, "files", 
                                     "pseudobulk_all_bcell_de_capture_correct.csv"))

all_plots <- lapply(names(gene_lists), function(x){
  print(x)
  make_de_plots(genes_plot = gene_lists[[x]], 
                save_name = paste(x, "b_cells", sep = "_"), 
                vsd = vsd_bcells, seurat_object = seurat_data,
                fig_dir = file.path(fig_dir, "gene_list_heatmaps", "all_bcells"))  
})



celltype_mapping <- c("Naive_1" = "Naive",
                      "Naive_2" = "Naive",
                      "Naive_3" = "Naive",
                      "Early_memory" = "Memory",
                      "Memory_IgE_IgG" = "Memory",
                      "Resting_memory" = "Memory",
                      "Activated_memory" = "Memory",
                      "Memory_IgA" = "Memory")

vsd_memory <- readRDS(file.path(save_dir, "rda_obj", 
                                "psuedobulk_capture_batch_corrected_Memory_object.rds"))
all_res_memory <- read.csv(file.path(save_dir, "files", 
                                     "pseudobulk_Memory_de_capture_correct.csv"))

seurat_memory <- subset(seurat_data, 
                        subset = RNA_combined_celltype %in%
                          names(celltype_mapping[celltype_mapping == "Memory"]))

all_plots <- lapply(names(gene_lists), function(x){
  print(x)
  genes_to_plot <- gene_lists[[x]]
  if(x == "plasma_cell_differentiation"){
    genes_to_plot <- genes_to_plot[genes_to_plot != "SDC1"]
  }
  make_de_plots(genes_plot = genes_to_plot, 
                save_name = paste(x, "memory", sep = "_"), 
                vsd = vsd_memory, seurat_object = seurat_memory,
                fig_dir = file.path(fig_dir, "gene_list_heatmaps", "memory_bcells"))  
})

vsd_naive <- readRDS(file.path(save_dir, "rda_obj", 
                               "psuedobulk_capture_batch_corrected_Memory_object.rds"))
all_res_naive <- read.csv(file.path(save_dir, "files", 
                                    "pseudobulk_Naive_de_capture_correct.csv"))

seurat_naive <- subset(seurat_data, 
                        subset = RNA_combined_celltype %in%
                          names(celltype_mapping[celltype_mapping == "Naive"]))

all_plots <- lapply(names(gene_lists), function(x){
  make_de_plots(genes_plot = gene_lists[[x]], 
                save_name = paste(x, "naive", sep = "_"), 
                vsd = vsd_naive, seurat_object = seurat_data,
                fig_dir = file.path(fig_dir, "gene_list_heatmaps", "naive_bcells"))  
})
