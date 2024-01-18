library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)
library(muscat)
library(SummarizedExperiment)

source(here("src/scripts/muscat_plotting_functions.R"))

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
fig_dir <- file.path(save_dir, "images", "figures")

ifelse(!dir.exists(fig_dir), dir.create(fig_dir), FALSE)

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))

vsd <- readRDS(file.path(save_dir, "rda_obj", 
                         "psuedobulk_capture_batch_corrected_object.rds"))

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


# Any diabetes antigen and no binding
name_mapping <- c("INS-tet" = "diabetes_antigen",
                  "GAD-tet" = "diabetes_antigen",
                  "IA2-tet" = "diabetes_antigen",
                  "TET-tet" = "Tet_antigen",
                  "Negative" = "Negative",
                  "Doublet" = "other",
                  "DNA-tet" = "other")

seurat_data$test_id <- name_mapping[as.character(seurat_data$tet_hash_id)]

# Make heatmaps ----------------------------------------------------------------
# Set factor for figuring out size of pdf is labeling
divide_by <- 5

# Put your gene list here, either manual or read from a file
# gene_list <- c("gene1", "gene2", "gene3")
# gene_list <- read.csv("/path/to/list")$gene
gene_list <- c("CD79A", "CD79B", "LYN", "SYK", "BTK", "BLNK")

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


# Make a meta dataframe for coloring. One column will need to perfectly 
# match the colnames in the heatmap. As long as you want the sample,
# status, and antigen on the top as colors, you won't need to change this
seurat_data$meta_col <- paste(seurat_data$sample, seurat_data$test_id,
                                sep = "_")
meta_df <- seurat_data[[c("sample", "test_id", "Status", "meta_col")]]
meta_ave <- meta_df
rownames(meta_ave) <- NULL
meta_ave <- distinct(meta_ave)
rownames(meta_ave) <- meta_ave$meta_col

meta_ave$sample <- factor(meta_ave$sample, levels = sample_order)
meta_ave$test_id <- factor(meta_ave$test_id, levels = antigen_order)

# IF YOU DON'T WANT TO PLOT ALL SAMPLES, SUBSET THE META DF HERE:
# Uncomment this and comment the next two lines.
# meta_ave <- meta_ave %>%
#   dplyr::arrange(sample, test_id) %>%
#   dplyr::filter(Status %in% c("no", "nd"))

meta_ave <- meta_ave %>%
  dplyr::arrange(sample, test_id) 

color_list <- list("sample" = sample_colors[names(sample_colors) %in% 
                                              unique(meta_df$sample)],
                   "test_id" = antigen_colors,
                   "Status" = status_colors)

# This is the heatmap with your genes and any samples that were in the
# seurat object
pdf(file.path(fig_dir, paste0(save_name, "_no_nd_all.pdf")))
print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave,
                             gene_list = gene_list,
                             meta_col = "meta_col", plot_meta_col = FALSE,
                             max_val = 2.5, min_val = -2.5,
                             cluster_rows = TRUE, cluster_cols = FALSE,
                             coloring = color_list, plot_rownames = FALSE))
dev.off()

fig_height <- round(length(gene_list) / divide_by)

# This is the same heatmap as above but with gene names
pdf(file.path(fig_dir, paste0(save_name, "_no_nd_all_name.pdf")),
    width = 8, height = fig_height)
print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave,
                             gene_list = gene_list,
                             meta_col = "meta_col", plot_meta_col = FALSE,
                             max_val = 2.5, min_val = -2.5,
                             cluster_rows = TRUE, cluster_cols = FALSE,
                             coloring = color_list, plot_rownames = TRUE))
dev.off()

# This is how to make the heatmap with just one of the antigens
meta_ave_2 <- meta_ave %>%
  dplyr::filter(test_id == "diabetes_antigen")

pdf(file.path(fig_dir, paste0(save_name, "_no_nd_all_genes_diabetes_plot.pdf")))
print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave_2,
                             gene_list = gene_list,
                             meta_col = "meta_col", plot_meta_col = FALSE,
                             max_val = 2.5, min_val = -2.5,
                             cluster_rows = TRUE, cluster_cols = FALSE,
                             coloring = color_list, plot_rownames = FALSE))
dev.off()

# Same as above but with gene names labeled
fig_height <- round(length(gene_list) / divide_by)
pdf(file.path(fig_dir, paste0(save_name, "_no_nd_all_genes_diabetes_plot_name.pdf")),
    width = 8, height = fig_height)
print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave_2,
                             gene_list = gene_list,
                             meta_col = "meta_col", plot_meta_col = FALSE,
                             max_val = 2.5, min_val = -2.5,
                             cluster_rows = TRUE, cluster_cols = FALSE,
                             coloring = color_list, plot_rownames = TRUE))
dev.off()
