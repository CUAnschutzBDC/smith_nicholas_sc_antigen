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

make_de_plots <- function(all_res, save_name, vsd, seurat_object, fig_dir){
  divide_by <- 5
  genes_plot <- lapply(unique(all_res$contrast), function(x){
    return_res <- all_res %>%
      dplyr::filter(contrast == x & abs(logFC) > 0.25)
    
    return(return_res$gene)
  })
  
  names(genes_plot) <- unique(all_res$contrast)
  
  pdf(file.path(fig_dir, paste0("upset_plot_pseudobulk_", 
                                save_name, "_de.pdf")))
  print(upset(fromList(genes_plot), order.by = "freq"))
  dev.off()
  
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
  
  
  pdf(file.path(fig_dir, paste0(save_name, "_no_nd_all.pdf")))
  print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave,
                               gene_list = genes_plot$`NO - ND`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = FALSE))
  dev.off()
  
  fig_height <- round(length(genes_plot$`NO - ND`) / divide_by)
  
  pdf(file.path(fig_dir, paste0(save_name, "_no_nd_all_name.pdf")),
      width = 8, height = fig_height)
  print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave,
                               gene_list = genes_plot$`NO - ND`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = TRUE))
  dev.off()
  
  pdf(file.path(fig_dir, paste0(save_name, "_aab_nd_all.pdf")))
  print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave,
                               gene_list = genes_plot$`AAB - ND`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = FALSE))
  dev.off()
  
  fig_height <- round(length(genes_plot$`AAB - ND`) / divide_by)
  pdf(file.path(fig_dir, paste0(save_name, "_aab_nd_all_name.pdf")),
      width = 8, height = fig_height)
  print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave,
                               gene_list = genes_plot$`AAB - ND`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = TRUE))
  dev.off()
  
  
  pdf(file.path(fig_dir, paste0(save_name, "_no_aab_all.pdf")))
  print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave,
                               gene_list = genes_plot$`NO - AAB`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = FALSE))
  dev.off()
  
  fig_height <- round(length(genes_plot$`NO - AAB`) / divide_by)
  pdf(file.path(fig_dir, paste0(save_name, "_no_aab_all_name.pdf")),
      width = 8, height = fig_height)
  print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave,
                               gene_list = genes_plot$`NO - AAB`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = TRUE))
  dev.off()
  
  
  
  meta_ave_2 <- meta_ave %>%
    dplyr::filter(test_id == "diabetes_antigen")
  
  pdf(file.path(fig_dir, paste0(save_name, "_no_nd_all_genes_diabetes_plot.pdf")))
  print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave_2,
                               gene_list = genes_plot$`NO - ND`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = FALSE))
  dev.off()
  
  fig_height <- round(length(genes_plot$`NO - ND`) / divide_by)
  pdf(file.path(fig_dir, paste0(save_name, "_no_nd_all_genes_diabetes_plot_name.pdf")),
      width = 8, height = fig_height)
  print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave_2,
                               gene_list = genes_plot$`NO - ND`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = TRUE))
  dev.off()
  
  pdf(file.path(fig_dir, paste0(save_name, "_aab_nd_all_diabetes_plot.pdf")))
  print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave_2,
                               gene_list = genes_plot$`AAB - ND`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = FALSE))
  dev.off()
  
  fig_height <- round(length(genes_plot$`AAB - ND`) / divide_by)
  pdf(file.path(fig_dir, paste0(save_name, "_aab_nd_all_diabetes_plot_name.pdf")),
      width = 8, height = fig_height)
  print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave_2,
                               gene_list = genes_plot$`AAB - ND`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = TRUE))
  dev.off()
  
  pdf(file.path(fig_dir, paste0(save_name, "_no_aab_all_diabetes_plot.pdf")))
  print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave_2,
                               gene_list = genes_plot$`NO - AAB`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = FALSE))
  dev.off()
  
  fig_height <- round(length(genes_plot$`NO - AAB`) / divide_by)
  
  pdf(file.path(fig_dir, paste0(save_name, "_no_aab_all_diabetes_plot_name.pdf")),
      width = 8, height = fig_height)
  print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave_2,
                               gene_list = genes_plot$`NO - AAB`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = TRUE))
  dev.off()
  
  # new_meta <- meta_df %>%
  #   dplyr::select(sample, Status) %>%
  #   dplyr::distinct()
  # 
  # rownames(new_meta) <- new_meta$sample
  # 
  # new_meta$sample <- factor(new_meta$sample, levels = sample_order)
  # 
  # new_meta <- new_meta %>%
  #   arrange(sample)
  # 
  # plot_heatmap(seurat_object = seurat_object, gene_list = genes_plot$`NO - ND`,
  #              meta_df = new_meta, average_expression = TRUE, 
  #              meta_col = "sample", cluster_rows = TRUE, plot_rownames = FALSE,
  #              color_list = color_list)
  # dev.off()
  # plot_heatmap(seurat_object = seurat_object, gene_list = genes_plot$`NO - ND`,
  #              meta_df = new_meta, average_expression = TRUE, 
  #              meta_col = "sample", cluster_rows = TRUE, plot_rownames = FALSE,
  #              color_list = color_list, cluster_cols = TRUE)
  # 
  # plot only the diabetes antigen
  subset_meta <- meta_ave %>%
    dplyr::filter(test_id == "diabetes_antigen") 
  
  subset_meta$sample <- factor(subset_meta$sample, levels = sample_order)
  
  subset_meta <- subset_meta %>%
    dplyr::arrange(sample)
  
  
  pdf(file.path(fig_dir, paste0(save_name,
                                               "_no_nd_diabetes_antigen.pdf")))
  
  print(make_corrected_heatamp(vsd = vsd, meta_ave = subset_meta,
                               gene_list = genes_plot$`NO_diabetes_antigen - ND_diabetes_antigen`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = FALSE))
  
  dev.off()
  
  fig_height <- round(length(genes_plot$`NO_diabetes_antigen - ND_diabetes_antigen`) / divide_by)
  
  
  pdf(file.path(fig_dir, paste0(save_name,
                                               "_no_nd_diabetes_antigen_name.pdf")),
      width = 8, height = fig_height)
  
  print(make_corrected_heatamp(vsd = vsd, meta_ave = subset_meta,
                               gene_list = genes_plot$`NO_diabetes_antigen - ND_diabetes_antigen`,
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
                paste0(save_name, "_no_nd_diabetes_antigen_with_tet.pdf")))
  
  print(make_corrected_heatamp(vsd = vsd, meta_ave = subset_meta,
                               gene_list = genes_plot$`NO_diabetes_antigen - ND_diabetes_antigen`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = FALSE))
  
  grid::grid.newpage()
  
  print(make_corrected_heatamp(vsd = vsd, meta_ave = subset_meta2,
                               gene_list = genes_plot$`NO_diabetes_antigen - ND_diabetes_antigen`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = FALSE))
  
  dev.off()
  
  fig_height <- round(length(genes_plot$`NO_diabetes_antigen - ND_diabetes_antigen`) / divide_by)
  
  pdf(file.path(fig_dir,
                paste0(save_name, "_no_nd_diabetes_antigen_with_tet_name.pdf")),
      width = 8, height = fig_height)
  
  print(make_corrected_heatamp(vsd = vsd, meta_ave = subset_meta,
                               gene_list = genes_plot$`NO_diabetes_antigen - ND_diabetes_antigen`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = TRUE))
  
  grid::grid.newpage()
  
  print(make_corrected_heatamp(vsd = vsd, meta_ave = subset_meta2,
                               gene_list = genes_plot$`NO_diabetes_antigen - ND_diabetes_antigen`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = TRUE))
  
  dev.off()
  
}

fig_dir <- file.path(save_dir, "images", "testing_heatmaps")

ifelse(!dir.exists(fig_dir), dir.create(fig_dir), FALSE)

# De genes made in 07_de_analysis.R
vsd_bcells <- readRDS(file.path(save_dir, "rda_obj", 
                                "psuedobulk_capture_batch_corrected_object.rds"))
all_res_bcells <- read.csv(file.path(save_dir, "files", 
                                     "pseudobulk_all_bcell_de_capture_correct.csv"))

make_de_plots(all_res = all_res_bcells, save_name = "all_cells", 
              vsd = vsd_bcells, seurat_object = seurat_data, fig_dir = fig_dir)

so_no_ts <- subset(seurat_data, subset = sample != "JH310-12_TS")
# De genes made in 07_de_analysis.R
vsd_bcells <- readRDS(file.path(save_dir, "rda_obj", 
                                "no_ts_psuedobulk_capture_batch_corrected_object.rds"))
all_res_bcells <- read.csv(file.path(save_dir, "files", 
                                     "no_ts_pseudobulk_all_bcell_de_capture_correct.csv"))

make_de_plots(all_res = all_res_bcells, save_name = "no_ts", 
              vsd = vsd_bcells, seurat_object = so_no_ts, fig_dir = fig_dir)
