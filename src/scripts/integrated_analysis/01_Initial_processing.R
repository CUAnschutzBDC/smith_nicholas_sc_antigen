# Next figure out how to add this to the snakemake pipeline
# Should run all other scripts first
# Maybe add in a "merged_scripts" section
# In getting scripts, if sample is merged, return all other scripts?

library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(here)
library(openxlsx)
library(ADTnorm)


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "merged"

sample_info <- args[[4]]
#sample_info <- here("files/sample_info.tsv")

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

samples_use <- sample_info[sample_info$sample != sample, ]$sample

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

vars.to.regress <- NULL

# Set directories
save_dir <- file.path(results_dir, "R_analysis", sample)

# Make directories
ifelse(!dir.exists(save_dir), dir.create(save_dir), FALSE)

ifelse(!dir.exists(file.path(save_dir, "images")),
       dir.create(file.path(save_dir, "images")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "files")),
       dir.create(file.path(save_dir, "files")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
       dir.create(file.path(save_dir, "rda_obj")), FALSE)

all_objs <- lapply(samples_use, function(x){
  so <- readRDS(file.path(results_dir, "R_analysis", x, 
                          "rda_obj", "seurat_processed.rds"))
  
  return(so)
})

# Merge Seurat objects
seurat_data <- merge(all_objs[[1]], all_objs[2:length(all_objs)],
                     add.cell.ids = samples_use)

# Remove any cells that aren't b cells
remove_celltypes <- c("CD14.Mono", "CD16.Mono", "CD4.TCM",
                      "CD8.Naive", "CD8.TCM", "CD8.TEM",
                      "gdT", "NK", "pDC")

`%notin%` <- Negate(`%in%`)
seurat_data <- subset(seurat_data, subset = RNA_celltype %notin% remove_celltypes)

seurat_data$old_status <- seurat_data$Status

status_mapping <- c("no" = "T1D",
                    "nd" = "ND",
                    "aab stage 1" = "AAB",
                    "aab stage 2" = "AAB")

seurat_data$Status <- status_mapping[seurat_data$old_status]

rm(all_objs)

seurat_data$sample <- seurat_data$orig.ident

seurat_data$sample <- factor(seurat_data$sample)

DefaultAssay(seurat_data) <- "RNA"

# Renormalize
seurat_data <- seurat_data %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

DefaultAssay(seurat_data) <- "CLR_ADT"

all_adts <- GetAssayData(seurat_data, assay = "CLR_ADT", slot = "counts") 
all_adts <- base::t(as.matrix(all_adts))

all_adts <- all_adts[ , !grepl("HTO", colnames(all_adts))]

seurat_data$celltype_cluster <- paste0(seurat_data$RNA_celltype,
                                       seurat_data$RNA_cluster,
                                       sep = "_")

sample_info <- seurat_data[[]] %>%
  dplyr::select(sample, celltype_cluster, nCount_ADT,
                RNA_celltype) %>%
  dplyr::mutate(batch = sample, 
                celltype_l1 = RNA_celltype,
                celltype_l2 = celltype_cluster)

# Run ADTnorm

# cell_x_adt_norm <- ADTnorm(
#   cell_x_adt = all_adts, 
#   cell_x_feature = sample_info, 
#   save_outpath = file.path(save_dir, "images", "ADTnorm"), 
#   study_name = "all_adts",
#   bimodal_marker = NULL, ## setting it to NULL will trigger ADTnorm to try different settings to find biomodal peaks for all the markers.
#   trimodal_marker = NULL, ## 
#   save_intermediate_fig = TRUE,
#   #exclude_zeroes = TRUE,
#   brewer_palettes = "Set1"
# )

bimodal_markers <- c("CD21", "IgD", "IgM", "CD27.1", "CXCR5.1",
                     "INS-tet", "GAD-tet",
                     "DNA-tet") 
trimodal_markers <- c()

process_markers <- colnames(all_adts)

process_markers <- process_markers[!process_markers %in% c("IgM")]

# CD21 failed in "JH313-15_NF", "JH313-15_JB", "JH313-15_DB", "JH313-15_OF"


cell_x_adt_norm <- ADTnorm(
  cell_x_adt = all_adts, 
  cell_x_feature = sample_info, 
  save_outpath = file.path(save_dir, "images", "ADTnorm"), 
  study_name = "all_adts",
  marker_to_process = process_markers, ## setting it to NULL by default will process all available markers in cell_x_adt.
  bimodal_marker = bimodal_markers, ## setting it to NULL will trigger ADTnorm to try different settings to find biomodal peaks for all the markers.
  trimodal_marker = trimodal_markers, ## 
  positive_peak = list(ADT = c("CXCR5.1", "CXCR5.1"),
                       sample = c("JH310-12_BH", "JH313-15_JB")),
  save_intermediate_fig = TRUE,
  #exclude_zeroes = TRUE,
  brewer_palettes = "Set1"
)

igm_norm <- ADTnorm(
  cell_x_adt = all_adts, 
  cell_x_feature = sample_info, 
  save_outpath = file.path(save_dir, "images", "ADTnorm"), 
  study_name = "all_adts",
  marker_to_process = "IgM", ## setting it to NULL by default will process all available markers in cell_x_adt.
  bimodal_marker = NULL, ## setting it to NULL will trigger ADTnorm to try different settings to find biomodal peaks for all the markers.
  trimodal_marker = "IgM", ## 
  # positive_peak = list(ADT = c("IgM", "IgM"),
  #                      sample = c("JH310-12_CP", "JH313-15_NF")),
  save_intermediate_fig = TRUE,
  brewer_palettes = "Set1",
  bw_smallest_bi = 1.5,
  detect_outlier_valley = TRUE,
  lower_peak_thres = 0.01,
  shoulder_valley = TRUE
)


all_norm <- cbind(cell_x_adt_norm, igm_norm)

#all_adt_obj <- CreateAssayObject(data = base::t(as.matrix(cell_x_adt_norm)))
all_adt_obj <- CreateAssayObject(data = base::t(as.matrix(all_norm)))

seurat_data[["ADT_norm"]] <- all_adt_obj

all_plots <- featDistPlot(seurat_data, geneset = colnames(all_norm), 
                          sep_by = "sample", combine = FALSE, 
                          assay = "DSB_ADT")

names(all_plots) <- colnames(all_norm)

DefaultAssay(seurat_data) <- "ADT_norm"

all_plots2 <- featDistPlot(seurat_data, geneset = colnames(all_norm), 
                           sep_by = "sample", combine = FALSE, 
                           assay = "ADT_norm")

names(all_plots2) <- colnames(all_norm)

combined_plots <- lapply(colnames(all_norm), function(x){
  plot_one <- all_plots[[x]] + 
    ggplot2::ggtitle("dsb_normalized")
  plot_two <- all_plots2[[x]] +
    ggplot2::ggtitle("ADTnorm")
  cowplot::plot_grid(plot_one, plot_two)
})

pdf(file.path(save_dir, "images", "ADTnorm", "all_violins.pdf"), width = 12)
print(combined_plots)

dev.off()


DefaultAssay(seurat_data) <- "RNA"

saveRDS(seurat_data, file = file.path(save_dir, "rda_obj",
                                      "seurat_adtnorm.rds"))

# Colors -----------------------------------------------------------------------
# final_colors <- c("Resting_memory" = "#924bdb", # Resting memory
#                   "Naive_1" = "#69ba3d", # Naive 1
#                   "Naive_2" = "#9a43a4", # Naive 2
#                   "Memory_IgE_IgG" = "#bf9b31", # Memory IgE/IgG1
#                   "Naive_3" = "#6477ce", # Naive 3
#                   "Memory_IgA" = "#d15131", # Memory IA
#                   "Early_memory" = "#4c9e8e", # Early Memory
#                   "BND2" = "#cc4570", #Bnd2
#                   "DN2" = "#648d4f", # DN2
#                   "Activated_memory" = "#985978", # Activated memory
#                   "Activated_naive" = "#a06846", # Activated naive
#                   "B.intermediate" = "#00008b",
#                   "CD14.Mono" = "#e0205a",
#                   "pDC" = "#ffb6d3",
#                   "Plasmablast" = "#ffac14",
#                   "CD8.TEM" = "#000000")

final_colors <- c("Resting_memory" = "#924bdb", # Resting memory
                  "Naive" = "#69ba3d", # Naive 1
                  "Activated_memory" = "#a06846", # Activated memory
                  "Translational_Intermediate" = "#4c9e8e",
                  "Plasmablast" = "#bf9b31") 

tetramer_full_colors <- MetBrewer::met.brewer(name = "Monet", n = 12,
                                         type = "continuous")

tetramer_full_colors <- tetramer_full_colors[c(1:9)]

names(tetramer_full_colors) <- c("INS-tet", "GAD-tet", "IA2-tet",
                                 "Islet_Reactive", "Negative",
                                 "DNA-tet", "TET-tet", "Other_Multi_Reactive",
                                 "Islet_Multi_Reactive")


sample_colors <- MetBrewer::met.brewer(name = "Archambault", n = 16,
                                       type = "continuous")


all_samples <- unique(seurat_data$sample)

names(sample_colors) <- all_samples

status_colors <- MetBrewer::met.brewer(name = "Hokusai1", n = 7)

status_colors <- status_colors[c(2, 4, 7)]

names(status_colors) <- c("T1D", "AAB", "ND")



all_colors <- list("cell_type_colors" = final_colors,
                   "tetramer_colors" = tetramer_full_colors,
                   "sample_colors" = sample_colors,
                   "status_colors" = status_colors)

saveRDS(all_colors, file = file.path("files/all_colors.rds"))
