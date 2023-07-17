library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(scuttle)
library(here)
library(openxlsx)
library(djvdj)


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "JH310-12_AP"

sample_info <- args[[4]]
#sample_info <- here("files/sample_info.tsv")

results_dir <- args[[2]]
#results_dir <- here("results")

sample_metadata <- args[[3]]
#sample_metadata <- here("files/Deidentified_donor_metadata.xlsx")
sample_metadata <- openxlsx::readWorkbook(sample_metadata, detectDates = TRUE)
colnames(sample_metadata) <- make.names(colnames(sample_metadata))
sample_metadata <- sample_metadata[!(is.na(sample_metadata$Sample.Name)),]
sample_metadata <- sample_metadata[sample_metadata$Sample.Name == sample,] %>%
  dplyr::select(!X.Key)

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

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

mt_pattern <- "^MT-" # "^MT-" for human, "^mt-" for mice

# Create seurat object ---------------------------------------------------------
seurat_object <- create_seurat_object(sample = sample,
                                      count_path = results_dir,
                                      ADT = ADT, hashtag = HTO,
                                      tenx_structure = "multi7",
                                      min_features = 200
                                      )

# Add mitochondrial percent
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,
                                                      pattern = mt_pattern)

sample_vect <- sample_metadata[1,,drop = TRUE]

seurat_object <- AddMetaData(seurat_object, metadata = sample_vect)

# Use scuttle for cutoffs ------------------------------------------------------
se <- as.SingleCellExperiment(seurat_object)
is.mito <- grep(mt_pattern, rownames(se))
per.cell <- perCellQCMetrics(se, subsets=list(Mito=is.mito))

qc.stats <- perCellQCFilters(per.cell,
                             sum.field = "sum",
                             detected.field = "detected",
                             sub.fields = "subsets_Mito_percent")

all_info <- cbind(per.cell, qc.stats)
colnames(all_info) <- paste0("cell_qc_", colnames(all_info))
all_info <- data.frame(all_info)
seurat_object <- AddMetaData(seurat_object, metadata = all_info)

rna_qual <- featDistPlot(seurat_object,
                          geneset = c("nFeature_RNA", "nCount_RNA",
                                      "percent.mt"),
                          plot_type = "violin", col_by = "cell_qc_discard")

if(ADT){
  adt_qual <- featDistPlot(seurat_object,
                            geneset = c("nFeature_ADT", "nCount_ADT",
                                        "percent.mt"),
                            plot_type = "violin", col_by = "cell_qc_discard")
}

pdf(file.path(save_dir, "images", "quality_plots.pdf"))
plot(rna_qual)
if(ADT){
  plot(adt_qual)
}
dev.off()

# Save before moving on
saveRDS(seurat_object, file = file.path(save_dir, "rda_obj",
                                        "seurat_unfilt.rds"))


seurat_object <- subset(seurat_object, subset = cell_qc_discard == FALSE)

# Normalization
# Single Cell Transform normalization
seurat_object <- SCTransform(seurat_object, vars.to.regress = vars.to.regress,
                             verbose = FALSE)

# Default normalization
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object) %>% 
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = vars.to.regress)

# Add in cell cycle
seurat_object <- CellCycleScoring(seurat_object, 
                                  g2m.features = cc.genes$g2m.genes,
                                  s.features = cc.genes$s.genes)


# Check against dropkick cells -------------------------------------------------
unfilt_path <- file.path(results_dir, sample, "outs", "multi", "count",
                         "raw_feature_bc_matrix")

full_seurat <- Read10X(unfilt_path)

full_seurat <- CreateSeuratObject(counts = full_seurat[["Gene Expression"]], 
                                  project = sample, min.cells = 0, 
                                  min.features = 0)

# Add mitochondrial percent
full_seurat[["percent.mt"]] <- PercentageFeatureSet(full_seurat,
                                                    pattern = mt_pattern)


dropkick_cells <- read.csv(file.path(save_dir, "files", "dropkick_cells.csv"),
                           row.names = 1)

full_seurat <- AddMetaData(full_seurat, dropkick_cells)

scuttle_meta <- seurat_object[[]] %>%
  dplyr::select(cell_qc_low_n_features, cell_qc_high_subsets_Mito_percent,
                cell_qc_low_lib_size) %>%
  dplyr::mutate(cell_qc_cell = "cell_qc_true")

full_seurat <- AddMetaData(full_seurat, scuttle_meta)

full_seurat$dropkick_label[is.na(full_seurat$dropkick_label)] <- "dropkick_false"
full_seurat$dropkick_label[full_seurat$dropkick_label == "True"] <- "dropkick_true"


full_seurat$cell_qc_cell[is.na(full_seurat$cell_qc_cell)] <- "cell_qc_false"

cm <- confusionMatrix(full_seurat$dropkick_label, full_seurat$cell_qc_cell)
cm <- cm / rowSums(cm)

comparison_heatmap <- pheatmap::pheatmap(cm)

full_seurat$cell_labels <- ifelse(full_seurat$dropkick_label == "dropkick_true" 
                                  & full_seurat$cell_qc_cell == "cell_qc_true",
                                  "both_cell", 
                                  ifelse(full_seurat$dropkick_label == "dropkick_true",
                                         "dropkick_only",
                                         ifelse(full_seurat$cell_qc_cell == "cell_qc_true",
                                                "cell_qc_cell", "no_cell")))


barplot_one <- scAnalysisR::stacked_barplots(full_seurat, 
                                             meta_col = "cell_labels")

rna_qual_two <- featDistPlot(full_seurat,
                             geneset = c("nFeature_RNA", "nCount_RNA",
                                         "percent.mt"),
                             plot_type = "violin", col_by = "cell_labels")

full_seurat <- subset(full_seurat, subset = cell_labels != "no_cell")

barplot_two <- scAnalysisR::stacked_barplots(full_seurat, 
                                             meta_col = "cell_labels")

pdf(file.path(save_dir, "images", "dropkick_vs_cellqc.pdf"))

print(barplot_one)
plot(rna_qual_two)
print(barplot_two)
grid::grid.newpage()
print(comparison_heatmap)

dev.off()

rm(full_seurat)

# VDJ --------------------------------------------------------------------------
# Read in VDJ data
vdj_dir <- file.path(results_dir, sample, "outs/per_sample_outs",
                     sample, "vdj_b")

seurat_object <- import_vdj(input = seurat_object,
                            vdj_dir = vdj_dir,
                            filter_paired = FALSE,
                            include_mutations = TRUE)

# Tetramers --------------------------------------------------------------------
# First pull out the tetramers
all_adts <- GetAssayData(seurat_object, assay = "ADT",
                         slot = "counts")

# First remove the incorrect insulin and other unnecessary ADTs
all_sums <- rowSums(all_adts)

# 1000 is arbitrary but should work to identify the actual ADTs. This should
# be checked and adjusted.
keep_adts <- all_sums[all_sums > 1000] 

all_adts <- all_adts[names(keep_adts),]

rownames(all_adts) <- gsub("INS-tet[a|b]", "INS-tet", rownames(all_adts))

seurat_object[["ADT"]] <- CreateAssayObject(counts = all_adts)

seurat_object <- NormalizeData(seurat_object, assay = "ADT",
                               normalization.method = "CLR", margin = 2)

all_tetramers <- all_adts[grepl("tet", rownames(all_adts)),]

seurat_object[["TET"]] <- CreateAssayObject(counts = all_tetramers)

seurat_object <- NormalizeData(seurat_object, assay = "TET",
                               normalization.method = "CLR", margin = 2)

if(HTO){
  # Demultiplex
  seurat_object <- HTODemux(seurat_object, assay = "HTO",
                            positive.quantile = 0.90)


  # Plots
  ridge_p <- RidgePlot(seurat_object, assay = "HTO",
                       features = rownames(stoc_data[["HTO"]])[1:2], ncol = 2)
  scatter_p <- FeatureScatter(seurat_object,
                              feature1 = "hto_Hashtag-STOC-86-008",
                              feature2 = "hto_Hashtag-STOC-86-009")
  Idents(seurat_object) <- "HTO_classification.global"
  vln_p <- VlnPlot(seurat_object, features = "nCount_RNA",
                   pt.size = 0.1, log = TRUE)

  pdf(file.path(save_dir, "images", "HTO_plots.pdf"))
  ridge_p
  scatter_p
  vln_p
  dev.off()
}

saveRDS(seurat_object, file = file.path(save_dir, "rda_obj",
                                        "seurat_start.rds"))
