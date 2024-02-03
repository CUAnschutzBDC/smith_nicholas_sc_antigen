library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(DropletUtils)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log


cell_types <- "RNA_celltype"
clusters <- "RNA_cluster"
pval <- 0.05
logfc <- 0.5

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

#args <- c("JH313-15_PH", here("results"), "", here("files/sample_info.tsv"))

sample <- args[[1]]
sample <- gsub("__.*", "", sample)

results_dir <- args[[2]]

sample_info <- args[[4]]

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident

ADT_pca <- sample_info$adt_pca
run_adt_umap <- sample_info$adt_umap

RNA_pcs <- sample_info$PCs
resolution <- sample_info$resolution


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

# Make violin plots of ins, ia2, gad, tet, dna
# Make using the log 
# Separate by cell type and hash call

celltype_mapping <- c("Activated_memory" = "Bcell",
                      "CD14.Mono" = "Other",
                      "DN2" = "Bcell",
                      "BND2" = "Bcell",
                      "pDC" = "Other",
                      "Naive_1" = "Bcell",
                      "Plasmablast" = "Bcell",
                      "Resting_memory" = "Bcell",
                      "undetermined" = "Other",
                      "CD8.TCM" = "Other",
                      "Early_Memory" = "Bcell",
                      "Naive_3" = "Bcell",
                      "Early_memory" = "Bcell",
                      "CD8.TEM" = "Other",
                      "Naive_2" = "Bcell",
                      "Memory_IgE_IgG" = "Bcell",
                      "B.intermediate" = "Bcell",
                      "Memory_IgA" = "Bcell",
                      "CD4.TCM" = "Other",
                      "CD8.Naive" = "Other",
                      "gdT" = "Other",
                      "Activated_naive" = "Bcell",
                      "CD16.Mono" = "Other",
                      "NK" = "Other")

seurat_data$grouped_celltype <- celltype_mapping[seurat_data$RNA_celltype]
  
#table(is.na(seurat_data$grouped_celltype), seurat_data$RNA_celltype)

all_plots <- featDistPlot(seurat_data, geneset = c("INS-tet", "GAD-tet",
                                                   "IA2-tet", "TET-tet", 
                                                   "DNA-tet"),
                          sep_by = "grouped_celltype", assay = "TET",
                          combine = FALSE)

all_plots2 <- featDistPlot(seurat_data, geneset = c("INS-tet", "GAD-tet",
                                                   "IA2-tet", "TET-tet", 
                                                   "DNA-tet"),
                          sep_by = "grouped_celltype", assay = "SCAR_TET",
                          combine = FALSE)

all_plots3 <- featDistPlot(seurat_data, geneset = c("INS-tet", "GAD-tet",
                                                    "IA2-tet", "TET-tet", 
                                                    "DNA-tet"),
                           sep_by = "grouped_celltype", assay = "SCAR_TET_LOG",
                           combine = FALSE)


all_plots4 <- featDistPlot(seurat_data, geneset = c("INS-tet", "GAD-tet",
                                                    "IA2-tet", "TET-tet", 
                                                    "DNA-tet"),
                           sep_by = "grouped_celltype", assay = "SCAR_TET_LIBRA",
                           combine = FALSE)



all_plots5 <- featDistPlot(seurat_data, geneset = c("INS-tet", "GAD-tet",
                                                    "IA2-tet", "TET-tet", 
                                                    "DNA-tet"),
                           sep_by = "grouped_celltype", assay = "SCAR_TET_PROPORTIONS",
                           combine = FALSE)

all_plots5 <- lapply(all_plots5, function(x){
  x <- x +
    ggplot2::ylim(-1, 50)
})

# Save all of these as plots
final_plots <- lapply(names(all_plots), function(i){
  plot_1 <- all_plots[[i]] +
    ggplot2::ggtitle("tet_clr")
  
  plot_2 <- all_plots2[[i]] +
    ggplot2::ggtitle("scar_tet_clr")
  
  plot_3 <- all_plots3[[i]] +
    ggplot2::ggtitle("scar_tet_log")
  
  plot_4 <- all_plots4[[i]] +
    ggplot2::ggtitle("scar_tet_libra")
  
  plot_5 <- all_plots5[[i]] +
    ggplot2::ggtitle("scar_tet_proportions")
  
  cowplot::plot_grid(plot_1, plot_2, plot_3, plot_4, plot_5,
                     nrow = 3, ncol = 2)
  
})

pdf(file.path(save_dir, "images", "tetramers_b_cell_vs_other.pdf"),
    height = 12, width = 8)
print(final_plots)

dev.off()

seurat_b <- subset(seurat_data, subset = grouped_celltype == "Bcell")
seurat_other <- subset(seurat_data, subset = grouped_celltype == "Other")

all_data_b <- GetAssayData(seurat_b, assay = "SCAR_TET_LOG", slot = "data") %>% t
all_data_other <- GetAssayData(seurat_other, assay = "SCAR_TET_LOG", slot = "data") %>% t

find_cutoffs <- function(b_data, other_data, tetramer, quartile = 0.95){
  cutoff <- quantile(other_data[ , tetramer], quartile)
  
  
  mean_value <- mean(other_data[ , tetramer])
  sd_value <- sd(other_data[ , tetramer])
  
  # Define cutoff (e.g., mean + 3 * sd)
  cutoff2 <- mean_value + 3 * sd_value
  
  
  # passing cutoffs
  passing_cutoff <- b_data[ , tetramer] >= cutoff

  if(cutoff == 0){
    div_by <- 0.1
  } else {
    div_by <- cutoff
  }
  
  new_proportion <- b_data[, tetramer] / cutoff
  return(list(cutoff_tf = passing_cutoff, new_propotions = new_proportion))
}

all_tetramers <- colnames(all_data_b)

all_cutoffs <- lapply(all_tetramers, function(x){
  passing_cutoff <- find_cutoffs(b_data = all_data_b, 
                                 other_data = all_data_other,
                                 tetramer = x)
  
  return(passing_cutoff)
})

names(all_cutoffs) <- all_tetramers

all_proportions <- lapply(all_cutoffs, function(x){
  return(x$new_propotions)
})

all_proportions <- data.frame(all_proportions)

missing_row_names <- setdiff(colnames(seurat_data), rownames(all_proportions))

# Add rows of zeros for missing row names
if (length(missing_row_names) > 0) {
  zero_mat <- matrix(0, ncol = ncol(all_proportions), 
                     nrow = length(missing_row_names))
  rownames(zero_mat) <- missing_row_names
  colnames(zero_mat) <- colnames(all_proportions)
  all_proportions <- rbind(all_proportions, zero_mat)
  all_proportions <- all_proportions[order(match(colnames(seurat_data),
                                                 rownames(all_proportions))),]
}

all_proportions <- t(all_proportions)

all_cutoffs <- lapply(all_cutoffs, function(x){
  return(x$cutoff_tf)
})

all_cutoffs <- data.frame(all_cutoffs)


seurat_data[["NEW_TET_PROPORTIONS"]] <- CreateAssayObject(data = all_proportions)


get_column_name <- function(row) {
  true_columns <- names(row)[row]
  if (length(true_columns) == 0) {
    return("Negative")
  } else {
    return(paste(true_columns, collapse = "_"))
  }
}
full_tet_name <- apply(all_cutoffs, 1, get_column_name)


get_column_name2 <- function(row) {
  true_columns <- names(row)[row]
  if (length(true_columns) == 0) {
    return("Negative")
  } else if (length(true_columns) == 1){
    return(true_columns)
  } else if(length(true_columns) > 1 & 
            any(grepl("DNA.tet|TET.tet", true_columns))){
    return("Other_Multi_Reactive")
  } else {
    return("Islet_Multi_Reactive")
  }
}
tet_name <- apply(all_cutoffs, 1, get_column_name2)

new_data_frame <- data.frame("tet_name_cutoff" = tet_name, 
                             "full_tet_name_cutoff" = full_tet_name)

seurat_data <- AddMetaData(seurat_data, metadata = new_data_frame)

seurat_data$tet_name_cutoff[is.na(seurat_data$tet_name_cutoff)] <- "Non_b_cell"

seurat_data$full_tet_name_cutoff[is.na(seurat_data$full_tet_name_cutoff)] <- "Non_b_cell"

graphics.off()
pdf(file.path(save_dir, "images", "tetramer_new_vs_old_cutoff.pdf"))

cm1 <- confusionMatrix(seurat_data$tet_name_cutoff, seurat_data$scar_hash_id)
cm1 <- cm1 / rowSums(cm1)
print(pheatmap::pheatmap(cm1, main = "previous_vs_new_id"))

cm2 <- confusionMatrix(seurat_data$tet_name_cutoff, seurat_data$libra_tet_hash_id)
cm2 <- cm2 / rowSums(cm2)
print(pheatmap::pheatmap(cm2, main = "previous_vs_libra_id"))

dev.off()


all_plots6 <- featDistPlot(seurat_data, geneset = c("INS-tet", "GAD-tet",
                                                    "IA2-tet", "TET-tet", 
                                                    "DNA-tet"),
                           sep_by = "grouped_celltype",
                           col_by = "tet_name_cutoff",
                           assay = "TET",
                           combine = FALSE)

all_plots7 <- featDistPlot(seurat_data, geneset = c("INS-tet", "GAD-tet",
                                                    "IA2-tet", "TET-tet", 
                                                    "DNA-tet"),
                           sep_by = "grouped_celltype",
                           col_by = "tet_name_cutoff",
                           assay = "SCAR_TET",
                           combine = FALSE)

all_plots8 <- featDistPlot(seurat_data, geneset = c("INS-tet", "GAD-tet",
                                                    "IA2-tet", "TET-tet", 
                                                    "DNA-tet"),
                           sep_by = "grouped_celltype",
                           col_by = "tet_name_cutoff",
                           assay = "SCAR_TET_LOG",
                           combine = FALSE)

all_plots9 <- featDistPlot(seurat_data, geneset = c("INS-tet", "GAD-tet",
                                                    "IA2-tet", "TET-tet", 
                                                    "DNA-tet"),
                           sep_by = "grouped_celltype",
                           col_by = "tet_name_cutoff",
                           assay = "SCAR_TET_LIBRA",
                           combine = FALSE)

all_plots10 <- featDistPlot(seurat_data, geneset = c("INS-tet", "GAD-tet",
                                                    "IA2-tet", "TET-tet", 
                                                    "DNA-tet"),
                           sep_by = "grouped_celltype",
                           col_by = "tet_name_cutoff",
                           assay = "SCAR_TET_PROPORTIONS",
                           combine = FALSE)

all_plots10 <- lapply(all_plots10, function(x){
  x <- x +
    ggplot2::ylim(-1, 50)
})

# Save all of these as plots
final_plots2 <- lapply(names(all_plots), function(i){
  plot_1 <- all_plots6[[i]] +
    ggplot2::ggtitle("tet_clr")
  
  plot_2 <- all_plots7[[i]] +
    ggplot2::ggtitle("scar_tet_clr")
  
  plot_3 <- all_plots8[[i]] +
    ggplot2::ggtitle("scar_tet_log")
  
  plot_4 <- all_plots9[[i]] +
    ggplot2::ggtitle("scar_tet_libra")
  
  plot_5 <- all_plots10[[i]] +
    ggplot2::ggtitle("scar_tet_proportions")
  
  cowplot::plot_grid(plot_1, plot_2, plot_3, plot_4, plot_5,
                     nrow = 3, ncol = 2)
  
})

pdf(file.path(save_dir, "images", "tetramers_b_cell_vs_other_new_cutoff.pdf"),
    height = 12, width = 8)
print(final_plots2)

dev.off()

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))