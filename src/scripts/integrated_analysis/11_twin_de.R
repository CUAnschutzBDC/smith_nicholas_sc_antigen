library(here)
library(scAnalysisR)
library(pheatmap)
library(tidyverse)
library(MAST)
library(SingleCellExperiment)
library(Seurat)
library(cowplot)
library(SCOPfunctions)
library(KEGGREST)
library(org.Hs.eg.db)

source(here("src/scripts/muscat_plotting_functions.R"))

all_colors <- readRDS(file = file.path("files/all_colors.rds"))


final_colors <- all_colors$cell_type_colors

final_colors <- final_colors[names(final_colors) != "Translational_Intermediate"]

tetramer_colors <- all_colors$tetramer_colors

sample_colors <- all_colors$sample_colors


status_colors <- all_colors$status_colors

normalization_method <- "log" # can be SCT or log
# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

#args <- commandArgs(trailingOnly = TRUE)

args <- c("merged", here("results"), "", here("files/sample_info.tsv"))

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "merged"

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- args[[4]]
#sample_info <- here("files/sample_info.tsv")

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
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))
DefaultAssay(seurat_data) <- "RNA"

seurat_data <- DietSeurat(object = seurat_data, assays = "RNA")

seurat_data$Status <- factor(seurat_data$Status, levels = c("ND", "AAB", "T1D"))

sisters <- c("109", "108", "107")

#sisters <- c("109", "107")

seurat_data <- subset(seurat_data, subset = sample %in% twins)

# All B cells ------------------------------------------------------------------

## Rename columns for compatability --------------------------------------------

test_mapping <- c("DNA.tet" = "Other",
                  "GAD.tet" = "Islet_Reactive",
                  "IA2.tet" = "Islet_Reactive",
                  "INS.tet" = "Islet_Reactive",
                  "Islet_Multi_Reactive" = "Islet_Reactive",
                  "Negative" = "Other",
                  "Other_Multi_Reactive" = "Other",
                  "TET.tet" = "Other")

seurat_data$test_id <- test_mapping[seurat_data$tet_name_cutoff]

# Name columns for muscat sample_id, group_id, cluster_id
seurat_data$sample_id <- paste(seurat_data$sample, seurat_data$test_id,
                               sep = "_")

seurat_data$test_full_new_status <- paste(seurat_data$Status, 
                                          seurat_data$test_id,
                                          sep = "_")


Idents(seurat_data) <- "test_full_new_status"

tests_run <- list(c("T1D_Islet_Reactive", "ND_Islet_Reactive"))


all_markers <- lapply(tests_run, function(x){
  markers <- FindMarkers(seurat_data, test.use = "MAST", 
                         ident.1 = x[[1]],
                         ident.2 = x[[2]])
  
  markers$ident.1 <- x[[1]]
  markers$ident.2 <- x[[2]]
  markers$gene <- rownames(markers)
  
  return(markers)
  
})

all_markers <- do.call(rbind, all_markers)

markers_sig <- all_markers[all_markers$p_val_adj < 0.05,]

table(markers_sig$ident.1, markers_sig$ident.2)

markers_sig$cluster <- paste(markers_sig$ident.1,
                             markers_sig$ident.2,
                             sep = "_")
# 
# cluster_markers <- lapply(unique(seurat_data$RNA_celltype), function(x){
#   subset_seurat <- subset(seurat_data, subset = RNA_celltype == x)
#   
#   markers <- FindMarkers(seurat_data, test.use = "MAST", 
#                          ident.1 = "T1D_Islet_Reactive",
#                          ident.2 = "ND_Islet_Reactive")
#   
#   markers$ident.1 <- "T1D_Islet_Reactive"
#   markers$ident.2 <- "ND_Islet_Reactive"
#   markers$cluster <- x
#   markers$gene <- rownames(markers)
#   
#   return(markers)
# })
# 
# cluster_markers <- do.call(rbind, cluster_markers)
# 
# cluster_markers_sig <- cluster_markers[cluster_markers$p_val_adj < 0.05,]


# TODO save markers
write.csv(markers_sig, file.path(save_dir, "files", "twin_mast_de.csv"))

all_markers <- read.csv(file.path(save_dir, "files", "mast_de.csv"))

original_de <- all_markers[all_markers$cluster !=
                             "T1D_Islet_Reactive_AAB_Islet_Reactive",]$gene

overlapping_markers <- intersect(markers_sig$gene, original_de)
# overlapping_cluster_markers <- intersect(cluster_markers_sig$gene,
#                                          all_markers$gene)

# Make heatmap of all DE
de_genes <- markers_sig$gene

de_genes <- unique(de_genes)

sample_order <- seurat_data[[]] %>%
  dplyr::select(sample, Status) %>%
  dplyr::distinct() %>%
  dplyr::arrange(Status)

new_sample_order <- sisters

seurat_data$sample <- factor(seurat_data$sample, levels = new_sample_order)

heatmap_data <- plot_heatmap(seurat_data, gene_list = original_de,
                             colors = sample_colors, meta_col = "sample",
                             average_expression = TRUE, assay = "RNA",
                             plot_rownames = FALSE, cluster_rows = TRUE, 
                             return_data = TRUE)
blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", 
                "#2DB7A3", "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")

sample_info <- sample_order %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("sample")

rownames(sample_info) <- make.names(rownames(sample_info))

annotation_colors <- list("Status" = status_colors)


# This seems like a reasonalbe way to do it... 
# some things to consider
# Try making this plot but with the full dataset
# Make it so it's the same order as the full dataset
new_counts <- t(scale(t(heatmap_data$counts), center = FALSE, scale = TRUE))

heatmap <- pheatmap::pheatmap(new_counts, cluster_rows = TRUE, 
                              cluster_cols = FALSE, show_rownames = FALSE, 
                              show_colnames = TRUE, annotation_col = sample_info, 
                              annotation_colors = annotation_colors, 
                              color = blueYellow, border_color = NA, 
                              clustering_method = "complete", silent = TRUE)

ifelse(!dir.exists(file.path(save_dir, "images", "de")),
       dir.create(file.path(save_dir, "images", "de")), 
       FALSE)

pdf(file.path(save_dir, "images", "de", "de_heatmap_average.pdf"),
    width = 8, height = 12)

print(heatmap)

dev.off()


gost_res <- scAnalysisR::run_gost(seurat_de = markers_sig, organism = "hsapiens")

# TODO save gost
gost_output <- gost_res$gost_output$result

ifelse(!dir.exists(file.path(save_dir, "files", "de")),
       dir.create(file.path(save_dir, "files", "de")), 
       FALSE)

save_gost(gost_output = gost_output, save_dir_plots = file.path(save_dir, "images", "de"),
          save_dir_text = file.path(save_dir, "files", "de"))

kegg_output <- gost_output[gost_output$source == "KEGG", ]


gost_output[grepl("KEGG:04662", gost_output$term_id),]

gost_output[grepl("KEGG:04612", gost_output$term_id),]

gost_output[grepl("KEGG:05169", gost_output$term_id),]


# Make heatmaps of important pathways
path <- keggLink("pathway", "hsa")

all_genes <- path[grepl("hsa04662", path)]

all_gene_id <- keggConv( "ncbi-geneid", names(all_genes))

all_gene_id <- gsub("ncbi-geneid:", "", all_gene_id)

gene_symbols <- mapIds(org.Hs.eg.db, keys = all_gene_id, 
                       keytype = "ENTREZID", column = "SYMBOL")

gene_symbols <- gene_symbols[!grepl("LOC[0-9]+", gene_symbols)]

kegg_output[grepl("infection", kegg_output$term_name),]


gost_res <- scAnalysisR::run_gost(seurat_de = markers_sig,
                                  organism = "hsapiens", evcodes = TRUE)

save_gost(gost_output = gost_res, save_dir_plots = file.path(save_dir, "images", "de"),
          save_dir_text = file.path(save_dir, "files", "de"))

gost_output2 <- gost_res$gost_output$result


kegg_output2 <- gost_output2[gost_output2$source == "KEGG", ]

kegg_output2[grepl("infection", kegg_output2$term_name),]

find_overlap <- function(query){
  test_kegg <- kegg_output2[kegg_output2$query == query, ]
  
  gene_list <- split(test_kegg[["intersection"]], seq_len(nrow(test_kegg)))
  names(gene_list) <- test_kegg$term_name
  
  gene_list <- lapply(gene_list, function(x){
    strsplit(x, ",")
  })
  
  
  bcr_genes <- gene_list$`B cell receptor signaling pathway`
  epstein_bar_genes <- gene_list$`Epstein-Barr virus infection`
  
  other_infections <- gene_list[grepl("infection", names(gene_list))]
  
  
  bcr_overlap <- intersect(bcr_genes[[1]], epstein_bar_genes[[1]])
  percent_overlap <- length(bcr_overlap) / length(epstein_bar_genes[[1]]) * 100
  
  bcr_df <- data.frame(percent_overlap)
  rownames(bcr_df) <- "B cell receptor signaling pathway"
  
  all_overlap <- lapply(other_infections, function(x){
    overlap <- intersect(x[[1]], epstein_bar_genes[[1]])
    percent_overlap <- length(overlap) / length(epstein_bar_genes[[1]]) * 100
    return(percent_overlap)
    
  })
  
  all_overlap <- data.frame(unlist(all_overlap))
  
  colnames(all_overlap) <- "percent_overlap"
  
  all_overlap <- rbind(all_overlap, bcr_df)
  all_overlap$test <- query
  
  return(all_overlap)
  
}

total_overlaps <- lapply(unique(kegg_output2$query), function(x){
  find_overlap(x)
})


total_overlaps <- do.call(rbind, total_overlaps)

write.csv(total_overlaps, file.path(save_dir, "files", "de", "epstein_bar_gene_overlap.csv"))
