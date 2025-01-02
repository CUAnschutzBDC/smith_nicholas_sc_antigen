library(Seurat)
library(openxlsx)
library(here)
library(scAnalysisR)
library(tidyverse)

args <- c("merged", here("results"), "", here("files/sample_info.tsv"))

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "merged"

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- args[[4]]
#sample_info <- here("files/sample_info.tsv")

save_dir <- file.path(results_dir, "R_analysis", sample)

fig_dir <- file.path(save_dir, "images", "final_figures", "revision")
ifelse(!dir.exists(fig_dir), dir.create(fig_dir), FALSE)

markers_sig <- read.csv(file.path(save_dir, "files", "mast_de.csv"))

# Pull out sig genes
de_genes <- markers_sig[markers_sig$cluster !=
                          "T1D_Islet_Reactive_AAB_Islet_Reactive",]$gene

de_genes <- unique(de_genes)

gwas_genes <- openxlsx::readWorkbook(xlsxFile = here("files/t1d_gwas.xlsx"),
                                     colNames = TRUE, startRow = 3)

gwas_genes$gene <- gsub("_[0-9]+:.*", "", gwas_genes$Signal.name)

all_gwas <- unique(gwas_genes$gene)


seurat_data <- readRDS(file.path(save_dir, "rda_obj", 
                                 "seurat_processed_no_doublet.rds"))


hypergeometric <- function(de_genes, gwas_genes,
                           seurat_object){
  x <- length(intersect(gwas_genes, de_genes))
  m <- length(intersect(gwas_genes, rownames(seurat_object)))
  n <- length(setdiff(rownames(seurat_object), gwas_genes))
  total <- nrow(seurat_object)
  k <- length(de_genes)
  expected_num <- (m * k)/total
  representation <- x/expected_num
  p_val <- sum(dhyper(x:k, m, n, k))
  return_df <- data.frame(overlaps = x, 
                          gene_list_len = m, 
                          DE_length = k, 
                          background_len = n, 
                          total_genes = total,
                          expected_overlap = expected_num, 
                          overrepresentation = representation,
                          p_val = p_val)
}

hyper_df <- hypergeometric(de_genes, all_gwas, seurat_data)

write.csv(hyper_df, file.path(fig_dir, "gwas_hypergeometric.csv"))
intersect(de_genes, all_gwas)

seurat_data$sample_status <- paste(seurat_data$orig.ident, seurat_data$Status,
                                   sep = "_")
barplot_data <- stacked_barplots(seurat_data, meta_col = "tet_name_cutoff",
                            split_by = "Status")$data

write.csv(barplot_data, file.path(fig_dir, "status_binding_percents.csv"))

all_colors <- readRDS(file = file.path("files/all_colors.rds"))

status_levels <- c("ND", "AAB", "T1D")
status_colors <- all_colors$status_colors

seurat_data$Status <- factor(seurat_data$Status, levels = status_levels)

featDistPlot(seurat_data, geneset = "PTPN2", 
             sep_by = "Status",
             combine = FALSE, color = status_colors)[[1]] +
  ggplot2::theme_classic()
