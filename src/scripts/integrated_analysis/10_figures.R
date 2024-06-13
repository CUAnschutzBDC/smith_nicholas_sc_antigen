library(here)
library(scAnalysisR)
library(pheatmap)
library(tidyverse)
library(splitstackshape)
library(circlize)
library(viridis)
library(ggalluvial)
library(Seurat)
library(djvdj)
library(KEGGREST)
library(org.Hs.eg.db)
library(treemapify)
library(rstatix)
library(VennDiagram)

# Need to add VennDiagram to docker


source(here("src/scripts/integrated_analysis/figure_functions.R"))

# Naive - resting memory/activate naive/intermediate - ABC/memory - Plasmablast

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

seurat_data <- subset(seurat_data, subset = imcantation_isotype != "IGHE")

DefaultAssay(seurat_data) <- "RNA"

# Set up factors ---------------------------------------------------------------
status_levels <- c("ND", "AAB", "T1D")

seurat_data$Status <- factor(seurat_data$Status, levels = status_levels)

plotting_levels <- c("Negative", "DNA.tet", "TET.tet",
                     "Islet_Multi_Reactive", "IA2.tet", 
                     "GAD.tet", "INS.tet")

seurat_data$tet_name_cutoff <- factor(seurat_data$tet_name_cutoff, 
                                  levels = plotting_levels)

# Fix cell type names
seurat_data$final_celltype[is.na(seurat_data$final_celltype)] <- "Early_Memory"
seurat_data$final_celltype[seurat_data$final_celltype == "Resting_Memory"] <- "Early_Memory"

seurat_data$final_celltype <- factor(seurat_data$final_celltype,
                                     levels = c("Naive", "Early_Memory", "ABC",
                                                "Memory", "Plasmablast"))

image_dir <- file.path(save_dir, "images", "final_figures")

ifelse(!dir.exists(image_dir), dir.create(image_dir), FALSE)

# Set up long object -----------------------------------------------------------
sep_columns <- c("chains", "cdr3", "cdr3_length",
                 "cdr3_nt_length", "v_gene", "d_gene", "j_gene", "c_gene",
                 "reads", "umis", "productive", "full_length",
                 "v_ins", "v_del", "v_mis", "d_ins", "d_del",
                 "d_mis", "j_ins", "j_del", "j_mis", "c_ins", "c_del", "c_mis",
                 "all_ins", "all_del", "all_mis", "vd_ins", "vd_del", "dj_ins",
                 "dj_del", "v_mis_freq", "d_mis_freq", "j_mis_freq",
                 "c_mis_freq", "all_mis_freq")
keep_columns <- c("isotype", "final_celltype", "sample", "paired",
                  "clonotype_id", "Status", "tet_name_cutoff",
                  "all_chains")

all_info <- seurat_data[[]] %>%
  dplyr::mutate(all_chains = chains) %>%
  dplyr::select(dplyr::all_of(c(sep_columns, keep_columns)), n_chains) %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::filter(!is.na(n_chains))


all_info_split <- cSplit(all_info, sep_columns, sep = ";", direction = "long") %>%
  dplyr::filter(!is.na(chains))


# Colors -----------------------------------------------------------------------
all_colors <- readRDS(file = file.path("files/all_colors.rds"))


final_colors <- all_colors$cell_type_colors

final_colors <- final_colors[names(final_colors) != "Translational_Intermediate"]

tetramer_colors <- all_colors$tetramer_colors
names(tetramer_colors) <- make.names(names(tetramer_colors))

sample_colors <- all_colors$sample_colors


status_colors <- all_colors$status_colors

heavy <- unique(all_info_split[all_info_split$chains == "IGH",]$v_gene)
light <- unique(all_info_split[all_info_split$chains %in% c("IGL", "IGK"),]$v_gene)

palette1 <- colorRampPalette(colors = 
                               RColorBrewer::brewer.pal(name = "Set1", n = 9))

palette2 <- colorRampPalette(colors = 
                               RColorBrewer::brewer.pal(name = "Set2", n = 8))

heavy_colors <- palette1(length(heavy))
names(heavy_colors) <- heavy

light_colors <- palette2(length(light))
names(light_colors) <- light

all_colors <- c(heavy_colors, light_colors)

# Main figures -----------------------------------------------------------------

## Figure 1 --------------------------------------------------------------------

### 1C -------------------------------------------------------------------------
# Facs percents from Catherine
facs_percents <- read.table(here("files/PE_pos_sorted_percents.csv"),
                            header = TRUE, sep = ",")

# Rename status to match other figures
new_status <- c("healthy control" = "ND",
                "aab stage 1" = "AAB",
                "aab stage 2" = "AAB",
                "new-onset" = "T1D")

facs_percents$percent <- facs_percents$PE..sorted
facs_percents$old_status <- facs_percents$Status

facs_percents$Status <- new_status[facs_percents$old_status]

# Reorder to correct ordering
facs_percents$Status <- factor(facs_percents$Status, levels = unique(new_status))

facs_percents$type <- "PE_positive"

# Do t test
stat.test <- facs_percents %>%
  t_test(percent ~ Status) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "Status")

p1 <- ggplot2::ggplot(facs_percents, ggplot2::aes(x = Status, y = percent)) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = Status)) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.75)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggpubr::stat_pvalue_manual(
    stat.test, 
    label = "p.adj.signif"
  )

pdf(file.path(image_dir, "1C_ag_positive_facs.pdf"),
    width = 6, height = 4)
print(p1)
dev.off()


## Figure 2 --------------------------------------------------------------------
# Set up
# Make a column of cell type and cluster
seurat_data$celltype_cluster <- paste(seurat_data$final_celltype, 
                                      seurat_data$RNA_cluster,
                                      sep = "_")

# Merge some clusters based on gene expression
new_clusters <- c("Naive_0" = "Naive_0",
                  "Naive_1" = "Naive_0",
                  "Early_Memory_2" = "Early_Memory_1",
                  "Naive_3" = "Naive_2",
                  "Memory_4" = "Memory_3",
                  "Naive_5" = "Naive_4",
                  "Naive_6" = "Naive_0",
                  "ABC_7" = "ABC_5",
                  "Early_Memory_8" = "Early_Memory_6",
                  "Plasmablast_9" = "Plasmablast_8",
                  "Early_Memory_10" = "Early_Memory_7",
                  "Plasmablast_11" = "Plasmablast_8",
                  "Early_Memory_12" = "Early_Memory_1")

seurat_data$celltype_cluster <- new_clusters[seurat_data$celltype_cluster]


# Pick colors, Naive = Green, resting memory = purple/blue, plasmablast = orange,
# memory = brown
cluster_celltype_colors <- c("Naive_0" = "#58b44d",
                             "Early_Memory_1" = "#746ec8",
                             "Naive_2" = "#61732b",
                             "Memory_3" = "#9f7239",
                             "Naive_4" = "#b2b72e",
                             "ABC_5" = "#4bae8b",
                             "Early_Memory_6" = "#bd5abe",
                             "Early_Memory_7" = "#61a0d5",
                             "Plasmablast_8" = "#d89a40")

# Factor based on developmental order
seurat_data$celltype_cluster <- factor(seurat_data$celltype_cluster,
                                       levels = c("Naive_0", "Naive_2", 
                                                  "Naive_4",
                                                  "Early_Memory_1",
                                                  "Early_Memory_6", 
                                                  "Early_Memory_7",
                                                  "ABC_5",
                                                  "Memory_3",
                                                  "Plasmablast_8"))

# Reorder the colors to make sense
cluster_celltype_colors <- cluster_celltype_colors[order(match(names(cluster_celltype_colors),
                                                               levels(seurat_data$celltype_cluster)))]

### 2A -------------------------------------------------------------------------

# Make a UMAP of all cells colored by cell type
pdf(file.path(image_dir, "2A_rna_umap.pdf"), width = 8, height = 8)

print(plotDimRed(seurat_data, col_by = "celltype_cluster", 
                 plot_type = "rna_mnn.umap",
                 color = cluster_celltype_colors, ggrastr = TRUE))

dev.off()

### 2B -------------------------------------------------------------------------

# Celltype genes, provided by Catherine, ordered for aesthetics
gene_list <- c("CD38", "JCHAIN", "MZB1", "XBP1", "PRDM1", "IGHG1",
               "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "CD24", 
               "AHNAK", "CD27", "MS4A1", "CD19", "FGR", "FCRL5", 
               "HCK", "ITGAX", "TBX21", "CD86", "IRF4", "CD83", 
               "PAX5", "BACH2", "IL4R", "CD69", "CXCR4", "IGHM",
               "IGHD", "FCER2")

# Make a dot plot across the celltype cluster identities
pdf(file.path(image_dir, "2B_gene_dotplot.pdf"), width = 8, height = 10)
Idents(seurat_data) <- "celltype_cluster"
print(DotPlot(seurat_data, features = gene_list, 
              cols = c( "#D3D3D3", "#3416d9"),
              scale = TRUE, assay = "RNA") +
  ggplot2::coord_flip() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))

dev.off()

### 2C -------------------------------------------------------------------------

# Cell type by status
pdf(file.path(image_dir, "2C_status_celltype_barplot.pdf"),
    width = 8, height = 8)
barplot_data <- stacked_barplots(seurat_data, meta_col = "celltype_cluster",
                                 split_by = "Status",
                                 color = cluster_celltype_colors,
                                 return_values = TRUE)
barplot <- barplot_data$barplot +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

print(barplot)

dev.off()

# Get some stats
seurat_data$sample_status <- paste(seurat_data$orig.ident, seurat_data$Status,
                                   sep = "_")
barplot_data <- stacked_barplots(seurat_data, meta_col = "celltype_cluster",
                                 split_by = "sample_status",
                                 color = cluster_celltype_colors,
                                 return_values = TRUE)$data

barplot_data$sample <- gsub("_.*", "", barplot_data$split_by)
barplot_data$Status <- gsub(".*_", "", barplot_data$split_by)
barplot_data$Status <- factor(barplot_data$Status, 
                              levels = levels(seurat_data$Status))
barplot_data$meta_col <- factor(barplot_data$meta_col,
                                levels = levels(seurat_data$celltype_cluster))
barplot_data <- dplyr::ungroup(barplot_data)

all_t_tests <- lapply(levels(barplot_data$meta_col), function(celltype){
  barplot_celltype <- barplot_data[barplot_data$meta_col == celltype,]
  res <- t_test(barplot_celltype, percents ~ Status)
  res$celltype <- celltype
  return(res)

})

all_t_tests <- do.call(rbind, all_t_tests)

all_t_tests$p.adj <- NULL
all_t_tests$p.adj.signif <- NULL
all_t_tests$.y. <- NULL

all_t_tests <- adjust_pvalue(all_t_tests, p.col = "p",
                             output.col = "p_adj", method = "bonferroni")

write.csv(all_t_tests, file.path(image_dir, "cell_type_statistics.csv"))

colnames(barplot_data) <- c("celltype_cluster", "split_by", "Freq", "percents",
                            "sample", "Status")
barplot_data <- barplot_data[,colnames(barplot_data) != "split_by"]

write.csv(barplot_data, file.path(image_dir, "cell_type_percens.csv"))

## Figure 3 --------------------------------------------------------------------

### 3A -------------------------------------------------------------------------

# Barplot of antigen binding
pdf(file.path(image_dir, "3A_status_antigen_barplot.pdf"),
    width = 8, height = 8)
barplot <- stacked_barplots(seurat_data, meta_col = "tet_name_cutoff",
                            split_by = "Status",
                            color = tetramer_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

print(barplot)

dev.off()

# Get some stats
seurat_data$sample_status <- paste(seurat_data$orig.ident, seurat_data$Status,
                                   sep = "_")
barplot_data <- stacked_barplots(seurat_data, meta_col = "tet_name_cutoff",
                                 split_by = "sample_status",
                                 color = cluster_celltype_colors,
                                 return_values = TRUE)$data

barplot_data$sample <- gsub("_.*", "", barplot_data$split_by)
barplot_data$Status <- gsub(".*_", "", barplot_data$split_by)
barplot_data$Status <- factor(barplot_data$Status, 
                              levels = levels(seurat_data$Status))
barplot_data$meta_col <- factor(barplot_data$meta_col,
                                levels = levels(seurat_data$tet_name_cutoff))
barplot_data <- dplyr::ungroup(barplot_data)

all_t_tests <- lapply(levels(barplot_data$meta_col), function(tet){
  barplot_tet <- barplot_data[barplot_data$meta_col == tet,]
  res <- t_test(barplot_tet, percents ~ Status)
  res$tet <- tet
  
  return(res)
  
})

all_t_tests <- do.call(rbind, all_t_tests)

res_combined <- barplot_data %>%
  dplyr::mutate(type = ifelse(meta_col %in% c("Negative", "DNA.tet", 
                                              "TET.tet"),
                              meta_col, "IAR")) %>%
  dplyr::group_by(sample, type) %>%
  dplyr::mutate(percents_combined = sum(percents)) %>%
  dplyr::ungroup() %>%
  dplyr::select(sample, percents_combined, Status, type) %>%
  dplyr::distinct() %>%
  dplyr::filter(type == "IAR")

res <- t_test(res_combined, percents_combined ~ Status)
res$tet <- "IAR_combined"

all_t_tests <- rbind(all_t_tests, res)


all_t_tests$p.adj <- NULL
all_t_tests$p.adj.signif <- NULL
all_t_tests$.y. <- NULL

all_t_tests <- adjust_pvalue(all_t_tests, p.col = "p",
                              output.col = "p_adj", method = "bonferroni")

write.csv(all_t_tests, file.path(image_dir, "tetramer_statistics.csv"))

colnames(barplot_data) <- c("tetramer", "split_by", "Freq", "percents",
                            "sample", "Status")
barplot_data <- barplot_data[,colnames(barplot_data) != "split_by"]

write.csv(barplot_data, file.path(image_dir, "tetramer_percents.csv"))

colnames(res_combined) <- c("sample", "percents", "Status", "type") 
write.csv(res_combined, file.path(image_dir, "IAR_percents.csv"))


### 3B -------------------------------------------------------------------------

# Umaps by tetramer coloring
all_umaps <- lapply(unique(seurat_data$tet_name_cutoff), function(x){
  umap_return <- plotDimRed(seurat_data, col_by = "celltype_cluster", 
                            plot_type = "rna_mnn.umap",
                            color = cluster_celltype_colors, ggrastr = TRUE,
                            highlight_group = TRUE, group = x,
                            meta_data_col = "tet_name_cutoff")[[1]] +
    ggplot2::theme(legend.position = "none")
  
  return(umap_return)
})

names(all_umaps) <- unique(seurat_data$tet_name_cutoff)

full_umaps <- cowplot::plot_grid(all_umaps$INS.tet, all_umaps$IA2.tet,
                                 all_umaps$GAD.tet, 
                                 all_umaps$Islet_Multi_Reactive,
                                 all_umaps$TET.tet, all_umaps$DNA.tet,
                                 all_umaps$Negative,
                                 ncol = 4)

pdf(file.path(image_dir, "3B_tetramer_umap.pdf"),
    height = 8, width = 16)
print(full_umaps)
dev.off()

### 3C -------------------------------------------------------------------------


# Make a barplot of antigen 
make_barplot <- function(seurat_object, tetramer_colors_use,
                         meta_data_col = "tet_name_cutoff"){
  seurat_object$full_split <- paste(seurat_object$final_celltype,
                                    seurat_object$Status, sep = "__")
  
  barplot <- stacked_barplots(seurat_object, meta_col = meta_data_col,
                              split_by = "full_split",
                              color = tetramer_colors_use,
                              return_values = TRUE)
  
  barplot_data <- barplot$data
  barplot_data$celltype <- gsub("__.*", "", barplot_data$split_by)
  barplot_data$Status <- gsub(".*__", "", barplot_data$split_by)
  barplot_data[[meta_data_col]] <- barplot_data$meta_col
  
  barplot_data$Status <- factor(barplot_data$Status, 
                                levels = status_levels)
  
  barplot_data$celltype <- factor(barplot_data$celltype,
                                  levels = c("Naive", "Early_Memory", "ABC",
                                             "Memory", "Plasmablast"))
  
  barplot <- ggplot2::ggplot(barplot_data, aes(x = Status, y = percents,
                                               fill = !!sym(meta_data_col))) +
    ggplot2::geom_bar(position = "stack", stat = "identity") +
    ggplot2::facet_grid(~celltype, switch = "x") +
    ggplot2::theme(strip.placement = "outside",
                   strip.background = element_rect(fill = NA, color = "white"),
                   panel.spacing = unit(-.01,"cm"),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::scale_fill_manual(values = tetramer_colors_use)
  
  return(barplot)
  
}

`%notin%` <- Negate(`%in%`)
tetramer_colors_use <- tetramer_colors[names(tetramer_colors) %notin%
                                         c("Other_Multi_Reactive", 
                                           "Islet_Reactive")]

barplot <- make_barplot(seurat_object = seurat_data,
                        tetramer_colors_use = tetramer_colors_use)

pdf(file.path(image_dir, "3C_tetramer_stacked_barplot.pdf"),
    height = 8, width = 8)
print(barplot)
dev.off()

### 3D -------------------------------------------------------------------------

seurat_islet <- subset(seurat_data, 
                       subset = tet_name_cutoff == "Islet_Multi_Reactive")

# Figure out the percent that contains each type
# Want a data frame that has 
# count percent  tetramer status

all_tetramer <- c("INS.tet", "IA2.tet", "GAD.tet", "DNA.tet", "TET.tet")
all_meta <- seurat_islet[[c("full_tet_name_cutoff", "Status")]]
all_meta <- all_meta %>%
  dplyr::group_by(Status) %>%
  dplyr::add_count(name = "total_cells")

all_dfs <- lapply(all_tetramer, function(tetramer){
  count_data <- all_meta %>%
    dplyr::filter(grepl(tetramer, full_tet_name_cutoff)) %>%
    dplyr::group_by(Status) %>%
    dplyr::add_count(name = "tetramer_cells") %>%
    dplyr::mutate(tetramer = tetramer, 
                  percent = tetramer_cells / total_cells * 100) %>%
    dplyr::select(-full_tet_name_cutoff) %>%
    dplyr::distinct()
  
  return(count_data)
    
})

all_dfs <- do.call(rbind, all_dfs)

all_dfs$tetramer <- factor(all_dfs$tetramer, levels = plotting_levels)

# Now make a barplot
tet_only_colors <- tetramer_colors[all_tetramer]
barplot2 <- ggplot2::ggplot(all_dfs, ggplot2::aes(x = Status,
                                                 y = percent,
                                                 fill = tetramer)) +
  ggplot2::geom_bar(stat = "identity", position = "dodge") +
  ggplot2::scale_fill_manual(values = tet_only_colors)

pdf(file.path(image_dir, "3D_islet_multi_tetramer_barplot.pdf"),
    height = 8, width = 8)
print(barplot2)
dev.off()

### 3E -------------------------------------------------------------------------

# Make a violin plot that is all cells for each antigen reactivity across
# status. Here each box will be cells from one antigen reactivity so nothing
# should be below 1.

# Start by adding the binding information to the meta data
seurat_meta <- seurat_data[[]]

score_data <- GetAssayData(seurat_data, slot = "data", assay = "NEW_TET_PROPORTIONS") %>%
  as.matrix() %>%
  t() %>%
  data.frame

full_meta <- merge(seurat_meta, score_data, by = "row.names")[c("Sample.Name",
                                                                "Status",
                                                                "tet_name_cutoff",
                                                                "INS.tet",
                                                                "TET.tet",
                                                                "IA2.tet",
                                                                "GAD.tet",
                                                                "DNA.tet")]

# Now I want to select out each single reactivity and the matching score
# column
reactivities_test <- c("INS.tet", "IA2.tet", "GAD.tet", "DNA.tet", "TET.tet")

all_reactivities <- lapply(reactivities_test, function(x){
  new_meta <- full_meta[full_meta$tet_name_cutoff == x, 
                        c("Status", "tet_name_cutoff", x)]
  colnames(new_meta) <- c("Status", "tet_name_cutoff", "score")
  return(new_meta)
})

all_reactivities <- do.call(rbind, all_reactivities)

final_plot <- ggplot2::ggplot(data = all_reactivities,
                                  ggplot2::aes(x = Status,
                                               y = score,
                                               fill = Status)) +
  ggplot2::geom_violin(scale = "width", position = ggplot2::position_dodge(1)) +
  ggplot2::facet_grid(~tet_name_cutoff, switch = "x") +
  ggplot2::scale_fill_manual(values = status_colors) + 
  ggplot2::stat_summary(fun = median, geom = "point", size = 2,
                        position = ggplot2::position_dodge(1))

# I didn't really like that. I prefer to see the distribution of binding 
# across all cells regardles of the call. Let's add a line at 1 indicating
# the cutoff

all_reactivities <- lapply(reactivities_test, function(x){
  new_meta <- full_meta[ , 
                        c("Status", "tet_name_cutoff", x)]
  colnames(new_meta) <- c("Status", "tet_name_cutoff", "score")
  new_meta$tetramer <- x
  return(new_meta)
})

all_reactivities <- do.call(rbind, all_reactivities)

all_reactivities$tetramer <- factor(all_reactivities$tetramer,
                                    levels = c("INS.tet", "IA2.tet",
                                               "GAD.tet", "TET.tet",
                                               "DNA.tet"))

final_plot <- ggplot2::ggplot(data = all_reactivities,
                              ggplot2::aes(x = Status,
                                           y = score,
                                           fill = Status)) +
  ggplot2::geom_violin(scale = "width", position = ggplot2::position_dodge(1)) +
  ggplot2::facet_grid(~tetramer, switch = "x") +
  ggplot2::scale_fill_manual(values = status_colors) + 
  ggplot2::stat_summary(fun = median, geom = "point", size = 2,
                        position = ggplot2::position_dodge(1)) 
  #ggplot2::geom_hline(linetype = "dotted", yintercept = 1)


pdf(file.path(image_dir, "3E_violin_tetramer_scores.pdf"),
    height = 4, width = 10)
print(final_plot)
dev.off()

## Figure 4 --------------------------------------------------------------------
# Read in de genes
markers_sig <- read.csv(file.path(save_dir, "files", "mast_de.csv"))
twin_markers <- read.csv(file.path(save_dir, "files", "twins_mast_de.csv"))
sister_markers <- read.csv(file.path(save_dir, "files", "sisters_mast_de.csv"))

# Pull out sig genes
de_genes <- markers_sig[markers_sig$cluster !=
                          "T1D_Islet_Reactive_AAB_Islet_Reactive",]$gene

de_genes <- unique(de_genes)

### 4A -------------------------------------------------------------------------
# Make heatmap of all DE
new_sample_order <- c("110", "116", "108", "107", "113", "114", "118",
                      "106", "117", "115", "105", "111", "102", "112",
                      "119", "109")

seurat_data$sample <- factor(seurat_data$sample, levels = new_sample_order)

seurat_islet <- subset(seurat_data, subset = tet_name_cutoff %in%
                         c("Islet_Multi_Reactive", "IA2.tet",
                           "GAD.tet", "INS.tet"))

heatmap_data <- plot_heatmap(seurat_islet, gene_list = de_genes,
                             colors = sample_colors, meta_col = "sample",
                             average_expression = TRUE, assay = "RNA",
                             plot_rownames = FALSE, cluster_rows = TRUE, 
                             return_data = TRUE)
blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", 
                "#2DB7A3", "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")

sample_order <- seurat_data[[]] %>%
  dplyr::select(sample, Status) %>%
  dplyr::distinct() %>%
  dplyr::arrange(Status)

sample_info <- sample_order %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("sample")

rownames(sample_info) <- make.names(rownames(sample_info))

annotation_colors <- list("Status" = status_colors)

heatmap <- pheatmap::pheatmap(heatmap_data$z_score, cluster_rows = TRUE, 
                              cluster_cols = FALSE, show_rownames = FALSE, 
                              show_colnames = TRUE, annotation_col = sample_info, 
                              annotation_colors = annotation_colors, 
                              color = blueYellow, border_color = NA, 
                              clustering_method = "complete", silent = TRUE)

graphics.off()

pdf(file.path(image_dir, "4A_de_heatmap_average.pdf"),
    width = 8, height = 12)

print(heatmap)

dev.off()

### 4B -------------------------------------------------------------------------

# Use jaccard distance to make a dendrogram of gse terms
gse_res <- read.csv(file.path(save_dir, "files/de/all_GSE_results.csv"))

name_mapping <- c("T1D" = "T1D_Islet_Reactive_ND_Islet_Reactive",
                  "AAB" = "AAB_Islet_Reactive_ND_Islet_Reactive")

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}




keep_pathways <- list("pro_inflammatory_signaling" =
                        c("KEGG:04657", "KEGG:04668", "KEGG:04064"),
                      "bacterial_infection_and_signaling" = 
                        c("KEGG:05135", "KEGG:05131", "KEGG:05130",
                          "KEGG:05132"),
                      "general_signaling" = 
                        c("KEGG:04062", "KEGG:04910", "KEGG:04010",
                          "KEGG:05415", "KEGG:04664", "KEGG:04662",
                          "KEGG:04666"),
                      "viral_infection_and_signaling" =
                        c("KEGG:04620", "KEGG:05161", "KEGG:05160",
                          "KEGG:05167", "KEGG:05162", "KEGG:05169",
                          "KEGG:05164", "KEGG:04622", "KEGG:04621"),
                      "autoimmune_disease" = 
                        c("KEGG:05323", "KEGG:05321", "KEGG:05320"),
                      "antigen_processing_presentation" =
                        c("KEGG:04612", "KEGG:04145"))

pathway_colors <- list("pro_inflammatory_signaling" = "#00b7a7",
                      "bacterial_infection_and_signaling" =  "#721b3e",
                      "general_signaling" =  "#007e2f",
                      "viral_infection_and_signaling" = "#b86092",
                      "autoimmune_disease" =  "#a40000",
                      "antigen_processing_presentation" = "#ffcd12")

gse_res <- gse_res[gse_res$term_id %in% unlist(keep_pathways),]

hc_res_t1d <- make_histogram(query = name_mapping[[1]], 
                             gse_res = gse_res)


hc_res_aab <- make_histogram(query = name_mapping[[2]], 
                             gse_res = gse_res)

hist1 <- plot_histogram(hc = hc_res_t1d$full, ylim = c(2, -1),
                        title = "T1D vs ND pathways",
                        gse_res = hc_res_t1d$gse_res)

hist2 <- plot_histogram(hc = hc_res_aab$full, ylim = c(2.1, -1),
                        title = "AAB vs ND pathways",
                        gse_res = hc_res_aab$gse_res)

# Here I figure out the ranges of each so I can replot with keys
# that match (and only need 1 key rather than 2)

min_p_val <- min(-log(hist1$layers[[3]]$data$p_value),
                 -log(hist1$layers[[3]]$data$p_value)) - 0.5

max_p_val <- max(-log(hist1$layers[[3]]$data$p_value),
                 -log(hist1$layers[[3]]$data$p_value)) + 0.5


min_size <- min(hist1$layers[[3]]$data$intersection_size,
                hist1$layers[[3]]$data$intersection_size)


max_size <- max(hist1$layers[[3]]$data$intersection_size,
                hist1$layers[[3]]$data$intersection_size)

# Now that I know the ranges, I can sync them ane remove the legend from 
# the first plot
hist1_scaled <- hist1 + 
  ggplot2::scale_color_gradient(limits = c(min_p_val, max_p_val),
                                low = "blue", high = "red") +
  ggplot2::scale_size(limits = c(min_size, max_size)) +
  ggplot2::guides(color = "none", size = "none") 

hist2_scaled <- hist2 + 
  ggplot2::scale_color_gradient(limits = c(min_p_val, max_p_val),
                                low = "blue", high = "red") +
  ggplot2::scale_size(limits = c(min_size, max_size))

combined <- cowplot::plot_grid(hist1_scaled, hist2_scaled, nrow = 1,
                   ncol = 2, rel_widths = c(0.8, 1))

pdf(file.path(save_dir, "images", "final_figures",
              "4B_kegg_clustering.pdf"), 
    height = 8, width = 18) 

print(combined)


dev.off()

### 4C -------------------------------------------------------------------------

# DE volcano plots
# Make volcano for all DE genes --> rerun DE for all
# Color volcano by each set of genes
all_markers <- read.csv(file.path(save_dir, "files", "mast_de_all.csv"))

de_genes <- all_markers[all_markers$ident.1 == "T1D_Islet_Reactive" & 
                          all_markers$ident.2 == "ND_Islet_Reactive",]

all_volcanos <- lapply(names(keep_pathways), function(group_name){
  kegg_pathways <- keep_pathways[[group_name]]
  color <- pathway_colors[[group_name]]
  
  return_plot <- make_volcano(de_genes = de_genes, 
                              kegg_pathways = kegg_pathways,
                              group_name = group_name,
                              color = color)
  
  return(return_plot)
})

final_volcanos <- cowplot::plot_grid(plotlist = all_volcanos, 
                                     nrow = 2, ncol = 3)

pdf(file.path(save_dir, "images", "final_figures",
              "4C_Kegg_volcanos.pdf"), 
    height = 6, width = 12) 

print(final_volcanos)
dev.off()

# Interactive volcanos
plotly_volcanos <- lapply(names(keep_pathways), function(group_name){
  kegg_pathways <- keep_pathways[[group_name]]
  color <- pathway_colors[[group_name]]
  
  return_plot <- make_volcano(de_genes = de_genes, 
                              kegg_pathways = kegg_pathways,
                              group_name = group_name,
                              color = color,
                              plotly = TRUE)
  
  return_plot <- plotly::ggplotly(return_plot, tooltip = "text")
  
  htmlwidgets::saveWidget(return_plot, file.path(save_dir, "images",
                                                 "final_figures",
                                                 paste0(group_name, 
                                                 "_interactive_volcano.html")))
  return(return_plot)
})

### 4D -------------------------------------------------------------------------

sample_levels <- levels(seurat_data$sample)
seurat_data$sample <- as.character(seurat_data$sample)
seurat_data$sample <- factor(seurat_data$sample, levels = sample_levels)

# Umap highlighting only the twins
twin_umap <- plotDimRed(seurat_data, col_by = "sample", color = sample_colors,
                        plot_type = "pca.umap", highlight_group = TRUE,
                        group = c("107", "108", "109"),
                        meta_data_col = "sample", ggrastr = TRUE,
                        reorder_cells = TRUE)[[1]] +
  ggplot2::ylab("UMAP 2") +
  ggplot2::xlab("UMAP 1")


pdf(file.path(save_dir, "images", "final_figures",
              "4D_twin_UMAP.pdf"), 
    height = 6, width = 6) 

print(twin_umap)
dev.off()

### 4E -------------------------------------------------------------------------

# Pull out sig genes
twin_de_genes <- twin_markers[twin_markers$cluster !=
                          "T1D_Islet_Reactive_AAB_Islet_Reactive",]$gene

twin_de_genes <- unique(twin_de_genes)

sister_de_genes <- sister_markers[sister_markers$cluster !=
                                "T1D_Islet_Reactive_AAB_Islet_Reactive",]$gene

sister_de_genes <- unique(sister_de_genes)

de_genes <- markers_sig[markers_sig$cluster == 
                          "T1D_Islet_Reactive_ND_Islet_Reactive",]$gene

de_list <- list(all_de = unique(de_genes), twin_de = twin_de_genes,
                sister_de = sister_de_genes)

venn_plot <- VennDiagram::venn.diagram(x = de_list,
                          category.names = names(de_list),
                          output = TRUE,
                          filename = NULL,
                          col=c("#440154ff", '#21908dff', '#fde725ff'),
                          fill = c(alpha("#440154ff",0.3),
                                   alpha('#21908dff',0.3),
                                   alpha('#fde725ff',0.3)),)


pdf(file.path(image_dir, "4E_de_venn.pdf"), height = 6, width = 6)
grid::grid.draw(venn_plot)
dev.off()

pdf(file.path(image_dir, "4E_de_upset.pdf"), height = 6, width = 6)
print(UpSetR::upset(UpSetR::fromList(de_list), order.by = "freq"))
dev.off()

## Figure 5 --------------------------------------------------------------------
# Make a barplot of antigen 

### 5A -------------------------------------------------------------------------

# This function could be easily modified to be much more flexible and work
# in both locations where the function has the same name
make_barplot <- function(seurat_object, isotype_colors,
                         meta_data_col = "imcantation_isotype"){
  seurat_object$full_split <- paste(seurat_object$Status,
                                    seurat_object$tet_name_cutoff, sep = "__")
  
  barplot <- stacked_barplots(seurat_object, meta_col = meta_data_col,
                              split_by = "full_split",
                              color = isotype_colors,
                              return_values = TRUE)
  
  barplot_data <- barplot$data
  barplot_data$Status <- gsub("__.*", "", barplot_data$split_by)
  barplot_data$tet_name_cutoff <- gsub(".*__", "", barplot_data$split_by)
  barplot_data[[meta_data_col]] <- barplot_data$meta_col
  
  barplot_data$tet_name_cutoff <- factor(barplot_data$tet_name_cutoff,
                                levels = plotting_levels)

  barplot_data$Status <- factor(barplot_data$Status,
                                  levels = levels(seurat_data$Status))
  
  barplot <- ggplot2::ggplot(barplot_data, aes(x = tet_name_cutoff, y = percents,
                                               fill = !!sym(meta_data_col))) +
    ggplot2::geom_bar(position = "stack", stat = "identity") +
    ggplot2::facet_grid(~Status, switch = "x") +
    ggplot2::theme(strip.placement = "outside",
                   strip.background = element_rect(fill = NA, color = "white"),
                   panel.spacing = unit(-.01,"cm"),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::scale_fill_manual(values = isotype_colors)
  
  return(barplot)
  
}
# 5A stacked bar chart of isotype distribution by cell type
isotype <- subset(seurat_data, subset = imcantation_isotype != "")

colors <- MetBrewer::met.brewer(name = "Peru1", n = 6)
colors2 <- MetBrewer::met.brewer(name = "Peru2", n = 7)
isotype_colors <- c(colors[c(5, 6, 3, 4)], colors2[4:7])
names(isotype_colors) <- c("IGHM", "IGHD",  "IGHA1", "IGHA2",
                           "IGHG1", "IGHG2", "IGHG3", "IGHG4")
isotype$imcantation_isotype <- factor(isotype$imcantation_isotype,
                                      levels = names(isotype_colors))


barplot <- make_barplot(isotype, isotype_colors)

pdf(file.path(image_dir, "5A_stacked_bar_chart_isotype_status.pdf"),
    width = 8, height = 8)

print(barplot)

dev.off()
graphics.off()

### 5B -------------------------------------------------------------------------

# 5B CDR3 length
sep_columns <- c("chains", "cdr3", "cdr3_length",
                 "cdr3_nt_length", "v_gene", "d_gene", "j_gene", "c_gene",
                 "reads", "umis", "productive", "full_length",
                 "v_ins", "v_del", "v_mis", "d_ins", "d_del",
                 "d_mis", "j_ins", "j_del", "j_mis", "c_ins", "c_del", "c_mis",
                 "all_ins", "all_del", "all_mis", "vd_ins", "vd_del", "dj_ins",
                 "dj_del", "v_mis_freq", "d_mis_freq", "j_mis_freq",
                 "c_mis_freq", "all_mis_freq")
keep_columns <- c("isotype", "final_celltype", "sample", "paired",
                  "clonotype_id", "Status", "tet_name_cutoff",
                  "all_chains")

all_info <- seurat_data[[]] %>%
  dplyr::mutate(all_chains = chains) %>%
  dplyr::select(dplyr::all_of(c(sep_columns, keep_columns)), n_chains) %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::filter(!is.na(n_chains))


all_info_split <- cSplit(all_info, sep_columns, sep = ";", direction = "long") %>%
  dplyr::filter(!is.na(chains))


# CDR3 length

all_info_split$tet_name_cutoff <- factor(all_info_split$tet_name_cutoff,
                                         levels = plotting_levels)

heavy_data <- all_info_split %>%
  dplyr::filter(chains == "IGH")

boxplot_cdr3_len <- ggplot2::ggplot(heavy_data, ggplot2::aes(y = cdr3_length,
                                                             x = tet_name_cutoff,
                                                             fill = Status)) +
  ggplot2::geom_boxplot(size = 0.1, outlier.size = 0.5) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1))

pdf(file.path(image_dir, "5B_cdr3_length.pdf"),
    width = 8, height = 4)

print(boxplot_cdr3_len)

dev.off()
graphics.off()

### 5C -------------------------------------------------------------------------

# Find charge of CDR3 with alakazam
heavy_data <- all_info_split %>%
  dplyr::filter(chains == "IGH")

heavy_data$charge <- alakazam::charge(seq = heavy_data$cdr3)

# Barplot, y = charge, x = tetramer, color = status
h_charge <- ggplot2::ggplot(heavy_data, ggplot2::aes(y = charge,
                                                     x = tet_name_cutoff,
                                                     fill = Status)) +
  ggplot2::geom_boxplot(size = 0.1, outlier.size = 0.5) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1)) +
  ggplot2::ggtitle("Heavy chain SMH")

pdf(file.path(image_dir, "5C_heavy_charge.pdf"),
    width = 8, height = 4)

print(h_charge)

dev.off()
graphics.off()

# Repeat with light
light_data <- all_info_split %>%
  dplyr::filter(chains %in% c("IGK", "IGL"))

light_data$charge <- alakazam::charge(seq = light_data$cdr3)

# Barplot, y = charge, x = tetramer, color = status
l_charge <- ggplot2::ggplot(light_data, ggplot2::aes(y = charge,
                                                     x = tet_name_cutoff,
                                                     fill = Status)) +
  ggplot2::geom_boxplot(size = 0.1, outlier.size = 0.5) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1)) +
  ggplot2::ggtitle("Heavy chain SMH")


pdf(file.path(image_dir, "5C_light_charge.pdf"),
    width = 8, height = 4)

print(l_charge)

dev.off()
graphics.off()

### 5D -------------------------------------------------------------------------

# JH gene usage by status

# First find the percent J gene usage by status
j_data <- heavy_data[ , c("j_gene", "sample", "Status")]

j_data <- j_data %>%
  dplyr::group_by(sample) %>% 
  dplyr::add_count(name = "sample_count") %>% # Total cells in sample
  dplyr::ungroup() %>%
  dplyr::group_by(sample, j_gene) %>%
  dplyr::add_count(name = "j_count") %>% # Total cells per j gene per sample
  dplyr::distinct() %>% # One line per sample/j gene
  dplyr::mutate(percent_j = j_count / sample_count * 100)


# Do t test
stat.test <- j_data %>%
  group_by(j_gene) %>%
  t_test(percent_j ~ Status) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "Status")


# Make a plot
j_plot <- ggplot2::ggplot(j_data, ggplot2::aes(x = Status,
                                               y = percent_j)) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = Status)) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::facet_grid(~j_gene, switch = "x") +
  ggpubr::stat_pvalue_manual(
    stat.test, 
    label = "p.adj.signif"
  )

pdf(file.path(image_dir, "5D_jh_usage.pdf"),
    width = 6, height = 6)

print(j_plot)

dev.off()
graphics.off()

### 5E -------------------------------------------------------------------------

# 5E SHM VH and VL

density_h_smh <- ggplot2::ggplot(heavy_data, ggplot2::aes(y = Status,
                                                          x = all_mis_freq,
                                                          fill = Status)) +
  ggridges::geom_density_ridges() +
  ggplot2::facet_grid(~tet_name_cutoff, switch = "x") +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1)) +
  ggplot2::ggtitle("Heavy chain SMH") +
  ggplot2::xlim(-0.02, 0.2) +
  ggplot2::coord_flip()

boxplot_h_smh <- ggplot2::ggplot(heavy_data, ggplot2::aes(y = (all_mis_freq),
                                                          x = tet_name_cutoff,
                                                          fill = Status)) +
  ggplot2::geom_boxplot(size = 0.1, outlier.size = 0.5) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1)) +
  ggplot2::ggtitle("Heavy chain SMH")

violin_h_smh <- ggplot2::ggplot(heavy_data, ggplot2::aes(y = (all_mis_freq),
                                                          x = tet_name_cutoff,
                                                          fill = Status)) +
  ggplot2::geom_violin(scale = "width", position = ggplot2::position_dodge(1)) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1)) +
  ggplot2::ggtitle("Heavy chain SMH") +
  ggplot2::stat_summary(fun = median, geom = "point", size = 2,
                        position = ggplot2::position_dodge(1)) +
  ggplot2::ylim(c(-0.01, 0.2))

light_data <- all_info_split %>%
  dplyr::filter(chains %in% c("IGK", "IGL"))

boxplot_l_smh <- ggplot2::ggplot(light_data, ggplot2::aes(y = all_mis_freq,
                                                          x = tet_name_cutoff,
                                                          fill = Status)) +
  ggplot2::geom_boxplot(size = 0.1, outlier.size = 0.5) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1)) +
  ggplot2::ggtitle("Light chain SMH")

violin_l_smh <- ggplot2::ggplot(light_data, ggplot2::aes(y = all_mis_freq,
                                                          x = tet_name_cutoff,
                                                          fill = Status)) +
  ggplot2::geom_violin(scale = "width", position = ggplot2::position_dodge(1)) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1)) +
  ggplot2::ggtitle("Light chain SMH") +
  ggplot2::stat_summary(fun = median, geom = "point", size = 2,
                        position = ggplot2::position_dodge(1)) +
  ggplot2::ylim(c(-0.01, 0.2))

density_l_smh <- ggplot2::ggplot(light_data, ggplot2::aes(y = Status,
                                                          x = all_mis_freq,
                                                          fill = Status)) +
  ggridges::geom_density_ridges() +
  ggplot2::facet_grid(~tet_name_cutoff, switch = "x") +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angl = 45, hjust = 1)) +
  ggplot2::ggtitle("Light chain SMH") +
  ggplot2::xlim(-0.02, 0.2) +
  ggplot2::coord_flip()


h_and_l <- cowplot::plot_grid(boxplot_h_smh, boxplot_l_smh,
                              nrow = 2, ncol = 1)

h_and_l_violin <- cowplot::plot_grid(violin_h_smh, violin_l_smh,
                                     nrow = 2, ncol = 1)

h_and_l_density <- cowplot::plot_grid(density_h_smh, density_l_smh,
                                      nrow = 2, ncol = 1)

pdf(file.path(image_dir, "5E_heavy_ligh_smh.pdf"),
    width = 8, height = 8)

print(h_and_l)

dev.off()
graphics.off()

pdf(file.path(image_dir, "5E_heavy_ligh_smh_density.pdf"),
    width = 8, height = 8)

print(h_and_l_density)

dev.off()
graphics.off()

pdf(file.path(image_dir, "5E_heavy_ligh_smh_violin.pdf"),
    width = 8, height = 8)

print(h_and_l_violin)

dev.off()
graphics.off()

## Figure 6 --------------------------------------------------------------------
### 6A -------------------------------------------------------------------------

# Shannon index
seurat_data$sample_celltype <- paste(seurat_data$sample,
                                     seurat_data$final_celltype,
                                     sep = "_")
seurat_data <- calc_diversity(seurat_data,
                              data_col    = "final_clone",
                              cluster_col = "sample_celltype",
                              method      = abdiv::shannon
)

diversity <- seurat_data[[]] %>%
  dplyr::select(final_clone_shannon_diversity, final_celltype, sample,
                Status) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(final_clone_shannon_diversity))

diversity_p <- ggplot2::ggplot(diversity, ggplot2::aes(x = Status,
                                                       y = final_clone_shannon_diversity,
                                                       fill = Status)) + 
  ggplot2::geom_boxplot(size = 0.25, outlier.size = 0.5) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::facet_grid(~final_celltype, switch = "x") +
  ggpubr::stat_compare_means(ggplot2::aes(group = Status),
                             method = "t.test",
                             label = "p.signif",
                             comparisons = list(c("AAB", "ND"),
                                                c("T1D", "ND"),
                                                c("T1D", "AAB"))) +
  ggplot2::theme(strip.placement = "outside",
                 strip.background = element_rect(fill = NA, color = "white"),
                 panel.spacing = unit(-.01,"cm"),
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

pdf(file.path(image_dir, "6A_shannon_alpha_diversity.pdf"),
    width = 8, height = 8)

print(diversity_p)

dev.off()
graphics.off()

### 6B -------------------------------------------------------------------------

# Tree maps
all_data <- all_info_split[all_info_split$chains == "IGH", c("v_gene", "Status")]

all_plots <- lapply(levels(all_data$Status), function(x){
  df <- all_data %>%
    dplyr::filter(Status == x) %>%
    dplyr::group_by(v_gene) %>%
    dplyr::add_count(name = "v_count") %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(v_percent = v_count / sum(v_count) * 100) %>%
    dplyr::mutate(plot_v = ifelse(v_percent > 0.25, v_gene, ""))
  
  
  # Fill could be v_percent
  plot <- ggplot(df, aes(area = v_percent, fill = v_gene, label = plot_v)) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_text(colour = "white", place = "center", reflow = T) +
    # viridis::scale_fill_viridis(option = "inferno",
    #                             discrete = FALSE, n.breaks = 50) +
    ggplot2::scale_fill_manual(values = heavy_colors) +
    ggplot2::ggtitle(x) +
    ggplot2::theme(legend.position = "none", 
                   plot.title = ggplot2::element_text(hjust = 0.5))
  
  return(plot)
})

combined_plot <- cowplot::plot_grid(plotlist = all_plots, nrow = 1,
                                    ncol = 3)


pdf(file.path(image_dir, "6B_vdj_percent_all.pdf"),
    width = 16, height = 6)

print(combined_plot)

dev.off()

### 6C -------------------------------------------------------------------------


# Odds ratio
stats_aab <- file.path(save_dir, "files", "vdj_files",
                       "stats", "all_aab_nd.xlsx")

stats_t1d <- file.path(save_dir, "files", "vdj_files",
                       "stats", "all_t1d_nd.xlsx")

stats_aab_light <- file.path(save_dir, "files", "vdj_files",
                             "stats", "all_light_aab_nd.xlsx")

stats_t1d_light <- file.path(save_dir, "files", "vdj_files",
                             "stats", "all_light_t1d_nd.xlsx")

stats_aab_hl <- file.path(save_dir, "files", "vdj_files",
                          "stats", "all_heavy_light_aab_nd.xlsx")

stats_t1d_hl <- file.path(save_dir, "files", "vdj_files",
                          "stats", "all_heavy_light_t1d_nd.xlsx")


# Heavy
odds_ratio_aab <- openxlsx::readWorkbook(stats_aab, sheet = "odds_ratio")

odds_ratio_t1d <- openxlsx::readWorkbook(stats_t1d, sheet = "odds_ratio")

odds_ratio_t1d <- odds_ratio_t1d[order(odds_ratio_t1d$odds_ratio),]

odds_ratio_aab <- odds_ratio_aab[order(match(odds_ratio_aab$v_gene,
                                             odds_ratio_t1d$v_gene)),]

# Light
odds_ratio_aab_light <- openxlsx::readWorkbook(stats_aab_light,
                                               sheet = "odds_ratio")

odds_ratio_t1d_light <- openxlsx::readWorkbook(stats_t1d_light,
                                               sheet = "odds_ratio")

odds_ratio_t1d_light <- odds_ratio_t1d_light[
  order(odds_ratio_t1d_light$odds_ratio),]

odds_ratio_aab_light <- odds_ratio_aab_light[
  order(match(odds_ratio_aab_light$v_gene,
              odds_ratio_t1d_light$v_gene)),]

# Heavy and light
odds_ratio_aab_hl <- openxlsx::readWorkbook(stats_aab_hl,
                                            sheet = "odds_ratio")

odds_ratio_t1d_hl <- openxlsx::readWorkbook(stats_t1d_hl,
                                            sheet = "odds_ratio")

odds_ratio_t1d_hl <- odds_ratio_t1d_hl[
  order(odds_ratio_t1d_hl$odds_ratio),]

odds_ratio_t1d_hl <- odds_ratio_t1d_hl[
  order(match(odds_ratio_t1d_hl$v_gene,
              odds_ratio_t1d_hl$v_gene)),]


make_plot <- function(odds_ratio, status, chain){
  # Set the order
  odds_ratio$v_gene <- factor(odds_ratio$v_gene,
                              levels = odds_ratio$v_gene)
  col_mapping <- c("P>0.05" = "#D3D3D3",
                   "P<0.05" = as.character(status_colors[status]))
  
  odds_plot <- ggplot2::ggplot(odds_ratio, 
                               ggplot2::aes(x = odds_ratio,
                                            y = v_gene,
                                            color = padj_category)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = X95_conf_int_low,
                                        xmax = X95_conf_int_high)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
    ggplot2::scale_color_manual(values = col_mapping) +
    ggplot2::ggtitle(paste0(chain, " chain: ", status, " vs ND")) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5))
}

t1d <- make_plot(odds_ratio_t1d, status = "T1D", chain = "Heavy")
aab <- make_plot(odds_ratio_aab, status = "AAB", chain = "Heavy")

# Remove chains with a very high confidence interval
odds_ratio_t1d_light <- odds_ratio_t1d_light[odds_ratio_t1d_light$X95_conf_int_high < 10,]
odds_ratio_aab_light <- odds_ratio_aab_light[odds_ratio_aab_light$X95_conf_int_high < 10,]

t1d_light <- make_plot(odds_ratio_t1d_light, status = "T1D", chain = "Light")
aab_light <- make_plot(odds_ratio_aab_light, status = "AAB", chain = "Light")

full_plot <- cowplot::plot_grid(t1d, aab, t1d_light, aab_light,
                                nrow = 2, ncol = 2, rel_heights = c(1, 1.2))

heavy_plot <- cowplot::plot_grid(t1d, aab,
                                 nrow = 1, ncol = 2)

light_plot <- cowplot::plot_grid(t1d_light, aab_light,
                                 nrow = 1, ncol = 2)

pdf(file.path(image_dir, "6C_all_odds.pdf"),
    width = 15, height = 20)

print(full_plot)

dev.off()
graphics.off()


pdf(file.path(image_dir, "6C_heavy_odds.pdf"),
    width = 15, height = 10)

print(heavy_plot)

dev.off()
graphics.off()

pdf(file.path(image_dir, "6C_heavy_odds.pdf"),
    width = 15, height = 11)

print(light_plot)

dev.off()
graphics.off()

### 6D -------------------------------------------------------------------------

# Tree maps light
all_data <- all_info_split[all_info_split$chains %in% c("IGK", "IGL"), c("v_gene", "Status")]

all_plots <- lapply(levels(all_data$Status), function(x){
  df <- all_data %>%
    dplyr::filter(Status == x) %>%
    dplyr::group_by(v_gene) %>%
    dplyr::add_count(name = "v_count") %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(v_percent = v_count / sum(v_count) * 100) %>%
    dplyr::mutate(plot_v = ifelse(v_percent > 0.25, v_gene, ""))
  
  
  # Fill could be v_percent
  plot <- ggplot(df, aes(area = v_percent, fill = v_gene, label = plot_v)) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_text(colour = "white", place = "center", reflow = T) +
    # viridis::scale_fill_viridis(option = "inferno",
    #                             discrete = FALSE, n.breaks = 50) +
    ggplot2::scale_fill_manual(values = light_colors) +
    ggplot2::ggtitle(x) +
    ggplot2::theme(legend.position = "none", 
                   plot.title = ggplot2::element_text(hjust = 0.5))
  
  return(plot)
})

combined_plot <- cowplot::plot_grid(plotlist = all_plots, nrow = 1,
                                    ncol = 3)


pdf(file.path(image_dir, "6D_vdj_percent_light.pdf"),
    width = 16, height = 6)

print(combined_plot)

dev.off()


### 6F -------------------------------------------------------------------------

# Only odds increased heavy light pairs
# Get only significant chains
sig_hl_t1d <- odds_ratio_t1d_hl[odds_ratio_t1d_hl$p_adj < 0.05 &
                                  odds_ratio_t1d_hl$odds_ratio > 1,]
sig_hl_aab <- odds_ratio_aab_hl[odds_ratio_aab_hl$p_adj < 0.05 &
                                  odds_ratio_aab_hl$odds_ratio > 1,]


all_chains <- c(sig_hl_t1d$v_gene, sig_hl_aab$v_gene)
use_cells <- seurat_data[[]] %>%
  dplyr::filter(v_gene %in% all_chains)

use_data <- all_info_split[all_info_split$barcode %in%
                             rownames(use_cells), ]

all_data <- count_genes_heavy_light(starting_df = use_data,
                                    group_by = "Status",
                                    subset_counts = 0,
                                    color_list = status_colors)


pdf(file.path(image_dir, "6F_combined_heavy_light_significant.pdf"),
    width = 8, height = 8)

par(mfrow = c(1, 1))
make_circos_plot(circos_df = all_data$df,
                 color = all_data$color_list, 
                 grid_color = all_colors)

dev.off()
graphics.off()

### 6G -------------------------------------------------------------------------


all_data <- count_genes_heavy_light(starting_df = use_data,
                                    group_by = "tet_name_cutoff",
                                    subset_counts = 0,
                                    color_list = tetramer_colors)


pdf(file.path(image_dir, "6G_combined_heavy_light_significant_tet_color.pdf"),
    width = 8, height = 8)

par(mfrow = c(1, 1))
make_circos_plot(circos_df = all_data$df,
                 color = all_data$color_list, 
                 grid_color = all_colors)

dev.off()
graphics.off()

## Figure 7 --------------------------------------------------------------------
clone_data <- file.path(save_dir, "files", "v_gene_counting", 
                        "clone_expansion.xlsx")

clone_expanded <- openxlsx::readWorkbook(clone_data, sheet = "expanded_clones")

### 7A -------------------------------------------------------------------------

use_cells <- clone_expanded$barcode
use_data <- all_info_split[all_info_split$barcode %in%
                             use_cells, ]

use_data$tet_name_cutoff <- factor(use_data$tet_name_cutoff,
                                   levels = plotting_levels)

pdf(file.path(image_dir, "7A_expanded_clone_circos.pdf"),
    width = 32, height = 16)

par(mfrow = c(2, 4))
for(tet_group in levels(use_data$tet_name_cutoff)){
  sub_data <- use_data[use_data$tet_name_cutoff == tet_group, ]
  all_data <- count_genes_heavy_light(starting_df = sub_data,
                                      group_by = "Status",
                                      subset_counts = 0,
                                      color_list = status_colors)  
  
  make_circos_plot(circos_df = all_data$df,
                   color = all_data$color_list, 
                   grid_color = all_colors)
  
  title(main = tet_group, cex.main = 1.5)
}

dev.off()

### 7B -------------------------------------------------------------------------

# Clone plots
clone_expanded <- openxlsx::readWorkbook(clone_data, sheet = "expanded_clones")

clone_expanded$tet_name_cutoff <- ifelse(clone_expanded$tet_name_cutoff == "Other_Multi_Reactive" & 
                                           grepl("INS|GAD|IA2", clone_expanded$full_tet_name_cutoff),
                                         "Islet_Multi_Reactive", clone_expanded$tet_name_cutoff)


make_clone_barplot <- function(clone_df){
  clone_df <- clone_df %>%
    dplyr::group_by(Status, tet_name_cutoff) %>%
    dplyr::add_count(name = "number_of_clones") %>%
    dplyr::select(tet_name_cutoff, number_of_clones, Status) %>%
    dplyr::distinct()
  
  clone_df$Status <- factor(clone_df$Status,
                            levels = status_levels)
  
  clone_df$tet_name_cutoff <- factor(clone_df$tet_name_cutoff,
                                     levels = plotting_levels)
  
  
  tet_colors <- tetramer_colors
  names(tet_colors) <- make.names(names(tet_colors))
  
  clone_barplot <- ggplot2::ggplot(data = clone_df, 
                                   ggplot2::aes(x = Status, y = number_of_clones,
                                                fill = tet_name_cutoff)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::scale_fill_manual(values = tet_colors)
  
  return(clone_barplot)
  
}

expanded_barplot <- make_clone_barplot(clone_expanded)

pdf(file.path(image_dir, "7B_private_clones.pdf"),
    width = 8, height = 8)

print(expanded_barplot)

dev.off()
graphics.off()

### 7C -------------------------------------------------------------------------

# Public clone circos
clone_public <- openxlsx::readWorkbook(clone_data, sheet = "public_clones")


use_cells <- clone_public$barcode

use_data <- all_info_split[all_info_split$barcode %in%
                             use_cells, ]

use_data <- use_data[order(use_data$sample),]

all_data <- count_genes_heavy_light(starting_df = use_data,
                                    group_by = "sample",
                                    subset_counts = 0,
                                    color_list = sample_colors)


pdf(file.path(image_dir, "7C_public_circos_sample_color.pdf"),
    width = 8, height = 8)

par(mfrow = c(1, 1))
make_circos_plot(circos_df = all_data$df,
                 color = all_data$color_list, 
                 grid_color = all_colors)

dev.off()
graphics.off()

my_hist <- ggplot(use_data, aes(x = cdr3_length, fill = sample)) + 
  geom_bar() +
  ggplot2::scale_fill_manual(values = sample_colors)

legend <- cowplot::get_legend(my_hist)
pdf(file.path(image_dir, "7C_legened.pdf"),
    width = 8, height = 8)
grid.draw(legend)
dev.off()

### 7D -------------------------------------------------------------------------

all_data <- count_genes_heavy_light(starting_df = use_data,
                                    group_by = "tet_name_cutoff",
                                    subset_counts = 0,
                                    color_list = tetramer_colors)


pdf(file.path(image_dir, "7D_public_circos_tet_color.pdf"),
    width = 8, height = 8)

par(mfrow = c(1, 1))
make_circos_plot(circos_df = all_data$df,
                 color = all_data$color_list, 
                 grid_color = all_colors)

dev.off()
graphics.off()

### 7E -------------------------------------------------------------------------

# Public clones

# First just select out number of samples, Status, status percent
plot_public <- clone_public[, c("final_clone", "number_of_samples",
                                "Status", "status_percent", "v_gene",
                                "j_gene", "clone_id")]

# Remove duplicate rows
plot_public <- unique(plot_public)

# Set status levels
plot_public$Status <- factor(plot_public$Status, levels = c("T1D", "AAB", "ND"))

# Order by number of samples and then by status
ordered_samples <- unique(plot_public$number_of_samples)
ordered_samples <- ordered_samples[order(ordered_samples)]
all_samples <- lapply(ordered_samples, function(number){
  sample_df <- plot_public[plot_public$number_of_samples == number,]
  
  # Now add ordering based on the status
  all_clones <- lapply(unique(sample_df$final_clone), function(clone){
    clone_df <- sample_df[sample_df$final_clone == clone, ]
    clone_df <- clone_df[order(clone_df$Status),]
    all_status <- paste(unique(clone_df$Status), collapse = "_")
    
    # This is an arbitrary order based on what status is included
    status_mapping <- c("T1D" = 1,
                        "T1D_AAB" = 2,
                        "AAB" = 3,
                        "T1D_ND" = 4,
                        "AAB_ND" = 5,
                        "T1D_AAB_ND" = 6,
                        "ND" = 7)
    
    clone_df$order_num <- status_mapping[all_status]
    
    return(clone_df)
  })
  all_clones <- do.call(rbind, all_clones)
  all_clones <- all_clones[order(all_clones$order_num),]
  return(all_clones)
})


all_samples <- do.call(rbind, all_samples)

plot_public$final_clone <- factor(plot_public$final_clone,
                               levels = unique(all_samples$final_clone))


# The status percent when added together must equal the number of samples
plot_public$status_fill <- (plot_public$status_percent / 100) * plot_public$number_of_samples


plot_public_test <- plot_public %>%
  dplyr::group_by(final_clone) %>%
  dplyr::mutate(test = sum(status_fill))

bar_plot <- ggplot2::ggplot(plot_public, ggplot2::aes(x = final_clone, y = status_fill,
                                          fill = Status)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::ylab("Number of samples") +
  ggplot2::xlab("Clone id")
  
pdf(file.path(save_dir, "images", "final_figures",
              "7E_public_clone_barplot.pdf"), 
    height = 4, width = 12) 

print(bar_plot)

dev.off()

## Supplement ------------------------------------------------------------------

### Supp 1 ---------------------------------------------------------------------
# Violin plots of mito, features, genes across samples
quality_plot <- featDistPlot(seurat_data, geneset = c("nCount_RNA",
                                                      "nFeature_RNA", 
                                                      "percent.mt"),
                             sep_by = "sample", col_by = "Status",
                             color = status_colors,
                             combine = FALSE)

quality_plot$percent.mt <- quality_plot$percent.mt +
  ggplot2::ylim(0, 20)

final_quality <- cowplot::plot_grid(plotlist = quality_plot, nrow = 3)

pdf(file.path(image_dir, "Supp1_quality_plots.pdf"))
print(final_quality)
dev.off()

### Supp 2 ---------------------------------------------------------------------
# UMAP of samples before and after correction

umap_one <- plotDimRed(seurat_data, col_by = "sample",
                       plot_type = "pca.umap", color = sample_colors,
                       ggrastr = TRUE,
                       reorder_cells = TRUE)[[1]]

umap_two <- plotDimRed(seurat_data, col_by = "final_celltype",
                       plot_type = "pca.umap", color = final_colors,
                       ggrastr = TRUE, reorder_cells = TRUE)[[1]]


umap_three <- plotDimRed(seurat_data, col_by = "sample",
                         plot_type = "rna_mnn.umap", color = sample_colors,
                         ggrastr = TRUE, reorder_cells = TRUE)[[1]]

umap_four <- plotDimRed(seurat_data, col_by = "final_celltype",
                        plot_type = "rna_mnn.umap", color = final_colors,
                        ggrastr = TRUE, reorder_cells = TRUE)[[1]]

full_plot <- cowplot::plot_grid(umap_one, umap_two, umap_three, umap_four,
                                rel_widths = c(0.84, 1, 0.84, 1),
                                nrow = 2, ncol = 2)

pdf(file.path(image_dir, "Supp2_quality_umaps.pdf"),
    width = 12, height = 10)

print(full_plot)

dev.off()

### Supp 3 ---------------------------------------------------------------------
# Insr mRNA by Ag reactivity and Status

save_violin <- featDistPlot(seurat_data, geneset = "INSR",
                            sep_by = "tet_name_cutoff",
                            col_by = "Status", color = status_colors,
                            combine = FALSE)

pdf(file.path(image_dir, "Supp3_insr_expr.pdf"),
    height = 4, width = 8)
print(save_violin)
dev.off()

graphics.off()

### Supp 6 ---------------------------------------------------------------------
pdf(file.path(image_dir, "Supp6_tetramer_new_vs_old_cutoff.pdf"))

cm2 <- confusionMatrix(seurat_data$tet_name_cutoff, seurat_data$scar_libra_tet_hash_id)
cm2 <- cm2 / rowSums(cm2)
print(pheatmap::pheatmap(cm2, main = "previous_vs_libra_id"))

dev.off()

graphics.off()

### Supp 7 ---------------------------------------------------------------------
#	List of top 10 genes per cell type
Idents(seurat_data) <- seurat_data$celltype_cluster

all_markers <- FindAllMarkers(seurat_data, assay = "RNA", 
                              logfc.threshold = 0.25, max.cells.per.ident = 2000,
                              only.pos = TRUE)

save_excel <- openxlsx::createWorkbook()

for(clust in unique(all_markers$cluster)){
  write_markers <- all_markers %>%
    dplyr::filter(cluster == clust, p_val_adj < 0.05)
  
  openxlsx::addWorksheet(wb = save_excel, sheetName = clust)
  
  openxlsx::writeData(wb = save_excel, sheet = clust, x = write_markers)
}

openxlsx::saveWorkbook(wb = save_excel, 
                       file = file.path(image_dir, "cluster_markers.xlsx"),
                       overwrite = TRUE)

# Ran this to make sure the results didn't change with the downsampling,
# they did not.
# all_markers_short <- all_markers %>%
#   dplyr::filter(abs(avg_log2FC) > 1)

use_markers <- all_markers %>%
  dplyr::filter(!grepl("IGK|IGH|IGL|AL139020.1", gene)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 10, wt = avg_log2FC) %>%
  dplyr::arrange(cluster, desc(avg_log2FC))

graphics.off()

pdf(file.path(image_dir, "Supp7_marker_heatmap_average.pdf"),
    width = 8, height = 12)

print(plot_heatmap(seurat_data, gene_list = use_markers$gene,
                   colors = cluster_celltype_colors, meta_col = "celltype_cluster",
                   average_expression = TRUE, assay = "RNA"))

dev.off()


graphics.off()

### Supp clones ----------------------------------------------------------------