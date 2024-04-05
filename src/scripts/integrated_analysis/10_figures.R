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

# Colors -----------------------------------------------------------------------
all_colors <- readRDS(file = file.path("files/all_colors.rds"))


final_colors <- all_colors$cell_type_colors

final_colors <- final_colors[names(final_colors) != "Translational_Intermediate"]

tetramer_colors <- all_colors$tetramer_colors
names(tetramer_colors) <- make.names(names(tetramer_colors))

sample_colors <- all_colors$sample_colors


status_colors <- all_colors$status_colors

# Main figures -----------------------------------------------------------------

## Figure 1 --------------------------------------------------------------------

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
all_combinations <- combn(unique(new_status), m = 2)

all_p_val <- lapply(1:ncol(all_combinations), function(col){
  val1 <- all_combinations[1, col]
  val2 <- all_combinations[2, col]
  t_test_val <- t.test(facs_percents[facs_percents$Status == val1,]$percent,
                       facs_percents[facs_percents$Status == val2,]$percent)
  return_df <- data.frame("group1" = val1, "group2" = val2, 
                          "p_value" = t_test_val$p.value)
  return(return_df)
})

all_p_val <- do.call(rbind, all_p_val)

# Make boxplot
# Add in stats
p1 <- ggplot2::ggplot(facs_percents, ggplot2::aes(x = Status, y = percent)) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = Status)) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.75)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggpubr::stat_pvalue_manual(
    all_p_val, 
    y.position = 12, step.increase = 0.1,
    label = "p_value"
  )

pdf(file.path(image_dir, "1C_ag_positive_facs.pdf"),
    width = 6, height = 4)
print(p1)
dev.off()


## Figure 2 --------------------------------------------------------------------
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

# Make a UMAP of all cells colored by cell type
pdf(file.path(image_dir, "2A_rna_umap.pdf"), width = 8, height = 8)

print(plotDimRed(seurat_data, col_by = "celltype_cluster", 
                 plot_type = "rna_mnn.umap",
                 color = cluster_celltype_colors, ggrastr = TRUE))

dev.off()

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

# Cell type by status
pdf(file.path(image_dir, "2C_status_celltype_barplot.pdf"),
    width = 8, height = 8)
barplot <- stacked_barplots(seurat_data, meta_col = "celltype_cluster",
                            split_by = "Status",
                            color = cluster_celltype_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

print(barplot)

dev.off()

## Figure 3 --------------------------------------------------------------------

# Barplot of antigen binding
pdf(file.path(image_dir, "3A_status_antigen_barplot.pdf"),
    width = 8, height = 8)
barplot <- stacked_barplots(seurat_data, meta_col = "tet_name_cutoff",
                            split_by = "Status",
                            color = tetramer_colors) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

print(barplot)

dev.off()

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

seurat_islet <- subset(seurat_data, 
                       subset = tet_name_cutoff == "Islet_Multi_Reactive")

fine_tet_options <- unique(seurat_islet$full_tet_name_cutoff)

fine_tet_options <- fine_tet_options[order(fine_tet_options)]

fine_tet_colors <- grDevices::colorRampPalette(
  colors = RColorBrewer::brewer.pal(n = 9, name = "Set1"))(length(fine_tet_options))

names(fine_tet_colors) <- fine_tet_options

barplot2 <- make_barplot(seurat_islet, tetramer_colors_use = fine_tet_colors,
                         meta_data_col = "full_tet_name_cutoff")

pdf(file.path(image_dir, "3D_islet_multi_tetramer_stacked_barplot.pdf"),
    height = 8, width = 12)
print(barplot2)
dev.off()


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
twin_markers <- read.csv(file.path(save_dir, "files", "twin_mast_de.csv"))

# Make heatmap of all DE
de_genes <- markers_sig[markers_sig$cluster !=
                          "T1D_Islet_Reactive_AAB_Islet_Reactive",]$gene

de_genes <- unique(de_genes)

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

# Use jaccard distance to make a dendrogram of gse terms
gse_res <- read.csv(file.path(save_dir, "files/de/all_GSE_results.csv"))

name_mapping <- c("T1D" = "T1D_Islet_Reactive_ND_Islet_Reactive",
                  "AAB" = "AAB_Islet_Reactive_ND_Islet_Reactive")

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


#' Generates an hclust object based on the comparison (query) and a status.
#' 
#' This assumes that gene ontology was run using `gprofiler2::gost`. The
#' final output will be a dotplot where the size is the number of
#' intersectig genes and the color is the log of the p value. These
#' dots will be on the end of the dendrogram built on the similarity 
#' between kett terms. Currently this function only works for KEGG because
#' the genes in each term need to be downloaded and only KEGG is currently
#' supported
#' 
#' @param query The query from the gse_res, the gse_res will be subset
#' by this query
#' @param gse_res the results table returned by `gprofiler2::gost`
#' @param p_val_cutoff what p_val to use to subset terms
make_histogram <- function(query, gse_res, p_val_cutoff = 0.01){
  
  # Subset to only kegg terms and the query of interest
  gse_res <- gse_res[gse_res$source == "KEGG" & 
                       gse_res$query == query,]
  
  # Make a list of all terms
  gse_list <- gse_res$term_id
  names(gse_list) <- gse_res$term_name
  
  # Find all genes associated with the kegg term
  # This could probably go faster if all genes from all kegg terms are
  # first pulled and keggConv and mapIds are only run once.
  all_genes <- lapply(gse_list, function(x){
    genes <- get_gene_list(x)
    
    return(genes)
  })
  
  # Make a matrix of the correct size, it should be all gene sets by all gene
  # sets
  total_values <- length(all_genes) ^ 2
  final_jaccard_matrix <- Matrix::Matrix(data = rep(0, total_values),
                                         nrow = length(all_genes),
                                         ncol = length(all_genes),
                                         sparse = TRUE)
  
  # Find the jaccard distance for all pairwise gene sets
  for (x in 1:length(all_genes)){
    for (y in 1:length(all_genes)){
      
      save_val <- jaccard(a = all_genes[[x]],
                          b = all_genes[[y]])
      
      # Save to both relevant positions in the matrix so you don't
      # run twice
      final_jaccard_matrix[x, y] <- save_val
      final_jaccard_matrix[y, x] <- save_val
    }
  }
  
  colnames(final_jaccard_matrix) <- names(all_genes)
  rownames(final_jaccard_matrix) <- names(all_genes)
  
  # Make into a distance matrix and run hclust
  hc_full <- stats::hclust(stats::as.dist(1- final_jaccard_matrix),
                           method = "ward.D2")
  
  # Subset based on the desired p-value
  gse_small <- gse_res[gse_res$p_value < p_val_cutoff,]
  
  jaccard_small <- final_jaccard_matrix[rownames(final_jaccard_matrix)
                                        %in% gse_small$term_name,
                                        colnames(final_jaccard_matrix)
                                        %in% gse_small$term_name]
  
  # Repeat clustering on the smaller data set
  hc_small <- stats::hclust(stats::as.dist(1- jaccard_small),
                            method = "ward.D2")
  
  return(list(full = hc_full, small = hc_small, gse_res = gse_res))
  
}

#' Finds all gene IDs from a a kegg term
#' 
#' This takes a kegg term and a path of all kegg terms and maps
#' to ncbi gene ids
#' 
#' @param kegg_term The list of kegg terms to query. If multiple are provided
#' a list of genes from all terms combined will be returned.
get_gene_list <- function(kegg_term){
  path <- keggLink("pathway", "hsa")
  
  search_term <- gsub("KEGG:", "", kegg_term)
  
  tryCatch({
    # Find genes in term
    all_genes <- lapply(search_term, function(term){
      genes_in_term <- path[grepl(term, path)]
      return(genes_in_term)
    })
    all_genes <- unlist(all_genes)
    
    # Translate to gene id
    all_gene_id <- keggConv("ncbi-geneid", unique(names(all_genes)))
    
    all_gene_id <- gsub("ncbi-geneid:", "", all_gene_id)
    
    # Map to gene symbols
    gene_symbols <- mapIds(org.Hs.eg.db, keys = all_gene_id, 
                           keytype = "ENTREZID", column = "SYMBOL")  
    
  }, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
    gene_symbols <- "NA"
  })
  
  return(as.character(gene_symbols))
}

#' Makes a histogram based on the output from `make_histogram`.
#' 
#' This assumes that `make_histogram` has already been run. Will take
#' that output and make a dendrogram
#' 
#' @param hc The hclust object output from `make_h`
#' @param gse_res the results table returned by `gprofiler2::gost`
#' @param y_lim how to set the y-lim, this is important to modify if
#' you save the figure in different sizes, it moves where the key is 
#' relative to the text
#' @param title The title of the plot
plot_histogram <- function(hc, gse_res, ylim = c(2, -2), title = NULL){
  # Convert hclust to dendrogram
  dend <- as.dendrogram(hc)
  
  # Make dendrogram into an object that plays well with ggplot
  ddata <- ggdendro::dendro_data(dend, type = "rectangle")
  
  # Add gse information
  dot_data <- merge(ddata$labels, gse_res, by.x = "label", by.y = "term_name",
                    all.x = TRUE, all.y = FALSE)
  
  
  # Make the plot with the dendrogram passed to geom segment and the 
  # gse res as the geom_point
  p <- ggplot(ggdendro::segment(ddata)) + 
    ggplot2::geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    ggplot2::geom_text(data = ddata$labels, ggplot2::aes(x = x, y = y - 0.06,
                                                         label = label),
                       hjust = 0, angle = 0, size = 3) +
    ggplot2::geom_point(data = dot_data, ggplot2::aes(x = x, y = y,
                                                      color = -log(p_value),
                                                      size = intersection_size)) +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    ggplot2::coord_flip() + 
    ggplot2::scale_y_reverse(expand = c(0.2, 0)) +
    ggplot2::ylim(ylim) +
    ggplot2::theme(axis.line = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title = element_blank())
  
  if(!is.null(title)){
    p <- p +
      ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  
  return(p)
  
}

#' Makes a volcano
#' 
#' Makes a volcano plot based on a previously run DE analysis. Assumes 
#' columns named `p_val_adj`, `avg_log2FC` and `gene` are present. Also
#' takes a list of kegg pathways. Any gene in any of those pathways will be
#' colored on top of the volcano plot. Logic to fix the p-values was borrowed
#' from `EnhancedVolcano`
#' 
#' @param de_genes The de genes to plot
#' @param kegg_pathways A list of kegg pathways. Genes from these will
#' be highlighted
#' @param group_name Name for the kegg pathways. This will be the name
#' of the plot
#' @param color What color to color the kegg pathway genes
make_volcano <- function(de_genes, kegg_pathways, group_name,
                         color = NULL){
  # Pull out all genes to highlight on the volcano plot
  highlight_genes <- get_gene_list(kegg_pathways)
  
  if(is.null(color)){
    color <- "#C41E3A"
  }
  
  # Fix y values (taken from enhanced volcano)
  if (min(de_genes$p_val_adj, na.rm=TRUE) == 0) {
    # <= v1.2
    #warning(paste("One or more P values is 0.",
    #  "Converting to minimum possible value..."),
    #  call. = FALSE)
    #toptable[which(toptable[[y]] == 0), y] <- .Machine$double.xmin
    warning(paste('One or more p-values is 0.',
                  'Converting to 10^-1 * current',
                  'lowest non-zero p-value...'),
            call. = FALSE)
    
    de_genes[which(de_genes[["p_val_adj"]] == 0), "p_val_adj"] <- min(
      de_genes[which(de_genes[["p_val_adj"]] != 0), "p_val_adj"],
      na.rm = TRUE) * 10^-1
  }
  
  # Add in colors
  de_genes$color <- ifelse(de_genes$gene %in% highlight_genes,
                           "kegg_gene", "other_gene")
  
  de_genes$color <- factor(de_genes$color, levels = c("other_gene",
                                                      "kegg_gene"))
  
  de_genes <- de_genes[order(de_genes$color),]
  
  color_mapping <- c("kegg_gene" = color, "other_gene" = "#999999")
  
  
  # Make a volcano plot
  volcano_plot <- ggplot2::ggplot(de_genes, ggplot2::aes(x = avg_log2FC,
                                                         y = -log(p_val_adj),
                                                         color = color)) +
    ggplot2::geom_point() + 
    ggplot2::xlim(-2, 2) +
    ggplot2::scale_color_manual(values = color_mapping) +
    ggplot2::ggtitle(group_name) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  return(volcano_plot)
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

sample_levels <- levels(seurat_data$sample)
seurat_data$sample <- as.character(seurat_data$sample)
seurat_data$sample <- factor(seurat_data$sample, levels = sample_levels)

# Umap highlighting only the twins
twin_umap <- plotDimRed(seurat_data, col_by = "sample", color = sample_colors,
                        plot_type = "pca.umap", highlight_group = TRUE,
                        group = c("107", "108", "109"),
                        meta_data_col = "sample", ggrastr = TRUE)[[1]] +
  ggplot2::ylab("UMAP 2") +
  ggplot2::xlab("UMAP 1")


pdf(file.path(save_dir, "images", "final_figures",
              "4D_twin_UMAP.pdf"), 
    height = 6, width = 6) 

print(twin_umap)
dev.off()

# Overlaps of DE genes

## Figure 5 --------------------------------------------------------------------
# Make a barplot of antigen 
make_barplot <- function(seurat_object, isotype_colors,
                         meta_data_col = "imcantation_isotype"){
  seurat_object$full_split <- paste(seurat_object$celltype_cluster,
                                    seurat_object$tet_name_cutoff, sep = "__")
  
  barplot <- stacked_barplots(seurat_object, meta_col = meta_data_col,
                              split_by = "full_split",
                              color = isotype_colors,
                              return_values = TRUE)
  
  barplot_data <- barplot$data
  barplot_data$celltype <- gsub("__.*", "", barplot_data$split_by)
  barplot_data$tet_name_cutoff <- gsub(".*__", "", barplot_data$split_by)
  barplot_data[[meta_data_col]] <- barplot_data$meta_col
  
  barplot_data$Status <- factor(barplot_data$tet_name_cutoff,
                                levels = plotting_levels)

  barplot_data$celltype <- factor(barplot_data$celltype,
                                  levels = levels(seurat_data$celltype_cluster))
  
  barplot <- ggplot2::ggplot(barplot_data, aes(x = Status, y = percents,
                                               fill = !!sym(meta_data_col))) +
    ggplot2::geom_bar(position = "stack", stat = "identity") +
    ggplot2::facet_grid(~celltype, switch = "x") +
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

pdf(file.path(image_dir, "5A_stacked_bar_chart_isotype_celltype.pdf"),
    width = 12, height = 8)

print(barplot)

dev.off()
graphics.off()

## Figure 7 --------------------------------------------------------------------

clone_data <- file.path(save_dir, "files", "v_gene_counting", 
                        "clone_expansion.xlsx")


clone_expanded <- openxlsx::readWorkbook(clone_data, sheet = "expanded_clones")
clone_public <- openxlsx::readWorkbook(clone_data, sheet = "public_clones")

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

