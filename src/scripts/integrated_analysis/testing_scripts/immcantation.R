library(here)
library(scAnalysisR)
library(pheatmap)
library(tidyverse)
library(alakazam)
library(shazam)


# Set directories
base_dir <- here()
sample <- "merged"


save_dir <- file.path(base_dir, "results", "R_analysis", sample)


# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))

ifelse(!dir.exists(file.path(save_dir, "images", "vdj_analysis")),
       dir.create(file.path(save_dir, "images", "vdj_analysis")),
       FALSE)

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

immcantation_docs <- lapply(unique(seurat_data$sample), function(x){
  immcantation_path <- file.path(here("results", x, "outs", 
                                      "immcantation", "filtered_contig_igblast_db-pass.tsv"))
  
  vdj_data <- alakazam::readChangeoDb(immcantation_path)
  
  seurat_sample <- subset(seurat_data, subset = sample == x)
  
  seurat_meta <- seurat_sample[[]] %>%
    dplyr::select(sample, ID, Initials, Sex, Status, tet_hash_id,
                  scar_hash_id, RNA_combined_celltype, TET_classification) %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::mutate(barcode = gsub("_[0-9]*", "", barcode))
  
  vdj_data <- vdj_data %>%
    dplyr::mutate(barcode = gsub("_contig_[0-9]", "", sequence_id)) %>%
    dplyr::filter(barcode %in% seurat_meta$barcode)
  
  vdj_data <- merge(vdj_data, seurat_meta, by = "barcode")
  
  return(vdj_data)
  
})

all_meta <- seurat_data[[]]

all_meta$cell_id <- gsub("_[0-9]+", "", rownames(all_meta))
all_meta$full_id <- rownames(all_meta)

all_meta <- all_meta %>%
  dplyr::select(cell_id, sample, full_id)

# Merge with immcantation based on cell id and sample,
# Add meta data to SO

immcantation_docs <- do.call(rbind, immcantation_docs)

immcantation_docs$vj_call <- paste(immcantation_docs$v_call,
                                   immcantation_docs$j_call,
                                   sep = "_")

# Quantify usage at the gene level
gene <- countGenes(immcantation_docs, gene="v_call", 
                   groups="Status", mode="gene")

# Quantify usage at the gene level
gene2 <- countGenes(immcantation_docs, gene="v_call", 
                   groups=c("Status", "sample"), mode="gene")

# Try to find the ones with the greatest differences?
gene2_wide <- gene %>%
  dplyr::select(!seq_count) %>%
  tidyr::pivot_wider(names_from = Status, values_from = seq_freq) %>%
  rowwise() %>%
  dplyr::mutate_all(~replace_na(.,0)) %>%
  dplyr::mutate(no_a2 = abs(no - aab_stage_2),
                no_a1 = abs(no - aab_stage_1),
                no_nd = abs(no - nd),
                nd_a2 = abs(nd - aab_stage_2),
                nd_a1 = abs(nd - aab_stage_1),
                a1_a2 = abs(aab_stage_1 - aab_stage_2)) %>%
  dplyr::mutate(max_diff = max(no_a2, no_a1, no_nd, nd_a2, nd_a1, a1_a2,
                               na.rm = TRUE)) %>%
  dplyr::arrange(desc(max_diff))

individual_plot <- gene2 %>%
  dplyr::filter(gene %in% gene2_wide$gene[1:10])

status_plot <- gene %>%
  dplyr::filter(gene %in% gene2_wide$gene[1:10])


plot_ind <- ggplot(individual_plot, aes(x=gene, y=seq_freq)) +
  theme_bw() +
  ggtitle("IGHV Usage Individual") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Percent of repertoire") +
  xlab("") +
  scale_y_continuous() +
  scale_color_manual(values = status_colors) +
  geom_point(aes(color=Status), size=5, alpha=0.8)

plot_status <- ggplot(status_plot, aes(x=gene, y=seq_freq)) +
  theme_bw() +
  ggtitle("IGHV Usage Status") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Percent of repertoire") +
  xlab("") +
  scale_y_continuous() +
  scale_color_manual(values = status_colors) +
  geom_point(aes(color=Status), size=5, alpha=0.8)


final_plot <- cowplot::plot_grid(plot_ind, plot_status, nrow = 1, ncol = 2)

pdf(file.path(save_dir, "images", "vdj_analysis", "v_gene_usage_all.pdf"),
    height = 8, width = 10)

print(final_plot)

dev.off()

# Subset to just ins/ia2/gad individually and repeat
# Look at multibinders (might be interesting to subset by what type of multibinders they are)

# Split by hash id -------------------------------------------------------------
gene3 <- countGenes(immcantation_docs, gene="v_call", 
                    groups=c("Status", "sample", "TET_classification"),
                    mode="gene")

gene_test <- countGenes(immcantation_docs, gene="v_call", 
                        groups=c("Status", "TET_classification"),
                        mode="gene")

# Try to find the ones with the greatest differences?
gene3_wide <- gene_test %>%
  dplyr::mutate(status_classification = paste(Status, TET_classification,
                                              sep = "_")) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene, seq_freq, status_classification) %>%
  tidyr::pivot_wider(names_from = status_classification, values_from = seq_freq) %>%
  rowwise() %>%
  dplyr::mutate_all(~replace_na(.,0))

# Find genes that are different between the status within a tet ----------------

all_tets <- unique(immcantation_docs$TET_classification)
singlet_tets <- all_tets[!grepl("_", all_tets)]
doublet_tets <- all_tets[grepl("_", all_tets)]

all_comparisons <- lapply(singlet_tets, function(x){

  # Find all columns that have that tet
  tet_cols <- colnames(gene3_wide)[grepl(x, colnames(gene3_wide))]
  
  # Remove any doublets
  tet_cols <- tet_cols[!grepl("-.*-", tet_cols)]
  
  # Find all pairwise comparisons
  comparisons <- combn(x = tet_cols, m = 2)
  
  individual_comparisons <- lapply(1:ncol(comparisons), function(y){
    group1 <- comparisons[1,y]
    group2 <- comparisons[2,y]
    col_name <- paste(group1, group2, sep = "__")
    
    new_data <- data.frame(abs(gene3_wide[[group1]] - gene3_wide[[group2]]))
    colnames(new_data) <- col_name
    rownames(new_data) <- gene3_wide$gene
    
    return(new_data)
  })
  
  individual_comparisons <- do.call(cbind, individual_comparisons)
  
  all_columns <- colnames(individual_comparisons)
  
  individual_comparisons <- individual_comparisons %>%
    tibble::rownames_to_column("gene") %>%
    rowwise() %>%
    dplyr::mutate_all(~replace_na(.,0)) %>%
    dplyr::mutate(max_diff = max(across(all_of(all_columns)),
                                 na.rm = TRUE)) %>%
    dplyr::arrange(desc(max_diff))
  
  individual_plot <- gene3 %>%
    dplyr::filter(gene %in% individual_comparisons$gene[1:10]) %>%
    dplyr::filter(TET_classification == x)
  

  
  status_plot <- gene_test %>% 
    dplyr::filter(gene %in% individual_comparisons$gene[1:10]) %>%
    dplyr::filter(TET_classification == x)
  
  g3 <- ggplot(individual_plot, aes(x=gene, y=seq_freq)) +
    theme_bw() +
    ggtitle(paste0(x, " IGHV Usage Individual")) +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
    ylab("Percent of repertoire") +
    xlab("") +
    scale_y_continuous() +
    scale_color_manual(values = status_colors) +
    geom_point(aes(color=Status), size=5, alpha=0.8)
  
  g4 <- ggplot(status_plot, aes(x=gene, y=seq_freq)) +
    theme_bw() +
    ggtitle(paste0(x, " IGHV Usage Status")) +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
    ylab("Percent of repertoire") +
    xlab("") +
    scale_y_continuous() +
    scale_color_manual(values = status_colors) +
    geom_point(aes(color=Status), size=5, alpha=0.8)
  
  final_plot <- cowplot::plot_grid(g3, g4, nrow = 1, ncol = 2)
  
  
  return(final_plot)


})

# Plots of V genes that are picked based on differences between status
# after seprating out tetramers
pdf(file.path(save_dir, "images", "vdj_analysis", "v_gene_usage_tetramer_between_status.pdf"),
    height = 8, width = 10)

print(all_comparisons)

dev.off()


# Repeat by genes that are different between tets within status ----------------

all_tets <- unique(immcantation_docs$TET_classification)
singlet_tets <- all_tets[!grepl("_", all_tets)]
doublet_tets <- all_tets[grepl("_", all_tets)]

all_status <- unique(immcantation_docs$Status)

all_comparisons <- lapply(all_status, function(x){
  
  # Find all columns that have that tet
  status_cols <- colnames(gene3_wide)[grepl(x, colnames(gene3_wide))]
  
  # Remove any doublets
  status_cols <- status_cols[!grepl("-.*-", status_cols)]
  
  # Find all pairwise comparisons
  comparisons <- combn(x = status_cols, m = 2)
  
  individual_comparisons <- lapply(1:ncol(comparisons), function(y){
    group1 <- comparisons[1,y]
    group2 <- comparisons[2,y]
    col_name <- paste(group1, group2, sep = "__")
    
    new_data <- data.frame(abs(gene3_wide[[group1]] - gene3_wide[[group2]]))
    colnames(new_data) <- col_name
    rownames(new_data) <- gene3_wide$gene
    
    return(new_data)
  })
  
  individual_comparisons <- do.call(cbind, individual_comparisons)
  
  all_columns <- colnames(individual_comparisons)
  
  individual_comparisons <- individual_comparisons %>%
    tibble::rownames_to_column("gene") %>%
    rowwise() %>%
    dplyr::mutate_all(~replace_na(.,0)) %>%
    dplyr::mutate(max_diff = max(across(all_of(all_columns)),
                                 na.rm = TRUE)) %>%
    dplyr::arrange(desc(max_diff))
  
  individual_plot <- gene3 %>%
    dplyr::filter(gene %in% individual_comparisons$gene[1:10]) %>%
    dplyr::filter(Status == x, TET_classification %in% singlet_tets)

  status_plot <- gene_test %>% 
    dplyr::filter(gene %in% individual_comparisons$gene[1:10]) %>%
    dplyr::filter(Status == x, TET_classification %in% singlet_tets)
  
  g3 <- ggplot(individual_plot, aes(x=gene, y=seq_freq)) +
    theme_bw() +
    ggtitle(paste0(x, " IGHV Usage Individual")) +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
    ylab("Percent of repertoire") +
    xlab("") +
    scale_y_continuous() +
    scale_color_manual(values = tetramer_colors) +
    geom_point(aes(color=TET_classification), size=5, alpha=0.8)
  
  g4 <- ggplot(status_plot, aes(x=gene, y=seq_freq)) +
    theme_bw() +
    ggtitle(paste0(x, " IGHV Usage Status")) +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
    ylab("Percent of repertoire") +
    xlab("") +
    scale_y_continuous() +
    scale_color_manual(values = tetramer_colors) +
    geom_point(aes(color=TET_classification), size=5, alpha=0.8)
  
  final_plot <- cowplot::plot_grid(g3, g4, nrow = 1, ncol = 2)
  
  
  return(final_plot)
  
  
})

# Plots of v genes that are selected by differences in tetramers within a
# status.
pdf(file.path(save_dir, "images", "vdj_analysis", "v_gene_usage_tetramer_within_status.pdf"),
    height = 8, width = 10)

print(all_comparisons)

dev.off()

# Split by IgA, IgG, IgM -------------------------------------------------------


# Clone analysis ---------------------------------------------------------------
# Identify clones
# https://shazam.readthedocs.io/en/stable/vignettes/DistToNearest-Vignette/

# Filter out any cells with two heavy chains
immcantation_docs$cell_full <- paste(immcantation_docs$barcode,
                                     immcantation_docs$sample,
                                     sep = "_")

immcantation_filtered <- immcantation_docs %>%
  dplyr::group_by(cell_full, locus) %>%
  dplyr::add_count(name = "locus_count") %>%
  dplyr::filter(locus_count < 2 & locus == "IGH")


dist_sc <- distToNearest(immcantation_filtered, cellIdColumn="cell_full",
                         locusColumn="locus", 
                         VJthenLen=FALSE, onlyHeavy=FALSE,
                         fields = "sample", model = "ham", normalize = "len")

# Find threshold using gmm method
# output <- findThreshold(dist_sc$dist_nearest, method="gmm", model="gamma-gamma",
#                         progress = TRUE)
# 
# # Plot distance histogram, Gaussian fits, and optimum threshold
# plot(output, binwidth=0.02, title="GMM Method: gamma-gamma")
# 
# print(output)

output2 <- findThreshold(dist_sc$dist_nearest, method="gmm", model="gamma-norm",
                        progress = TRUE)

pdf(file.path(save_dir, "images", "shazam_distance.pdf"))

plot(output2, binwidth=0.02, title="GMM Method: gamma-gamma")

p1 <- ggplot(subset(dist_sc, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("Hamming distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=output2@threshold, color="firebrick", linetype=2)
plot(p1)

p2 <- ggplot(subset(dist_sc, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("Hamming distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=output2@threshold, color="firebrick", linetype=2) + 
  ggplot2::facet_grid(sample ~ ., scales = "free_y")
plot(p2)

dev.off()
write.table(immcantation_filtered, 
            file.path(save_dir, "files", "immcantation_combined.tsv"),
            sep = "\t")

write.table(output2@threshold,
            file.path(save_dir, "files", "immcantation_dist_threshold.tsv"),
            sep = "\t")
