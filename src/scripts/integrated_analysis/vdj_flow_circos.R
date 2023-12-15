# This document makes VDJ plots
# 1. Circos plots of V + J connections (in all B cells and separated by isotype,
# sample and tetramer)
# 2. Heatmaps of V + J connection (in all B cells and separated by isotype,
# sample and tetramer)
# 3. Barplots of the fraction of V+J pairs within samples and overall
# 4. Statistics showing the odds ratios and fishers exact test of the counts
# of each v gene per comparison and t test of all the proportions
# 5. Polar plots showing the fraction of each V gene per sample (in all
# B cells and separated by isotype, sample, and tetramer)
# Still to add --> 
# 1. Flow plots that show the same data as the circos plots
# 2. All plots across light V+J pairings as well
# 3. Flow plots, circos plots, and heatmaps showing heavy + light V pairings
library(here)
library(scAnalysisR)
library(pheatmap)
library(tidyverse)
library(splitstackshape)
library(circlize)
library(viridis)
library(plotly)
library(ggalluvial)

normalization_method <- "log" # can be SCT or log
# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

#sample <- args[[1]]
#sample <- gsub("__.*", "", sample)
sample <- "merged"

#results_dir <- args[[2]]
results_dir <- here("results")

#sample_info <- args[[4]]
sample_info <- here("files/sample_info.tsv")

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

# This can be added into the snakemake pipeline, will need to write a new
# rule
#all_groups <- c("smh", "significant_v", "smh_significant_v")
#all_groups <- c("all", "memory", "smh", "significant_v")

all_groups <- c("smh_3_plus", "odds_significant_all",
                "odds_significant_no", "odds_significant_aab1",
                "odds_significant_aab2", "smh_significant_v_all",
                "smh_significant_v_no", "smh_significant_v_aab1",
                "smh_significant_v_aab2")   


#cells_use <- "memory"
get_smh <- function(seurat_data){
  # Pull out v_gene information
  sep_columns <- c("chains", "cdr3", "cdr3_length",
                   "cdr3_nt_length", "v_gene", "d_gene", "j_gene", "c_gene",
                   "reads", "umis", "productive", "full_length",
                   "v_ins", "v_del", "v_mis", "d_ins", "d_del",
                   "d_mis", "j_ins", "j_del", "j_mis", "c_ins", "c_del", "c_mis",
                   "all_ins", "all_del", "all_mis", "vd_ins", "vd_del", "dj_ins",
                   "dj_del", "v_mis_freq", "d_mis_freq", "j_mis_freq",
                   "c_mis_freq", "all_mis_freq")
  keep_columns <- c("isotype", "RNA_combined_celltype", "sample", "paired",
                    "clonotype_id", "Status", "tet_hash_id", "all_chains")
  
  all_info <- seurat_data[[]] %>%
    dplyr::mutate(all_chains = chains) %>%
    dplyr::select(dplyr::all_of(c(sep_columns, keep_columns)), n_chains) %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::filter(!is.na(n_chains))
  
  
  all_info_split <- cSplit(all_info, sep_columns, sep = ";", direction = "long") %>%
    dplyr::filter(!is.na(chains))
  
  all_info_split <- all_info_split %>%
    dplyr::filter(isotype %in% c("IGHA", "IGHD", "IGHE", "IGHG"),
                  v_mis > 2)
  
  seurat_data <- subset(seurat_data, cells = unique(all_info_split$barcode))
  
  return(seurat_data)
}

get_sig_v <- function(seurat_data, tests_use, save_dir){
  # Read in odds results
  tests_use <- c("no_nd", "aab1_nd", "aab2_nd")
  all_v_genes <- lapply(tests_use, function(x){
    
    # Keep only genes that were used significantly more highly in 
    # the condition vs control
    odds_data <- openxlsx::read.xlsx(file.path(save_dir, "files", "vdj_files",
                                               "stats", paste0(x, ".xlsx")),
                                     sheet = "odds_ratio") %>%
      dplyr::filter(p_adj < 0.05, odds_ratio > 1)
    
    return(odds_data$v_gene)
    
  })
  
  all_v_genes <- unique(unlist(all_v_genes))
  
  seurat_meta <- seurat_data[[]]
  
  filtered_data <- seurat_meta[grepl(paste(all_v_genes, collapse = "|"), 
                                     seurat_meta$v_gene, 
                                     ignore.case = TRUE, perl = TRUE), ]
  
  seurat_data <- subset(seurat_data, cells = rownames(filtered_data))
  
}

image_save_dir <- file.path(save_dir, "images",
                            "vdj_analysis")

file_save_dir <- file.path(save_dir, "files",
                           "vdj_analysis")

ifelse(!dir.exists(image_save_dir), dir.create(image_save_dir), FALSE)

ifelse(!dir.exists(file_save_dir), dir.create(file_save_dir), FALSE)


for(cells_use in c(all_groups)){
  
  # Read in data
  seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))
  
  if(cells_use == "all"){
    vdj_dir <- file.path(image_save_dir, "heavy_light_plots")
    vdj_files <- file.path(file_save_dir, "heavy_light_files")
    seurat_data <- subset(seurat_data, 
                          subset = RNA_combined_celltype %in% 
                            c("Activated_memory", "Memory_IgA",
                              "Resting_memory", "B.intermediate",
                              "BND2", "Naive_1", "Naive_3", "Plasmablast"))
  } else if (cells_use == "memory") {
    vdj_dir <- file.path(image_save_dir, "memory_vdj_plots")
    vdj_files <- file.path(file_save_dir, "memory_vdj_files")
    seurat_data <- subset(seurat_data, 
                          subset = RNA_combined_celltype %in% 
                            c("Activated_memory", "Memory_IgA",
                              "Resting_memory"))
  } else if (grepl("odds_significant", cells_use)){
    vdj_dir <- file.path(image_save_dir, 
                         paste0(cells_use, "_vdj_plots"))
    vdj_files <- file.path(file_save_dir,
                           paste0(cells_use, "_vdj_files"))
    
    tests_use <- gsub("odds_significant_", "", cells_use)
    if(tests_use == "all"){
      tests_use <- c("no_nd", "aab1_nd", "aab2_nd")
    }
    
    seurat_data <- get_sig_v(seurat_data = seurat_data,
                             tests_use = tests_use,
                             save_dir = save_dir)
    
  } else if(cells_use == "smh_3_plus") {
    vdj_dir <- file.path(image_save_dir, 
                         paste0(cells_use, "_vdj_plots"))
    vdj_files <- file.path(file_save_dir, 
                           paste0(cells_use, "_vdj_files"))
    
    seurat_data <- get_smh(seurat_data)
    
    
  } else if(grepl("smh_significant_v", cells_use)) {
    vdj_dir <- file.path(image_save_dir, 
                         paste0(cells_use, "_v_vdj_plots"))
    vdj_files <- file.path(file_save_dir,
                           paste0(cells_use, "smh_3_plus_sig_v_vdj_files"))
    
    seurat_data <- get_smh(seurat_data)
    
    tests_use <- gsub("smh_significant_v_", "", cells_use)
    if(tests_use == "all"){
      tests_use <- c("no_nd", "aab1_nd", "aab2_nd")
    }
    
    seurat_data <- get_sig_v(seurat_data = seurat_data,
                             tests_use = tests_use,
                             save_dir = save_dir)
  }
  
  
  
  
  ifelse(!dir.exists(vdj_dir), dir.create(vdj_dir), FALSE)
  ifelse(!dir.exists(vdj_files), dir.create(vdj_files), FALSE)
  
  
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
  
  # Analysis ---------------------------------------------------------------------
  
  # The clones com from Chang-o define clones
  # based only on heavy or light
  # Based on the distance found in 07_run_immcantation using shazam
  # Note, these appear to be found within samples individually
  # clone_info <- read.table(file.path(save_dir, "define_clones",
  #                                    "immcantation_combined_clone-pass.tsv"),
  #                          sep = "\t", header = TRUE)
  # 
  # sample_clone_info <- read.table(file.path(save_dir, "define_clones",
  #                                           "immcantation_combined_clone-pass_sample.tsv"),
  #                                 sep = "\t", header = TRUE)
  # 
  # sample_clone_info <- sample_clone_info %>%
  #   dplyr::select(sequence_id, clone_id, locus, sample) %>%
  #   dplyr::rename(sample_clone_id = clone_id,
  #                 sample_locus = locus) %>%
  #   dplyr::mutate(sequence_id_sample = paste(sequence_id,
  #                                            sample, sep = "_")) %>%
  #   dplyr::select(-sequence_id, -sample)
  # 
  # clone_info <- clone_info %>%
  #   dplyr::mutate(sequence_id_sample = paste(sequence_id,
  #                                            sample, sep = "_")) %>%
  #   merge(sample_clone_info, by = "sequence_id_sample", all.x = TRUE)
  # 
  
  sep_columns <- c("chains", "cdr3", "cdr3_length",
                   "cdr3_nt_length", "v_gene", "d_gene", "j_gene", "c_gene",
                   "reads", "umis", "productive", "full_length",
                   "v_ins", "v_del", "v_mis", "d_ins", "d_del",
                   "d_mis", "j_ins", "j_del", "j_mis", "c_ins", "c_del", "c_mis",
                   "all_ins", "all_del", "all_mis", "vd_ins", "vd_del", "dj_ins",
                   "dj_del", "v_mis_freq", "d_mis_freq", "j_mis_freq",
                   "c_mis_freq", "all_mis_freq")
  keep_columns <- c("isotype", "RNA_combined_celltype", "sample", "paired",
                    "clonotype_id", "Status", "tet_hash_id", "all_chains")
  
  all_info <- seurat_data[[]] %>%
    dplyr::mutate(all_chains = chains) %>%
    dplyr::select(dplyr::all_of(c(sep_columns, keep_columns)), n_chains) %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::filter(!is.na(n_chains))
  
  
  all_info_split <- cSplit(all_info, sep_columns, sep = ";", direction = "long") %>%
    dplyr::filter(!is.na(chains))
  
  all_v <- unique(all_info_split$v_gene)
  all_j <- unique(all_info_split$j_gene)
  
  palette1 <- colorRampPalette(colors = 
                                 RColorBrewer::brewer.pal(name = "Set1", n = 9))
  
  palette2 <- colorRampPalette(colors = 
                                 RColorBrewer::brewer.pal(name = "Set2", n = 8))
  
  v_colors <- palette1(length(all_v))
  names(v_colors) <- all_v
  
  j_colors <- palette2(length(all_j))
  names(j_colors) <- all_j
  
  all_colors <- c(v_colors, j_colors)
  
  # Analysis ---------------------------------------------------------------------
  
  # Circos plots to make
  # 1. Heavy V+J pairs
  # 2. Light V+J pairs
  # 3. Heavy + Light V+J pairs
  # 4. Heavy VJ pairs between samples
  # 5. Light VJ pairs between samples
  # 6. Heavy + Light VJ pairs between samples
  # For each, also separate by chain...
  
  # Functions
  
  make_circos_plot <- function(save_name, circos_df, color = NULL,
                               grid_color = NULL){
    if(nrow(circos_df) > 0){
      pdf(save_name, height = 17, width = 17)
      
      chordDiagram(circos_df, annotationTrack = "grid", 
                   preAllocateTracks = 1, col = color,
                   grid.col = grid_color)
      
      # we go back to the first track and customize sector labels
      circos.track(track.index = 1, panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                    facing = "clockwise", niceFacing = TRUE, adj = c(-0.1, 0.5))
      }, bg.border = NA) # here set bg.border to NA is important
      
      dev.off()
      
      circos.clear()    
    }
    
  }
  
  
  make_alluvial_plot <- function(save_name, alluvial_df, color = NULL,
                                 plot_facet = FALSE){
    if(nrow(alluvial_df) > 0){
      min_label <- sum(alluvial_df$value) / 100
      
      if(is.null(color)){
        palette1 <- colorRampPalette(colors = 
                                       RColorBrewer::brewer.pal(name = "Set1", n = 9))
        
        color <- palette1(length(unique(alluvial_df$col_by)))
      }
      
      alluvial_df <- alluvial_df %>%
        dplyr::group_by(v_gene) %>%
        dplyr::mutate(total_v = sum(value)) %>%
        dplyr::group_by(j_gene) %>%
        dplyr::mutate(total_j = sum(value)) %>%
        dplyr::ungroup()
      
      v_levels <- alluvial_df %>%
        dplyr::select(v_gene, total_v) %>%
        dplyr::distinct() %>%
        dplyr::arrange(desc(total_v))
      
      j_levels <- alluvial_df %>%
        dplyr::select(j_gene, total_j) %>%
        dplyr::distinct() %>%
        dplyr::arrange(desc(total_j))
      
      alluvial_df$v_gene <- factor(alluvial_df$v_gene, levels = v_levels$v_gene)
      alluvial_df$j_gene <- factor(alluvial_df$j_gene, levels = j_levels$j_gene)
      
      pdf(save_name, height = 17, width = 17)
      
      save_plot <- ggplot2::ggplot(alluvial_df, ggplot2::aes(y = value,
                                                             axis1 = v_gene,
                                                             axis2 = j_gene, 
                                                             fill = col_by)) +
        ggalluvial::geom_flow(width = 1/12) +
        ggalluvial::geom_stratum(width = 1/12, fill = "black", color = "grey") +
        ggplot2::geom_text(stat = "stratum", 
                           ggplot2::aes(label = ggplot2::after_stat(stratum)),
                           size = 4, min.y = min_label, color = "white") +
        ggplot2::scale_x_discrete(limits = c("v_gene", "j_gene"),
                                  expand = c(.05, .05)) +
        ggplot2::scale_fill_manual(values = color) 
      
      print(save_plot)
      
      if(plot_facet){
        save_plot_2 <- lapply(unique(alluvial_df$col_by), function(x){
          
          alluvial_df$new_color <- ifelse(alluvial_df$col_by == x, x, "other")
          
          color <- c(color, "other" = "#FFFFFF")
          
          return_plot <- ggplot2::ggplot(alluvial_df, ggplot2::aes(y = value,
                                                                   axis1 = v_gene,
                                                                   axis2 = j_gene, 
                                                                   fill = new_color)) +
            ggalluvial::geom_flow(width = 1/12) +
            ggalluvial::geom_stratum(width = 1/12, fill = "black", color = "grey") +
            ggplot2::geom_text(stat = "stratum", 
                               ggplot2::aes(label = ggplot2::after_stat(stratum)),
                               size = 4, min.y = min_label, color = "white") +
            ggplot2::scale_x_discrete(limits = c("v_gene", "j_gene"),
                                      expand = c(.05, .05)) +
            ggplot2::scale_fill_manual(values = color) 
        })
        
        print(cowplot::plot_grid(plotlist = save_plot_2))
        
      }
      
      dev.off()
      
    }
    
  }
  
  
  count_genes <- function(starting_df, group_by = "all",
                          chain = "IGH",
                          keep_chains = c("IGH", "IGH;IGK",
                                          "IGH;IGK;IGK", "IGH;IGK;IGL",
                                          "IGH;IGL", "IGH;IGL;IGL"),
                          color_list = NULL, subset_counts = 0){
    if(group_by == "all"){
      select_cols <- c("v_gene", "j_gene")
    } else if (group_by == "sample") {
      select_cols <- c("v_gene", "j_gene", "sample")
    } else if (group_by == "Status") {
      select_cols <- c("v_gene", "j_gene", "Status")
    } else {
      stop("group_by must be 'all', 'sample', or 'Status'")
    }
    
    return_df <- starting_df %>%
      dplyr::filter(all_chains %in% keep_chains) %>%
      dplyr::filter(chains == chain) %>%
      dplyr::select(dplyr::all_of(select_cols))
    
    return_df <- return_df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(select_cols))) %>% 
      dplyr::add_count(name = "value") %>%
      dplyr::distinct() %>%
      dplyr::ungroup() %>%
      dplyr::filter(value >= subset_counts)
    
    if(group_by != "all"){
      return_df$color <- color_list[return_df[[group_by]]]
      color_values <- return_df$color
      return_df <- return_df %>%
        dplyr::select(-dplyr::all_of(group_by), -color)
    } else{
      color_values <- NULL
    }
    
    return(list("df" = return_df, "color_list" = color_values))
  }
  
  count_genes_heavy_light <- function(starting_df, group_by = "all",
                                      chain = c("IGK", "IGL"),
                                      color_list = NULL, subset_counts = 0){
    
    # Hard coding this because the function won't work if there are multiple
    # heavy and light chains.
    keep_chains <- c("IGH;IGK", "IGH;IGL")
    
    if(group_by == "all"){
      select_cols <- c()
    } else if (group_by == "sample") {
      select_cols <- c("sample")
    } else if (group_by == "Status") {
      select_cols <- c("Status")
    } else {
      stop("group_by must be 'all', 'sample', or 'Status'")
    }
    
    return_df <- starting_df %>%
      dplyr::filter(all_chains %in% keep_chains) %>%
      dplyr::select(dplyr::all_of(select_cols), chains, barcode, v_gene) %>%
      tidyr::pivot_wider(names_from = chains, values_from = v_gene)
    
    if(identical(sort(chain), sort(c("IGK", "IGL")))){
      return_df <- return_df %>%
        dplyr::mutate(IGK_IGL = ifelse(!is.na(IGK), IGK, IGL)) %>%
        dplyr::select(dplyr::all_of(select_cols), IGH, IGK_IGL)
      test_chain <- "IGK_IGL"
    } else {
      return_df <- return_df %>%
        dplyr::select(dplyr::all_of(c(select_cols, "IGH", chain)))
      test_chain <- chain
    }
    
    
    return_df <- return_df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(select_cols, test_chain, "IGH")))) %>% 
      dplyr::add_count(name = "value") %>%
      dplyr::distinct() %>%
      dplyr::ungroup() %>%
      dplyr::filter(value >= subset_counts)
    
    if(group_by != "all"){
      return_df$color <- color_list[return_df[[group_by]]]
      color_values <- return_df$color
      return_df <- return_df %>%
        dplyr::select(-dplyr::all_of(group_by), -color)
    } else{
      color_values <- NULL
    }
    
    return(list("df" = return_df, "color_list" = color_values))
  }
  
  graphics.off()
  # Circos -----------------------------------------------------------------------
  
  test_circos <- list("all" = list(directory = file.path(vdj_dir, "all_pairs"),
                                   subset_counts = 0),
                      "10_plus" = list(directory = file.path(vdj_dir, 
                                                             "vj_counts_10_plus"),
                                       subset_counts = 10))
  
  for(i in names(test_circos)){
    directory <- test_circos[[i]]$directory
    subset_counts <- test_circos[[i]]$subset_counts
    
    ifelse(!dir.exists(directory), dir.create(directory), FALSE)
    
    ## Heavy V+J pairs -------------------------------------------------------------
    # I'll start with 10x
    all_data <- count_genes(starting_df = all_info_split,
                            group_by = "all",
                            subset_counts = subset_counts)
    
    make_circos_plot(save_name = file.path(directory,
                                           "heavy_v_j_circos.pdf"),
                     circos_df = all_data$df,
                     color = all_data$color_list, 
                     grid_color = all_colors)
    
    alluvial_df <- all_data$df
    alluvial_df$col_by <- alluvial_df$v_gene
    
    make_alluvial_plot(save_name = file.path(directory,
                                             "heavy_v_j_flow.pdf"),
                       alluvial_df = alluvial_df,
                       color = all_colors)
    
    ### Repeat by sample ---------------------------------------------------------
    sample_data <- count_genes(starting_df = all_info_split,
                               group_by = "sample", color_list = sample_colors,
                               subset_counts = subset_counts)
    
    make_circos_plot(save_name = file.path(directory,
                                           "heavy_v_j_circos_sample_color.pdf"),
                     circos_df = sample_data$df,
                     color = sample_data$color_list, 
                     grid_color = all_colors)
    
    alluvial_df <- sample_data$df
    alluvial_df$color <- sample_data$color_list
    alluvial_df$col_by <- names(sample_data$color_list)
    
    make_alluvial_plot(save_name = file.path(directory,
                                             "heavy_v_j_flow_sample_color.pdf"),
                       alluvial_df = alluvial_df,
                       color =  sample_data$color_list[!duplicated(names(sample_data$color_list))],
                       plot_facet = TRUE)
    
    ### Repeat by status ---------------------------------------------------------
    status_data <- count_genes(starting_df = all_info_split,
                               group_by = "Status", color_list = status_colors,
                               subset_counts = subset_counts)
    
    make_circos_plot(save_name = file.path(directory,
                                           "heavy_v_j_circos_status_color.pdf"),
                     circos_df = status_data$df,
                     color = status_data$color_list, 
                     grid_color = all_colors)
    
    
    alluvial_df <- status_data$df
    alluvial_df$color <- status_data$color_list
    alluvial_df$col_by <- names(status_data$color_list)
    
    make_alluvial_plot(save_name = file.path(directory,
                                             "heavy_v_j_flow_status_color.pdf"),
                       alluvial_df = alluvial_df,
                       color =  status_data$color_list[!duplicated(names(status_data$color_list))],
                       plot_facet = TRUE)
    
    
    ## Separate by Status --------------------------------------------------------
    invisible(lapply(unique(all_info_split$Status), function(x){
      test_data <- all_info_split %>%
        dplyr::filter(Status == x)
      
      ## All samples
      all_data <- count_genes(starting_df = test_data,
                              group_by = "all",
                              subset_counts = subset_counts)
      
      make_circos_plot(save_name = file.path(directory,
                                             paste0(x, "_heavy_v_j_circos.pdf")),
                       circos_df = all_data$df,
                       color = all_data$color_list, 
                       grid_color = all_colors)
      
      alluvial_df <- all_data$df
      alluvial_df$col_by <- alluvial_df$v_gene
      
      make_alluvial_plot(save_name = file.path(directory,
                                               paste0(x, "heavy_v_j_flow.pdf")),
                         alluvial_df = alluvial_df,
                         color = all_colors)
      
      # Repeat by sample
      sample_data <- count_genes(starting_df = test_data,
                                 group_by = "sample", color_list = sample_colors,
                                 subset_counts = subset_counts)
      
      make_circos_plot(save_name = file.path(directory,
                                             paste0(x, "_heavy_v_j_circos_sample_color.pdf")),
                       circos_df = sample_data$df,
                       color = sample_data$color_list, 
                       grid_color = all_colors)
      
      alluvial_df <- sample_data$df
      alluvial_df$color <- sample_data$color_list
      alluvial_df$col_by <- names(sample_data$color_list)
      
      make_alluvial_plot(save_name = file.path(directory,
                                               paste0(x, "heavy_v_j_flow_sample_color.pdf")),
                         alluvial_df = alluvial_df,
                         color =  sample_data$color_list[!duplicated(names(sample_data$color_list))],
                         plot_facet = TRUE)
      
      # Repeat by status
      sample_data <- count_genes(starting_df = test_data,
                                 group_by = "Status", color_list = status_colors,
                                 subset_counts = subset_counts)
      
      make_circos_plot(save_name = file.path(directory,
                                             paste0(x, "_heavy_v_j_circos_status_color.pdf")),
                       circos_df = alluvial_df,
                       color = sample_data$color_list, 
                       grid_color = all_colors)
      
      alluvial_df <- sample_data$df
      alluvial_df$color <- sample_data$color_list
      alluvial_df$col_by <- names(sample_data$color_list)
      
      make_alluvial_plot(save_name = file.path(directory,
                                               paste0(x, "heavy_v_j_flow_status_color.pdf")),
                         alluvial_df = alluvial_df,
                         color =  sample_data$color_list[!duplicated(names(sample_data$color_list))],
                         plot_facet = TRUE)
    }))
    
    
    
    ## Isotype switched status ---------------------------------------------------
    if(grepl("smh", cells_use)){
      isotypes_use <-  c("IGHA", "IGHD", "IGHG") 
    } else {
      isotypes_use <- c("IGHA", "IGHD", "IGHG", "IGHM")
    }
    
    invisible(lapply(isotypes_use, function(x){
      test_data <- all_info_split %>%
        dplyr::filter(isotype == x)
      
      ## All samples
      all_data <- count_genes(starting_df = test_data,
                              group_by = "all",
                              subset_counts = subset_counts)
      
      make_circos_plot(save_name = file.path(directory,
                                             paste0(x, "_heavy_v_j_circos.pdf")),
                       circos_df = all_data$df,
                       color = all_data$color_list, 
                       grid_color = all_colors)
      
      alluvial_df <- all_data$df
      alluvial_df$col_by <- alluvial_df$v_gene
      
      make_alluvial_plot(save_name = file.path(directory,
                                               paste0(x, "heavy_v_j_flow.pdf")),
                         alluvial_df = alluvial_df,
                         color = all_colors)
      
      
      # Repeat by sample
      sample_data <- count_genes(starting_df = test_data,
                                 group_by = "sample", color_list = sample_colors,
                                 subset_counts = subset_counts)
      
      make_circos_plot(save_name = file.path(directory,
                                             paste0(x, "_heavy_v_j_circos_sample_color.pdf")),
                       circos_df = sample_data$df,
                       color = sample_data$color_list, 
                       grid_color = all_colors)
      
      alluvial_df <- sample_data$df
      alluvial_df$color <- sample_data$color_list
      alluvial_df$col_by <- names(sample_data$color_list)
      
      make_alluvial_plot(save_name = file.path(directory,
                                               paste0(x, "heavy_v_j_flow_sample_color.pdf")),
                         alluvial_df = alluvial_df,
                         color =  sample_data$color_list[!duplicated(names(sample_data$color_list))],
                         plot_facet = TRUE)
      
      # Repeat by status
      sample_data <- count_genes(starting_df = test_data,
                                 group_by = "Status", color_list = status_colors,
                                 subset_counts = subset_counts)
      
      make_circos_plot(save_name = file.path(directory,
                                             paste0(x, "_heavy_v_j_circos_status_color.pdf")),
                       circos_df = sample_data$df,
                       color = sample_data$color_list, 
                       grid_color = all_colors)
      
      alluvial_df <- sample_data$df
      alluvial_df$color <- sample_data$color_list
      alluvial_df$col_by <- names(sample_data$color_list)
      
      make_alluvial_plot(save_name = file.path(directory,
                                               paste0(x, "heavy_v_j_flow_status_color.pdf")),
                         alluvial_df = alluvial_df,
                         color =  sample_data$color_list[!duplicated(names(sample_data$color_list))],
                         plot_facet = TRUE)
    }))
    
    ## Tetramer ------------------------------------------------------------------
    tetramers_use <- c("DNA-tet", "Doublet", "GAD-tet", "IA2-tet",
                       "INS-tet", "Negative", "TET-tet")
    invisible(lapply(tetramers_use, function(x){
      test_data <- all_info_split %>%
        dplyr::filter(tet_hash_id == x)
      
      ## All samples
      all_data <- count_genes(starting_df = test_data,
                              group_by = "all",
                              subset_counts = subset_counts)
      
      make_circos_plot(save_name = file.path(directory,
                                             paste0(x, "_heavy_v_j_circos.pdf")),
                       circos_df = all_data$df,
                       color = all_data$color_list, 
                       grid_color = all_colors)
      
      alluvial_df <- all_data$df
      alluvial_df$col_by <- alluvial_df$v_gene
      
      make_alluvial_plot(save_name = file.path(directory,
                                               paste0(x, "heavy_v_j_flow.pdf")),
                         alluvial_df = alluvial_df,
                         color = all_colors)
      
      
      # Repeat by sample
      sample_data <- count_genes(starting_df = test_data,
                                 group_by = "sample", color_list = sample_colors,
                                 subset_counts = subset_counts)
      
      make_circos_plot(save_name = file.path(directory,
                                             paste0(x, "_heavy_v_j_circos_sample_color.pdf")),
                       circos_df = sample_data$df,
                       color = sample_data$color_list, 
                       grid_color = all_colors)
      
      alluvial_df <- sample_data$df
      alluvial_df$color <- sample_data$color_list
      alluvial_df$col_by <- names(sample_data$color_list)
      
      make_alluvial_plot(save_name = file.path(directory,
                                               paste0(x, "heavy_v_j_flow_sample_color.pdf")),
                         alluvial_df = alluvial_df,
                         color =  sample_data$color_list[!duplicated(names(sample_data$color_list))],
                         plot_facet = TRUE)
      
      # Repeat by status
      sample_data <- count_genes(starting_df = test_data,
                                 group_by = "Status", color_list = status_colors,
                                 subset_counts = subset_counts)
      
      make_circos_plot(save_name = file.path(directory,
                                             paste0(x, "_heavy_v_j_circos_status_color.pdf")),
                       circos_df = sample_data$df,
                       color = sample_data$color_list, 
                       grid_color = all_colors)
      
      alluvial_df <- sample_data$df
      alluvial_df$color <- sample_data$color_list
      alluvial_df$col_by <- names(sample_data$color_list)
      
      make_alluvial_plot(save_name = file.path(directory,
                                               paste0(x, "heavy_v_j_flow_status_color.pdf")),
                         alluvial_df = alluvial_df,
                         color =  sample_data$color_list[!duplicated(names(sample_data$color_list))],
                         plot_facet = TRUE)
    }))
    
  }
  graphics.off()
  
  
  # Heatmap ----------------------------------------------------------------------
  make_heatmap_info <- function(all_info_split, select_cols,
                                group_cols = NULL, type = NULL,
                                subset_counts = 0,
                                keep_chains = c("IGH", "IGH;IGK", "IGH;IGK;IGK",
                                                "IGH;IGK;IGL", "IGH;IGL", 
                                                "IGH;IGL;IGL"),
                                chains_use = c("IGH")){
    
    if(is.null(group_cols)){
      group_cols <- select_cols
    }
    all_info_split_heatmap <- all_info_split %>%
      dplyr::filter(all_chains %in% keep_chains) %>%
      dplyr::filter(chains %in% chains_use) %>%
      dplyr::select(dplyr::all_of(select_cols))
    
    all_info_split_heatmap <- all_info_split_heatmap %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>% 
      dplyr::add_count(name = "value") %>%
      dplyr::distinct()
    
    if(is.null(type)){
      all_info_split_heatmap <- all_info_split_heatmap %>%
        dplyr::mutate(j_group = j_gene) %>%
        ungroup()
    } else if(type == "sample"){
      all_info_split_heatmap <- all_info_split_heatmap %>%
        dplyr::mutate(j_group = paste(j_gene, sample, sep = "_")) %>%
        ungroup()
    } else if (type == "status") {
      all_info_split_heatmap <- all_info_split_heatmap %>%
        dplyr::mutate(j_group = paste(j_gene, Status, sep = "_")) %>%
        ungroup()
    } else {
      stop("type arugment only excepts NULL, status or sample")
    }
    
    new_info <- all_info_split_heatmap %>%
      dplyr::select(v_gene, j_group, value) %>%
      tidyr::pivot_wider(names_from = j_group, values_from = value) %>%
      tibble::column_to_rownames("v_gene")
    
    new_info[is.na(new_info)] <- 0
    
    # Find max for each row
    row_max <- apply(new_info, 1, max)
    
    new_info <- new_info[row_max > subset_counts, ]
    
    if (is.null(type)){
      sample_info <- NULL
    } else if(type %in% c("sample", "status")){
      select_cols_sample_info <- c("j_group", select_cols)
      select_cols_sample_info <- 
        select_cols_sample_info[select_cols_sample_info != "v_gene"]
      sample_info <- all_info_split_heatmap %>%
        dplyr::select(dplyr::all_of(select_cols_sample_info)) %>%
        dplyr::distinct() %>%
        tibble::column_to_rownames("j_group")
      
      if(!identical(rownames(sample_info), colnames(new_info))){
        new_info <- new_info[ , rownames(sample_info)]
      }  
    }
    
    return(list(heatmap_plot = new_info, sample_info = sample_info))
    
  }
  
  plot_heatmap <- function(new_info, sample_info = NULL,
                           save_name = NULL, coloring = NULL,
                           show_colnames = FALSE){
    pdf(save_name, height = 10, width = 10)
    
    if(nrow(new_info) > 1){
      color_function <- colorRampPalette(c("white", "pink", "red"))
      
      if(is_null(sample_info)){
        pheatmap(new_info, color = color_function(100),
                 show_colnames = show_colnames)
        
        grid::grid.newpage()
        new_info <- t(scale(t(new_info)))
        
        pheatmap(new_info, color = color_function(100),
                 show_colnames = show_colnames)    
      } else {
        pheatmap(new_info, color = color_function(100), annotation_col = sample_info,
                 annotation_colors = coloring, show_colnames = show_colnames)
        
        grid::grid.newpage()
        new_info <- t(scale(t(new_info)))
        
        pheatmap(new_info, color = color_function(100), annotation_col = sample_info,
                 annotation_colors = coloring, show_colnames = show_colnames)
      }
      
    } else {
      return(NULL)
    }
    
    
    dev.off()
  }
  
  for(i in names(test_circos)){
    directory <- test_circos[[i]]$directory
    subset_counts <- test_circos[[i]]$subset_counts
    
    heatmap_res <- make_heatmap_info(all_info_split, 
                                     select_cols = c("v_gene", "j_gene"),
                                     subset_counts = subset_counts)
    
    plot_heatmap(new_info = heatmap_res$heatmap_plot, sample_info = NULL,
                 save_name = file.path(directory, "heavy_v_j_heatmap.pdf"),
                 coloring = NULL, show_colnames = TRUE)
    
    ## Heatmap by sample ---------------------------------------------------------
    
    coloring <- list(sample = sample_colors,
                     j_gene = j_colors,
                     Status = status_colors)
    
    heatmap_res <- make_heatmap_info(all_info_split, 
                                     select_cols = c("v_gene", "j_gene",
                                                     "sample", "Status"),
                                     group_cols = c("v_gene", "j_gene", "sample"),
                                     type = "sample")
    
    plot_heatmap(new_info = heatmap_res$heatmap_plot,
                 sample_info = heatmap_res$sample_info,
                 save_name = file.path(directory, "heavy_v_j_heatmap_sample.pdf"),
                 coloring = coloring)
    
    ## Heatmap by status ---------------------------------------------------------
    coloring <- list(j_gene = j_colors,
                     Status = status_colors)
    
    heatmap_res <- make_heatmap_info(all_info_split, 
                                     select_cols = c("v_gene", "j_gene", "Status"),
                                     group_cols = c("v_gene", "j_gene", "Status"),
                                     type = "status")
    
    
    plot_heatmap(new_info = heatmap_res$heatmap_plot,
                 sample_info = heatmap_res$sample_info,
                 save_name = file.path(directory, "heavy_v_j_heatmap_status.pdf"),
                 coloring = coloring)
    
    ## Separate by Status --------------------------------------------------------
    invisible(lapply(unique(all_info_split$Status), function(x){
      test_data <- all_info_split %>%
        dplyr::filter(Status == x)
      
      ## All samples
      heatmap_res <- make_heatmap_info(test_data, 
                                       select_cols = c("v_gene", "j_gene"),
                                       subset_counts = subset_counts)
      
      plot_heatmap(new_info = heatmap_res$heatmap_plot, sample_info = NULL,
                   save_name = file.path(directory,
                                         paste0(x, "_heavy_v_j_heatmap.pdf")),
                   coloring = NULL, show_colnames = TRUE)
      
      # Repeat by sample
      coloring <- list(sample = sample_colors,
                       j_gene = j_colors,
                       Status = status_colors)
      
      heatmap_res <- make_heatmap_info(test_data, 
                                       select_cols = c("v_gene", "j_gene",
                                                       "sample", "Status"),
                                       group_cols = c("v_gene", "j_gene", "sample"),
                                       type = "sample")
      
      plot_heatmap(new_info = heatmap_res$heatmap_plot,
                   sample_info = heatmap_res$sample_info,
                   save_name = file.path(directory, 
                                         paste0(x, "_heavy_v_j_heatmap_sample.pdf")),
                   coloring = coloring)
      
      # Repeat by status
      coloring <- list(j_gene = j_colors,
                       Status = status_colors)
      
      heatmap_res <- make_heatmap_info(test_data, 
                                       select_cols = c("v_gene", "j_gene", "Status"),
                                       group_cols = c("v_gene", "j_gene", "Status"),
                                       type = "status")
      
      
      plot_heatmap(new_info = heatmap_res$heatmap_plot,
                   sample_info = heatmap_res$sample_info,
                   save_name = file.path(directory,
                                         paste0(x, "_heavy_v_j_heatmap_status.pdf")),
                   coloring = coloring)
    }))
    
    
    
    ## Isotype switched status ---------------------------------------------------
    if(grepl("smh", cells_use)){
      isotypes_use <-  c("IGHA", "IGHD", "IGHG") 
    } else {
      isotypes_use <- c("IGHA", "IGHD", "IGHG", "IGHM")
    }    
    invisible(lapply(isotypes_use, function(x){
      test_data <- all_info_split %>%
        dplyr::filter(isotype == x)
      
      ## All samples
      heatmap_res <- make_heatmap_info(test_data, 
                                       select_cols = c("v_gene", "j_gene"),
                                       subset_counts = subset_counts)
      
      plot_heatmap(new_info = heatmap_res$heatmap_plot, sample_info = NULL,
                   save_name = file.path(directory,
                                         paste0(x, "_heavy_v_j_heatmap.pdf")),
                   coloring = NULL,
                   show_colnames = TRUE)
      
      # Repeat by sample
      coloring <- list(sample = sample_colors,
                       j_gene = j_colors,
                       Status = status_colors)
      
      heatmap_res <- make_heatmap_info(test_data, 
                                       select_cols = c("v_gene", "j_gene",
                                                       "sample", "Status"),
                                       group_cols = c("v_gene", "j_gene", "sample"),
                                       type = "sample")
      
      plot_heatmap(new_info = heatmap_res$heatmap_plot,
                   sample_info = heatmap_res$sample_info,
                   save_name = file.path(directory, 
                                         paste0(x, "_heavy_v_j_heatmap_sample.pdf")),
                   coloring = coloring)
      
      # Repeat by status
      coloring <- list(j_gene = j_colors,
                       Status = status_colors)
      
      heatmap_res <- make_heatmap_info(test_data, 
                                       select_cols = c("v_gene", "j_gene", "Status"),
                                       group_cols = c("v_gene", "j_gene", "Status"),
                                       type = "status")
      
      
      plot_heatmap(new_info = heatmap_res$heatmap_plot,
                   sample_info = heatmap_res$sample_info,
                   save_name = file.path(directory,
                                         paste0(x, "_heavy_v_j_heatmap_status.pdf")),
                   coloring = coloring)
    }))
    
    ## Tetramer ------------------------------------------------------------------
    tetramers_use <- c("DNA-tet", "Doublet", "GAD-tet", "IA2-tet",
                       "INS-tet", "Negative", "TET-tet")
    invisible(lapply(tetramers_use, function(x){
      test_data <- all_info_split %>%
        dplyr::filter(tet_hash_id == x)
      
      ## All samples
      heatmap_res <- make_heatmap_info(test_data, 
                                       select_cols = c("v_gene", "j_gene"),
                                       subset_counts = subset_counts)
      
      plot_heatmap(new_info = heatmap_res$heatmap_plot, sample_info = NULL,
                   save_name = file.path(directory,
                                         paste0(x, "_heavy_v_j_heatmap.pdf")),
                   coloring = NULL, show_colnames = TRUE)
      
      # Repeat by sample
      coloring <- list(sample = sample_colors,
                       j_gene = j_colors,
                       Status = status_colors)
      
      heatmap_res <- make_heatmap_info(test_data, 
                                       select_cols = c("v_gene", "j_gene",
                                                       "sample", "Status"),
                                       group_cols = c("v_gene", "j_gene", "sample"),
                                       type = "sample")
      
      plot_heatmap(new_info = heatmap_res$heatmap_plot,
                   sample_info = heatmap_res$sample_info,
                   save_name = file.path(directory, 
                                         paste0(x, "_heavy_v_j_heatmap_sample.pdf")),
                   coloring = coloring)
      
      # Repeat by status
      coloring <- list(j_gene = j_colors,
                       Status = status_colors)
      
      heatmap_res <- make_heatmap_info(test_data, 
                                       select_cols = c("v_gene", "j_gene", "Status"),
                                       group_cols = c("v_gene", "j_gene", "Status"),
                                       type = "status")
      
      
      plot_heatmap(new_info = heatmap_res$heatmap_plot,
                   sample_info = heatmap_res$sample_info,
                   save_name = file.path(directory,
                                         paste0(x, "_heavy_v_j_heatmap_status.pdf")),
                   coloring = coloring)
    }))
    
  }
  graphics.off()
  
  
}
