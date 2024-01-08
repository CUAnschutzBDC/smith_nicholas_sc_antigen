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

# This can be added into the snakemake pipeline, will need to write a new
# rule
cells_use <- "all"
#cells_use <- "memory"

for(cells_use in c("all", "memory")){
  
  print(cells_use)
  
  # Read in data
  seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))
  
  if(cells_use == "all"){
    vdj_dir <- file.path(save_dir, "images", "vdj_plots")
    vdj_files <- file.path(save_dir, "files", "vdj_files")
    seurat_data <- subset(seurat_data, 
                          subset = RNA_combined_celltype %in% 
                            c("Activated_memory", "Memory_IgA",
                              "Resting_memory", "B.intermediate",
                              "BND2", "Naive_1", "Naive_3", "Plasmablast"))
  } else {
    vdj_dir <- file.path(save_dir, "images", "memory_vdj_plots")
    vdj_files <- file.path(save_dir, "files", "memory_vdj_files")
    seurat_data <- subset(seurat_data, 
                          subset = RNA_combined_celltype %in% 
                            c("Activated_memory", "Memory_IgA",
                              "Resting_memory"))
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
  
  print("circos")
  print("")
  for(i in names(test_circos)){
    print(i)
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
    isotypes_use <- c("IGHA", "IGHD", "IGHG", "IGHM")
    
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
  print("")
  print("heatmap")
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
    
    if(nrow(new_info) > 0){
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
    print(i)
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
    isotypes_use <- c("IGHA", "IGHD", "IGHG", "IGHM")
    
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
  
  # Barplots ---------------------------------------------------------------------
  # Here I want to make barplots that are the percent of the repertoire for 
  # each VJ pair. 
  # Steps
  # 1. Count all vj pairs within each grouping (sample, tetramer, isotype, status)
  # 2. Count all cells within each grouping (sample, tetramer, isotype, status)
  # 3. Use these two to find percent
  # 4. Plot percent -> try to keep colors consistent between all plots...?
  # 5. Repeat this subsetting to sample (group tetramer, isotype, status),
  # tetramer (group sample, isotype, status), isotype (group sample, tetramer,
  # status), status (group sample, tetramer, isotype) 
  # I'm not certain on the subsetting front at the moment
  keep_chains <- c("IGH", "IGH;IGK", "IGH;IGK;IGK",
                   "IGH;IGK;IGL", "IGH;IGL", 
                   "IGH;IGL;IGL")
  
  test_set <- all_info_split %>%
    dplyr::filter(chains %in% c("IGH")) %>% 
    dplyr::filter(all_chains %in% keep_chains) %>%
    dplyr::select(sample, v_gene, j_gene, Status) %>%
    dplyr::mutate(vj_gene = paste(v_gene, j_gene, sep = "_")) %>%
    dplyr::group_by(sample, vj_gene) %>%
    dplyr::add_count(name = "sample_vj_count") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sample) %>%
    dplyr::add_count(name = "sample_count") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(vj_gene) %>%
    dplyr::add_count(name = "full_vj_count") %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(frequency_vj_within_sample = sample_vj_count / sample_count * 100,
                  frequency_vj_by_sample = sample_vj_count / full_vj_count * 100)
  
  
  pdf(file.path(vdj_dir, "vj_barplots_test.pdf"),
      height = 10, width = 20)
  # print(ggplot2::ggplot(test_set, ggplot2::aes(x = sample, y = frequency_vj_within_sample,
  #                                        fill = vj_gene)) +
  #   ggplot2::geom_bar(stat = "identity", position = "stack"))
  
  print(ggplot2::ggplot(test_set, ggplot2::aes(x = vj_gene, y = sample_vj_count,
                                               fill = sample)) +
          ggplot2::geom_bar(stat = "identity", position = "stack") +
          ggplot2::scale_fill_manual(values = sample_colors) + 
          ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
          ggplot2::ggtitle("VJ count colored by sample"))
  
  print(ggplot2::ggplot(test_set, ggplot2::aes(x = vj_gene, y = sample_vj_count,
                                               fill = Status)) +
          ggplot2::geom_bar(stat = "identity", position = "stack") +
          ggplot2::scale_fill_manual(values = status_colors) +
          ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggplot2::ggtitle("VJ count colored by status"))
  
  
  print(ggplot2::ggplot(test_set, ggplot2::aes(x = vj_gene, y = frequency_vj_by_sample,
                                               fill = sample)) +
          ggplot2::geom_bar(stat = "identity", position = "stack") +
          ggplot2::scale_fill_manual(values = sample_colors) +
          ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggplot2::ggtitle("Percent of each vj gene by sample"))
  
  print(ggplot2::ggplot(test_set, ggplot2::aes(x = vj_gene, y = frequency_vj_by_sample,
                                               fill = Status)) +
          ggplot2::geom_bar(stat = "identity", position = "stack") +
          ggplot2::scale_fill_manual(values = status_colors) +
          ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggplot2::ggtitle("Percent of each vj gene by sample"))
  
  
  # Select any genes where one status is highly represented
  subset_data <- test_set %>%
    dplyr::group_by(Status, vj_gene) %>%
    dplyr::mutate(percent_status = sum(frequency_vj_by_sample)) %>%
    dplyr::filter(percent_status > 50)
  
  subset_plot <- test_set %>%
    dplyr::filter(vj_gene %in% unique(subset_data$vj_gene))
  
  print(ggplot2::ggplot(subset_plot, ggplot2::aes(x = vj_gene, y = frequency_vj_by_sample,
                                                  fill = sample)) +
          ggplot2::geom_bar(stat = "identity", position = "stack") +
          ggplot2::scale_fill_manual(values = sample_colors) +
          ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("VJ genes with at least 50% of cells from one status"))
  
  print(ggplot2::ggplot(subset_plot, ggplot2::aes(x = vj_gene, y = frequency_vj_by_sample,
                                                  fill = Status)) +
          ggplot2::geom_bar(stat = "identity", position = "stack") +
          ggplot2::scale_fill_manual(values = status_colors) +
          ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("VJ genes with at least 50% of cells from one status"))
  
  
  # Select any genes where one status is not represented
  subset_data <- test_set %>%
    dplyr::select(Status, vj_gene) %>%
    dplyr::distinct() %>%
    dplyr::group_by(vj_gene) %>%
    dplyr::add_count(name = "status_represented") %>%
    dplyr::filter(status_represented < 4)
  
  subset_plot <- test_set %>%
    dplyr::filter(vj_gene %in% unique(subset_data$vj_gene))
  
  print(ggplot2::ggplot(subset_plot, ggplot2::aes(x = vj_gene, y = frequency_vj_by_sample,
                                                  fill = sample)) +
          ggplot2::geom_bar(stat = "identity", position = "stack") +
          ggplot2::scale_fill_manual(values = sample_colors) +
          ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("VJ genes not seen in all status"))
  
  print(ggplot2::ggplot(subset_plot, ggplot2::aes(x = vj_gene, y = frequency_vj_by_sample,
                                                  fill = Status)) +
          ggplot2::geom_bar(stat = "identity", position = "stack") +
          ggplot2::scale_fill_manual(values = status_colors) +
          ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("VJ genes not see in all status"))
  
  dev.off()
  
  
  # Make a data frame
  # Need
  # 1. Percent of the vj gene (Here, adding across samples should be 100%)
  # 2. Percent of the sample (Here adding across vj genes should be 100%)
  # 3. Number of samples
  # 4. Number of samples per status
  # Do for V and VJ
  
  vdj_data <- all_info_split %>%
    dplyr::filter(chains %in% c("IGH")) %>% 
    dplyr::filter(all_chains %in% keep_chains) %>%
    dplyr::select(sample, v_gene, j_gene, Status) %>%
    dplyr::mutate(vj_gene = paste(v_gene, j_gene, sep = "_")) %>%
    dplyr::group_by(sample, vj_gene) %>%
    dplyr::add_count(name = "sample_vj_count") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sample) %>%
    dplyr::add_count(name = "sample_count") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(vj_gene) %>%
    dplyr::add_count(name = "full_vj_count") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Status, vj_gene) %>%
    dplyr::add_count(name = "status_vj_count") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Status) %>%
    dplyr::add_count(name = "status_count") %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(frequency_vj_within_sample = sample_vj_count / sample_count * 100, 
                  frequency_vj_by_sample = sample_vj_count / full_vj_count * 100,
                  frequency_vj_within_status = status_vj_count / status_count * 100,
                  frequency_vj_by_status = status_vj_count / full_vj_count * 100) %>%
    dplyr::group_by(vj_gene) %>%
    dplyr::add_count(name = "total_samples") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(vj_gene, Status) %>%
    dplyr::add_count(name = "samples_per_condition")
  
  save_data <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = save_data, sheetName = "vj_information")
  openxlsx::writeData(wb = save_data, sheet = "vj_information", x= vdj_data)
  
  # Make this into a function that builds v_data based on whatever subset df,
  # and also makes the lists below
  v_data <- all_info_split %>%
    dplyr::filter(chains %in% c("IGH")) %>% 
    dplyr::filter(all_chains %in% keep_chains) %>%
    dplyr::select(sample, v_gene, Status) %>%
    dplyr::group_by(sample, v_gene) %>%
    dplyr::add_count(name = "sample_v_count") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sample) %>%
    dplyr::add_count(name = "sample_count") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(v_gene) %>%
    dplyr::add_count(name = "full_v_count") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Status, v_gene) %>%
    dplyr::add_count(name = "status_v_count") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Status) %>%
    dplyr::add_count(name = "status_count") %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::mutate(frequency_v_within_sample = sample_v_count / sample_count * 100, 
                  frequency_v_by_sample = sample_v_count / full_v_count * 100,
                  frequency_v_within_status = status_v_count / status_count * 100,
                  frequency_v_by_status = status_v_count / full_v_count * 100) %>%
    dplyr::group_by(v_gene) %>%
    dplyr::add_count(name = "total_samples") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(v_gene, Status) %>%
    dplyr::add_count(name = "samples_per_condition")
  
  openxlsx::addWorksheet(wb = save_data, sheetName = "v_information")
  openxlsx::writeData(wb = save_data, sheet = "v_information", x= v_data)
  
  openxlsx::saveWorkbook(wb = save_data, 
                         file = file.path(vdj_files, "vj_proportions.xlsx"),
                         overwrite = TRUE)
  
  
  # frequency within sample, sum = 100 --> frequency_vj_within_sample
  # frequency within vj all combined, sum = 100 --> frequency_vj_by_sample percent 
  # of that v gene contributed by the sample
  
  # Want frequency_vj_within_sample for polar plot and stats
  
  # Stats ------------------------------------------------------------------------
  # Think about this paper
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9377262/
  # Figure 1
  
  # For the stats
  # Try two ways
  # 1. Find proportions across status (ignoring sample)
  # 2. Use this to compute odds ratio and fisher's exact test --> here the 
  # contengency table will be x = nd no y = v present v absent
  # 3. Use the same contengency table to find an odds ratio
  # 4. Find proportions across samples
  # 5. Perform a t test including the sample information
  
  
  all_samples_no <- v_data %>%
    ungroup() %>%
    dplyr::filter(Status == "no")%>%
    dplyr::select(sample) %>%
    dplyr::distinct()
  
  all_samples_nd <- v_data %>%
    ungroup() %>%
    dplyr::filter(Status == "nd")%>%
    dplyr::select(sample) %>%
    dplyr::distinct()
  
  all_samples_aab1 <- v_data %>%
    ungroup() %>%
    dplyr::filter(Status == "aab_stage_1")%>%
    dplyr::select(sample) %>%
    dplyr::distinct()
  
  all_samples_aab2 <- v_data %>%
    ungroup() %>%
    dplyr::filter(Status == "aab_stage_2")%>%
    dplyr::select(sample) %>%
    dplyr::distinct()
  
  add_zero_vals <- function(all_samples, t_vals, status){
    if(all(all_samples$sample %in% t_vals$sample)){
      return(t_vals)
    } else {
      # Make a data frame with zeros
      missing_samples <- all_samples[!all_samples$sample %in% t_vals$sample,]
      missing_df <- data.frame(sample = missing_samples$sample,
                               Status = status,
                               frequency_v_within_sample = 0)
      
      return_df <- rbind(t_vals, missing_df)
      return(return_df)
    }
  }
  
  add_odds_zero <- function(status_one, status_two, odds_data, v_gene,
                            full_df){
    if(all(c(status_one, status_two) %in% odds_data$Status)){
      return(odds_data)
    } else if(!status_one %in% odds_data$Status){
      status_count_df <- full_df %>%
        dplyr::filter(Status == status_one)
      status_count <- unique(status_count_df$status_count)
      if(length(status_count) != 1){
        stop("Something went wrong with counting")
      }
      new_df <- data.frame(v_gene = v_gene,
                           Status = status_one,
                           status_v_count = 0,
                           status_count = status_count)
      return_df <- rbind(odds_data, new_df)
      return(return_df)
    } else if(!status_two %in% odds_data$Status){
      status_count_df <- full_df %>%
        dplyr::filter(Status == status_two)
      status_count <- unique(status_count_df$status_count)
      if(length(status_count) != 1){
        stop("Something went wrong with counting")
      }
      new_df <- data.frame(v_gene = v_gene,
                           Status = status_two,
                           status_v_count = 0,
                           status_count = status_count)
      return_df <- rbind(odds_data, new_df)
      return(return_df)
    }
  }
  
  v_j_statistics <- function(v_data, status_one, status_two,
                             all_samples_one, all_samples_two,
                             min_v_genes = 5){
    # Let's start with no vs nd
    contengency <- v_data %>%
      dplyr::filter(Status %in% c(status_one, status_two))
    
    v_gene_list <- unique(contengency$v_gene)
    
    all_res <- lapply(v_gene_list, function(x){
      v_gene_check <- contengency %>%
        dplyr::filter(v_gene == x)
      
      # Odds ratio
      # change to status_v_count and status_count
      odds_data <- v_gene_check %>%
        dplyr::select(v_gene, Status, status_v_count, status_count) %>%
        dplyr::distinct()
      
      odds_data <- add_odds_zero(status_one = status_one, status_two = status_two,
                                 odds_data = odds_data, v_gene = x,
                                 full_df = contengency)
      
      a <- as.numeric(odds_data[odds_data$Status == status_one,
                                "status_v_count"])
      b <- as.numeric(odds_data[odds_data$Status == status_two,
                                "status_v_count"])
      total_a <- as.numeric(odds_data[odds_data$Status == status_one,
                                      "status_count"])
      
      total_b <- as.numeric(odds_data[odds_data$Status == status_two,
                                      "status_count"])
      c <- total_a - a
      d <- total_b - b
      
      or <- (a*d) / (b*c)
      
      # Fisher's exact
      my_mat <- matrix(c(a, c, b, d), ncol = 2)
      
      fisher_result <- fisher.test(x = my_mat)    
      
      # T test - use all samples
      one_vals <- v_gene_check[v_gene_check$Status == status_one,]
      two_vals <- v_gene_check[v_gene_check$Status == status_two,]
      
      one_vals <- one_vals %>%
        ungroup() %>%
        dplyr::select(sample, Status, frequency_v_within_sample)
      
      two_vals <- two_vals %>%
        ungroup() %>% 
        dplyr::select(sample, Status, frequency_v_within_sample)
      
      one_vals <- add_zero_vals(all_samples = all_samples_one,
                                t_vals = one_vals,
                                status = status_one)
      
      two_vals <- add_zero_vals(all_samples = all_samples_two,
                                t_vals = two_vals,
                                status = status_two)
      
      t_test_res <- t.test(x = one_vals$frequency_v_within_sample, 
                           y = two_vals$frequency_v_within_sample, 
                           alternative = "two.sided")    
      
      
      one_vals <- one_vals %>%
        ungroup() %>%
        dplyr::select(sample, Status, frequency_v_within_sample)
      
      two_vals <- two_vals %>%
        ungroup() %>% 
        dplyr::select(sample, Status, frequency_v_within_sample)
      
      t_test_df <- rbind(one_vals, two_vals)
      
      # Only return if enough of the v gene are present in either group
      if(sum(a, b) >= min_v_genes){
        return(list("odds_data" = odds_data, 
                    "odds_ratio" = or,
                    "fisher_result" = fisher_result,
                    "t_test_df" = t_test_df,
                    "t_test_res" = t_test_res,
                    "v_gene" = x))      
      } else {
        return(list("odds_data" = odds_data, 
                    "odds_ratio" = NULL,
                    "fisher_result" = NULL,
                    "t_test_df" = t_test_df,
                    "t_test_res" = NULL,
                    "v_gene" = x))
      }
      
      
      
    })
    
    names(all_res) = v_gene_list
    
    # Pull out t test data
    all_t_test <- lapply(all_res, function(x){
      v_gene <- x$v_gene
      t_test_res <- x$t_test_res
      
      # Only return t test information if there were enough counts of the
      # v gene
      if(!is.null(t_test_res)){
        p_value <- t_test_res$p.value
        t <- t_test_res$statistic
        ci_low <- t_test_res$conf.int[[1]]
        ci_high <- t_test_res$conf.int[[2]] 
        
        return_data <- data.frame("v_gene" = v_gene,
                                  "p_value" = p_value,
                                  "t_stat" = t,
                                  "95_conf_int_low" = ci_low,
                                  "95_conf_int_high" = ci_high,
                                  "status_one" = status_one,
                                  "status_two" = status_two)    
        
        
        return(return_data)
      } else {
        return(NULL)
      }
    })
    
    
    all_t_test <- do.call(rbind, all_t_test)
    
    # Pull out t test data for box plot
    all_t_test_plot <- lapply(all_res, function(x){
      v_gene <- x$v_gene
      t_test_res <- x$t_test_df
      if(!is.null(t_test_res)){
        t_test_res$v_gene <- v_gene
        t_test_res$status_one <- status_one
        t_test_res$status_two <- status_two
        
        return(t_test_res)      
      } else {
        return(NULL)
      }
      
    })
    
    
    all_t_test_plot <- do.call(rbind, all_t_test_plot)
    
    
    # Pull out odds ratio
    all_odds_ratio <- lapply(all_res, function(x){
      v_gene <- x$v_gene
      odds_ratio <- x$odds_ratio
      if(!is.null(odds_ratio)){
        fisher_p <- x$fisher_result$p.value
        ci_low <- x$fisher_result$conf.int[[1]]
        ci_high <- x$fisher_result$conf.int[[2]]
        return_data <- data.frame("v_gene" = v_gene,
                                  "p_value" = fisher_p,
                                  "odds_ratio" = odds_ratio,
                                  "95_conf_int_low" = ci_low,
                                  "95_conf_int_high" = ci_high,
                                  "status_one" = status_one,
                                  "status_two" = status_two)     
        
        return(return_data)      
      } else {
        return(NULL)
      }
      
    })
    
    
    all_odds_ratio <- do.call(rbind, all_odds_ratio)
    
    # Pull out odds ratio data for plots
    all_or_plot <- lapply(all_res, function(x){
      v_gene <- x$v_gene
      odds_res <- x$odds_data
      odds_res$v_gene <- v_gene
      odds_res$status_one <- status_one
      odds_res$status_two <- status_two
      
      return(odds_res)
    })
    
    
    all_or_plot <- do.call(rbind, all_or_plot)
    
    all_or_plot$fraction <- all_or_plot$status_v_count / all_or_plot$status_count
    
    return(list(all_t_test = all_t_test,
                all_t_test_plot = all_t_test_plot,
                all_odds_ratio = all_odds_ratio,
                all_or_plot = all_or_plot))
    
  }
  
  # Make a plot of odds ratios
  
  # Plots want - nd vs all for all subsets
  # subsets include: 
  # 1. All tetramer positive cells - don't clump diabetes reactive
  # 2. All subsets of class switching
  
  # Build comparisons
  comparison_builder <- function(starting_df){
    # Make this into a function that builds v_data based on whatever subset df,
    # and also makes the lists below
    v_data <- starting_df %>%
      dplyr::filter(chains %in% c("IGH")) %>% 
      dplyr::filter(all_chains %in% keep_chains) %>%
      dplyr::select(sample, v_gene, Status) %>%
      dplyr::group_by(sample, v_gene) %>%
      dplyr::add_count(name = "sample_v_count") %>%
      dplyr::ungroup() %>%
      dplyr::group_by(sample) %>%
      dplyr::add_count(name = "sample_count") %>%
      dplyr::ungroup() %>%
      dplyr::group_by(v_gene) %>%
      dplyr::add_count(name = "full_v_count") %>%
      dplyr::ungroup() %>%
      dplyr::group_by(Status, v_gene) %>%
      dplyr::add_count(name = "status_v_count") %>%
      dplyr::ungroup() %>%
      dplyr::group_by(Status) %>%
      dplyr::add_count(name = "status_count") %>%
      dplyr::ungroup() %>%
      dplyr::distinct() %>%
      dplyr::mutate(frequency_v_within_sample = sample_v_count / sample_count * 100, 
                    frequency_v_by_sample = sample_v_count / full_v_count * 100,
                    frequency_v_within_status = status_v_count / status_count * 100,
                    frequency_v_by_status = status_v_count / full_v_count * 100) %>%
      dplyr::group_by(v_gene) %>%
      dplyr::add_count(name = "total_samples") %>%
      dplyr::ungroup() %>%
      dplyr::group_by(v_gene, Status) %>%
      dplyr::add_count(name = "samples_per_condition")
    
    return(v_data)
  }
  
  
  all_samples_no <- all_info_split %>%
    ungroup() %>%
    dplyr::filter(Status == "no")%>%
    dplyr::select(sample) %>%
    dplyr::distinct()
  
  all_samples_nd <- all_info_split %>%
    ungroup() %>%
    dplyr::filter(Status == "nd")%>%
    dplyr::select(sample) %>%
    dplyr::distinct()
  
  all_samples_aab1 <- all_info_split %>%
    ungroup() %>%
    dplyr::filter(Status == "aab_stage_1")%>%
    dplyr::select(sample) %>%
    dplyr::distinct()
  
  all_samples_aab2 <- all_info_split %>%
    ungroup() %>%
    dplyr::filter(Status == "aab_stage_2")%>%
    dplyr::select(sample) %>%
    dplyr::distinct()
  
  
  all_tests <- list("no_nd" = c("no", "nd"),
                    "aab1_nd" = c("aab_stage_1", "nd"),
                    "aab2_nd" = c("aab_stage_2", "nd"))
  
  sample_mapping <- list("nd" = all_samples_nd,
                         "no" = all_samples_no,
                         "aab_stage_1" = all_samples_aab1,
                         "aab_stage_2" = all_samples_aab2)
  
  build_tests <- lapply(c("all", "antigen", "isotype"), function(x){
    if(x == "all"){
      v_df <- list(comparison_builder(starting_df = all_info_split))
      names(v_df) <- "all"
    } else if(x == "antigen") {
      v_df <- lapply(unique(all_info_split$tet_hash_id), function(x){
        subset_df <- all_info_split %>%
          dplyr::filter(tet_hash_id == x)
        return(comparison_builder(starting_df = subset_df))
      })
      names(v_df) <- unique(all_info_split$tet_hash_id)
    } else if(x == "isotype"){
      isotype_use <- c("IGHA", "IGHD", "IGHG", "IGHM")
      v_df <- lapply(isotype_use, function(x){
        subset_df <- all_info_split %>%
          dplyr::filter(isotype == x)
        return(comparison_builder(starting_df = subset_df))
      })
      names(v_df) <- isotype_use
    }
    
    # Now we have the v df. We want to build a list that includes that
    # information as well.
    all_test_full <- lapply(names(v_df), function(x){
      individual_test <- lapply(all_tests, function(y){
        return_data <- list(y, v_df[[x]])
      })
      names(individual_test) <- paste(x, names(individual_test), sep = "_")
      return(individual_test)
    })
    
    all_test_full <- do.call("c", all_test_full)
    
    all_stats <- lapply(all_test_full, function(x){
      status_one <- x[[1]][[1]]
      status_two <- x[[1]][[2]]
      all_samples_one <- sample_mapping[[status_one]]
      all_samples_two <- sample_mapping[[status_two]]
      v_j_statistics(x[[2]], status_one = status_one,
                     status_two = status_two,
                     all_samples_one = all_samples_one,
                     all_samples_two = all_samples_two)  
    })
    
    #names(all_stats) <- names(all_tests)
    
    # Do full p-value correction for the stats tests
    all_odds <- lapply(names(all_stats), function(x){
      return_df <- all_stats[[x]]$all_odds_ratio
      return_df$test <- x
      return(return_df)
    })
    
    all_odds <- do.call(rbind, all_odds)
    
    all_odds$p_adj <- p.adjust(p = all_odds$p_value,
                               method = "bonferroni")
    
    all_t <- lapply(names(all_stats), function(x){
      return_df <- all_stats[[x]]$all_t_test
      return_df$test <- x
      return(return_df)
    })
    
    all_t <- do.call(rbind, all_t)
    
    all_t$p_adj <- p.adjust(p = all_t$p_value,
                            method = "bonferroni")
    
    save_dir_stats_files <- file.path(vdj_files, "stats")
    ifelse(!dir.exists(save_dir_stats_files), dir.create(save_dir_stats_files),
           FALSE)
    
    save_dir_stats <- file.path(vdj_dir, "stats")
    ifelse(!dir.exists(save_dir_stats), dir.create(save_dir_stats),
           FALSE)
    
    all_plots <- lapply(names(all_test_full), function(x){
      all_t_test_plot <- all_stats[[x]]$all_t_test_plot %>%
        dplyr::arrange(v_gene)
      # Make a plot of t-test
      t_plot <- ggplot2::ggplot(all_t_test_plot,
                                ggplot2::aes(x = v_gene,
                                             y = frequency_v_within_sample,
                                             fill = Status)) +
        ggplot2::geom_boxplot() +
        ggplot2::scale_fill_manual(values = status_colors) +
        ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))  
      
      # Save data
      pdf(file.path(save_dir_stats, paste0(x, "_t_test.pdf")),
          width = 12, height = 8)
      print(t_plot)
      dev.off()
      
      odds_data <- all_odds %>%
        dplyr::filter(test == x)
      
      odds_data$padj_category <- ifelse(odds_data$p_adj < 0.05, "P<0.05",
                                        "P>0.05")
      col_mapping <- c("P>0.05" = "#D3D3D3",
                       "P<0.05" = "#AA4A44")
      
      colors <- col_mapping[odds_data$padj_category]
      
      odds_plot <- ggplot2::ggplot(odds_data, 
                                   ggplot2::aes(x = odds_ratio,
                                                y = v_gene,
                                                color = padj_category)) +
        ggplot2::geom_point() +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = X95_conf_int_low,
                                            xmax = X95_conf_int_high)) +
        ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
        ggplot2::scale_color_manual(values = col_mapping)
      
      
      # Save data
      pdf(file.path(save_dir_stats, paste0(x, "_odds_ratio.pdf")),
          width = 8, height = 12)
      print(odds_plot)
      dev.off()
      
      odds_heatmap_data <- all_stats[[x]]$all_or_plot %>% 
        dplyr::select(v_gene, Status, fraction) %>%
        tidyr::pivot_wider(names_from = Status, values_from = fraction) %>%
        dplyr::arrange(v_gene) %>%
        tibble::column_to_rownames("v_gene")
      
      
      pdf(file.path(save_dir_stats, paste0(x, "_fraction_heatmap.pdf")),
          width = 3, height = 12)
      pheatmap(odds_heatmap_data, cluster_cols = FALSE, cluster_rows = FALSE,
               color = viridis(n = 20))
      dev.off()
      
      graphics.off()
      
      # Save all data to excel files
      # save_dir_stats_files
      save_file <- openxlsx::createWorkbook()
      
      # Save data used to make the heatmap
      openxlsx::addWorksheet(wb = save_file, sheetName = "heatmap_fraction")
      openxlsx::writeData(wb = save_file, sheet = "heatmap_fraction",
                          x = all_stats[[x]]$all_or_plot)
      
      # Save data with odds ratios and fisher exact p-value
      openxlsx::addWorksheet(wb = save_file, sheetName = "odds_ratio")
      openxlsx::writeData(wb = save_file, sheet = "odds_ratio",
                          x = odds_data)
      
      # Save data used to make the box plot
      openxlsx::addWorksheet(wb = save_file, sheetName = "boxplot_data")
      openxlsx::writeData(wb = save_file, sheet = "boxplot_data",
                          x = all_t_test_plot)
      
      # Save data with t test values
      save_t <- all_t %>%
        dplyr::filter(test == x)
      openxlsx::addWorksheet(wb = save_file, sheetName = "t_test")
      openxlsx::writeData(wb = save_file, sheet = "t_test",
                          x = save_t)
      
      openxlsx::saveWorkbook(wb = save_file,
                             file = file.path(save_dir_stats_files,
                                              paste0(x, ".xlsx")),
                             overwrite = TRUE)
      
    })
  })
  
  
  
  # Polar plots ------------------------------------------------------------------
  make_polar_plot <- function(starting_data, save_name,
                              save_dir){
    v_data <- starting_data %>%
      dplyr::filter(chains %in% c("IGH")) %>% 
      dplyr::filter(all_chains %in% keep_chains) %>%
      dplyr::select(sample, v_gene, Status) %>%
      dplyr::group_by(sample, v_gene) %>%
      dplyr::add_count(name = "sample_v_count") %>%
      dplyr::ungroup() %>%
      dplyr::group_by(sample) %>%
      dplyr::add_count(name = "sample_count") %>%
      dplyr::ungroup() %>%
      dplyr::group_by(v_gene) %>%
      dplyr::add_count(name = "full_v_count") %>%
      dplyr::ungroup() %>%
      dplyr::group_by(Status, v_gene) %>%
      dplyr::add_count(name = "status_v_count") %>%
      dplyr::ungroup() %>%
      dplyr::group_by(Status) %>%
      dplyr::add_count(name = "status_count") %>%
      dplyr::ungroup() %>%
      dplyr::distinct() %>%
      dplyr::mutate(frequency_v_within_sample = sample_v_count / sample_count * 100, 
                    frequency_v_by_sample = sample_v_count / full_v_count * 100,
                    frequency_v_within_status = status_v_count / status_count * 100,
                    frequency_v_by_status = status_v_count / full_v_count * 100) %>%
      dplyr::group_by(v_gene) %>%
      dplyr::add_count(name = "total_samples") %>%
      dplyr::ungroup() %>%
      dplyr::group_by(v_gene, Status) %>%
      dplyr::add_count(name = "samples_per_condition")
    
    v_data_plot <- v_data %>%
      dplyr::select(v_gene, frequency_v_within_sample,
                    sample, Status)
    
    v_order <- unique(v_data_plot$v_gene)
    
    fig <- plot_ly(
      type = 'scatterpolar',
      mode = 'lines'
    )
    
    
    for(i in unique(v_data_plot$sample)){
      plot_data <- v_data_plot %>%
        dplyr::filter(sample == i)
      plot_data <- plot_data[order(match(plot_data$v_gene, v_order)),]
      theta_vals <- c(plot_data$v_gene, plot_data$v_gene[[1]])
      fig <- fig %>%
        add_trace(
          r = c(plot_data$frequency_v_within_sample, 
                plot_data$frequency_v_within_sample[[1]]),
          theta = theta_vals,
          name = i,
          line=list(color = sample_colors[i])
        )
      
    }
    
    htmlwidgets::saveWidget(widget = fig, 
                            file.path(save_dir, 
                                      paste0(save_name,
                                             "_polar_plot_sample.html")))
    
    fig <- plot_ly(
      type = 'scatterpolar',
      mode = 'lines'
    )
    
    for(i in unique(v_data_plot$sample)){
      plot_data <- v_data_plot %>%
        dplyr::filter(sample == i)
      status <- unique(plot_data$Status)
      plot_data <- plot_data[order(match(plot_data$v_gene, v_order)),]
      theta_vals <- c(plot_data$v_gene, plot_data$v_gene[[1]])
      fig <- fig %>%
        add_trace(
          r = c(plot_data$frequency_v_within_sample, 
                plot_data$frequency_v_within_sample[[1]]),
          theta = theta_vals,
          name = i,
          line=list(color = status_colors[status])
        )
      
    }
    
    htmlwidgets::saveWidget(widget = fig, 
                            file.path(save_dir, 
                                      paste0(save_name,
                                             "_polar_plot_status_color.html")))
    
    # By status
    v_data_plot <- v_data %>%
      dplyr::select(v_gene, frequency_v_within_status,
                    Status) %>%
      dplyr::distinct()
    
    v_order <- unique(v_data_plot$v_gene)
    
    fig <- plot_ly(
      type = 'scatterpolar',
      mode = 'lines'
    )
    
    
    for(i in unique(v_data_plot$Status)){
      plot_data <- v_data_plot %>%
        dplyr::filter(Status == i)
      plot_data <- plot_data[order(match(plot_data$v_gene, v_order)),]
      theta_vals <- c(plot_data$v_gene, plot_data$v_gene[[1]])
      fig <- fig %>%
        add_trace(
          r = c(plot_data$frequency_v_within_status, 
                plot_data$frequency_v_within_status[[1]]),
          theta = theta_vals,
          name = i,
          line=list(color = status_colors[i])
        )
      
    }
    
    htmlwidgets::saveWidget(widget = fig, 
                            file.path(save_dir, 
                                      paste0(save_name,
                                             "_polar_plot_status.html")))
    
    return(v_data)
    
    
  }
  
  polar_dir <- file.path(vdj_dir, "polar_plots")
  ifelse(!dir.exists(polar_dir), dir.create(polar_dir), FALSE)
  
  
  all_b_data <- make_polar_plot(starting_data = all_info_split,
                                save_name = "all_b",
                                save_dir = polar_dir)
  
  # Now do it for tetramer and isotype
  tetramers_use <- c("DNA-tet", "Doublet", "GAD-tet", "IA2-tet",
                     "INS-tet", "Negative", "TET-tet")
  tet_data <- lapply(unique(tetramers_use), function(x){
    starting_data <- all_info_split %>%
      dplyr::filter(tet_hash_id == x)
    
    make_polar_plot(starting_data = starting_data,
                    save_name = x,
                    save_dir = polar_dir)
    
  })
  names(tet_data) <- tetramers_use
  
  isotypes_use <- c("IGHA", "IGHD", "IGHG", "IGHM")
  isotype_data <- lapply(unique(isotypes_use), function(x){
    starting_data <- all_info_split %>%
      dplyr::filter(isotype == x)
    
    make_polar_plot(starting_data = starting_data,
                    save_name = x,
                    save_dir = polar_dir)
    
  })
  
  names(isotype_data) <- isotypes_use
  
  polar_plots <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(wb = polar_plots, sheetName = "all_b")
  openxlsx::writeData(wb = polar_plots, sheet = "all_b", x = all_b_data)
  
  invisible(lapply(names(tet_data), function(x){
    openxlsx::addWorksheet(wb = polar_plots, sheetName = x)
    openxlsx::writeData(wb = polar_plots, sheet = x, x = tet_data[[x]])
    
  }))
  
  invisible(lapply(names(isotype_data), function(x){
    openxlsx::addWorksheet(wb = polar_plots, sheetName = x)
    openxlsx::writeData(wb = polar_plots, sheet = x, x = isotype_data[[x]])
    
  }))
  
  openxlsx::saveWorkbook(wb = polar_plots, 
                         file = file.path(vdj_files, "polar_plot_data.xlsx"),
                         overwrite = TRUE)
  
}
