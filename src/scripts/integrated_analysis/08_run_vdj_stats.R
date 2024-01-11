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

# Colors -----------------------------------------------------------------------
all_colors <- readRDS(file = file.path("files/all_colors.rds"))


final_colors <- all_colors$cell_type_colors

tetramer_colors <- all_colors$tetramer_colors

sample_colors <- all_colors$sample_colors


status_colors <- all_colors$status_colors


for(cells_use in c("all", "memory")){
  
  print(cells_use)
  
  # Read in data
  seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))
  
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

  
  # Analysis ---------------------------------------------------------------------
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
  
  # Make a data frame
  # Need
  # 1. Percent of the vj gene (Here, adding across samples should be 100%)
  # 2. Percent of the sample (Here adding across vj genes should be 100%)
  # 3. Number of samples
  # 4. Number of samples per status
  # Do for V and VJ

  keep_chains <- c("IGH", "IGH;IGK", "IGH;IGK;IGK",
                   "IGH;IGK;IGL", "IGH;IGL", 
                   "IGH;IGL;IGL")
  
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
  
  
  all_samples_t1d <- v_data %>%
    ungroup() %>%
    dplyr::filter(Status == "T1D")%>%
    dplyr::select(sample) %>%
    dplyr::distinct()
  
  all_samples_nd <- v_data %>%
    ungroup() %>%
    dplyr::filter(Status == "ND")%>%
    dplyr::select(sample) %>%
    dplyr::distinct()
  
  all_samples_aab <- v_data %>%
    ungroup() %>%
    dplyr::filter(Status == "AAB")%>%
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
  
  
  all_samples_t1d <- all_info_split %>%
    ungroup() %>%
    dplyr::filter(Status == "T1D")%>%
    dplyr::select(sample) %>%
    dplyr::distinct()
  
  all_samples_nd <- all_info_split %>%
    ungroup() %>%
    dplyr::filter(Status == "ND")%>%
    dplyr::select(sample) %>%
    dplyr::distinct()
  
  all_samples_aab <- all_info_split %>%
    ungroup() %>%
    dplyr::filter(Status == "AAB")%>%
    dplyr::select(sample) %>%
    dplyr::distinct()
  
  
  all_tests <- list("t1d_nd" = c("T1D", "ND"),
                    "aab_nd" = c("AAB", "ND"))
  
  sample_mapping <- list("ND" = all_samples_nd,
                         "T1D" = all_samples_t1d,
                         "AAB" = all_samples_aab)
  
  build_tests <- lapply(c("all", "antigen", "isotype"), function(x){
    if(x == "all"){
      v_df <- list(comparison_builder(starting_df = all_info_split))
      names(v_df) <- "all"
    } else if(x == "antigen") {
      v_df <- lapply(unique(all_info_split$scar_hash_id), function(x){
        subset_df <- all_info_split %>%
          dplyr::filter(scar_hash_id == x)
        return(comparison_builder(starting_df = subset_df))
      })
      names(v_df) <- unique(all_info_split$scar_hash_id)
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
    all_test_full <- lapply(names(v_df), function(test_type){
      individual_test <- lapply(all_tests, function(y){
        return_data <- list(y, v_df[[test_type]])
      })
      names(individual_test) <- paste(test_type, names(individual_test), sep = "_")
      return(individual_test)
    })
    
    all_test_full <- do.call("c", all_test_full)
    
    all_stats <- lapply(all_test_full, function(test_use){
      status_one <- test_use[[1]][[1]]
      status_two <- test_use[[1]][[2]]
      all_samples_one <- sample_mapping[[status_one]]
      all_samples_two <- sample_mapping[[status_two]]
      v_j_statistics(test_use[[2]], status_one = status_one,
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
    
    all_plots <- lapply(names(all_test_full), function(y){
      all_t_test_plot <- all_stats[[y]]$all_t_test_plot %>%
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
      pdf(file.path(save_dir_stats, paste0(y, "_t_test.pdf")),
          width = 12, height = 8)
      print(t_plot)
      dev.off()
      
      odds_data <- all_odds %>%
        dplyr::filter(test == y)
      
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
      pdf(file.path(save_dir_stats, paste0(y, "_odds_ratio.pdf")),
          width = 8, height = 12)
      print(odds_plot)
      dev.off()
      
      odds_heatmap_data <- all_stats[[y]]$all_or_plot %>% 
        dplyr::select(v_gene, Status, fraction) %>%
        tidyr::pivot_wider(names_from = Status, values_from = fraction) %>%
        dplyr::arrange(v_gene) %>%
        tibble::column_to_rownames("v_gene")
      
      
      pdf(file.path(save_dir_stats, paste0(y, "_fraction_heatmap.pdf")),
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
                          x = all_stats[[y]]$all_or_plot)
      
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
        dplyr::filter(test == y)
      openxlsx::addWorksheet(wb = save_file, sheetName = "t_test")
      openxlsx::writeData(wb = save_file, sheet = "t_test",
                          x = save_t)
      
      openxlsx::saveWorkbook(wb = save_file,
                             file = file.path(save_dir_stats_files,
                                              paste0(y, ".xlsx")),
                             overwrite = TRUE)
      
    })
  })
}
  
