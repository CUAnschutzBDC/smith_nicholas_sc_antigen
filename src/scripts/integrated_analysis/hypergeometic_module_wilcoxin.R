library(here)
library(scAnalysisR)
library(tidyverse)
library(Seurat)

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

gene_lists <- readRDS(here("files/human_gene_lists.rda"))

gene_lists <- lapply(gene_lists, function(x){
  return(x$V1)
})

seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))
all_human_chemokines <- read.csv(here("files/all_human_chemokines.csv"))

gene_lists$chemokines <- all_human_chemokines$Systemic.name

gene_lists_use <- gene_lists[c("antigen_processing", "chemokines", 
                                "costimulation")]

de_table <- read.csv(file.path(save_dir, "files", 
                               "pseudobulk_all_bcell_de_processed_date.csv"))

# Colors -----------------------------------------------------------------------


status_colors <- MetBrewer::met.brewer(name = "Egypt", n = 4,
                                       type = "discrete")

names(status_colors) <- c("nd", "no", "aab_stage_1", "aab_stage_2")

# Hypergeometric ---------------------------------------------------------------

# Columns need cluster, p_val_adj, avg_log2FC

logfc_up <- de_table %>%
  dplyr::filter(logFC > 0) %>% # If the logfc is positive, we want the first from the contrast
  dplyr::mutate(cluster = gsub(" - .*", "", contrast))

logfc_down <- de_table %>%
  dplyr::filter(logFC < 0) %>% # If the logfc is positive, we want the first from the contrast
  dplyr::mutate(cluster = gsub(".* - ", "", contrast)) %>%
  dplyr::mutate(logFC = abs(logFC))

de_df <- rbind(logfc_up, logfc_down) %>%
  dplyr::rename(p_val_adj = p_adj.glb, avg_log2FC = logFC)

# Run hypergeometric
hypergeometric <- hypergeometric_test(seurat_object = seurat_data,
                                      gene_list = gene_lists_use,
                                      DE_table = de_df,
                                      DE_p_cutoff = 0.05,
                                      DE_lfc_cutoff = 0.5,
                                      correction_method = "fdr")

sample_info <- data.frame(cluster = unique(hypergeometric$cluster))
rownames(sample_info) <- sample_info$cluster


heatmap1 <- plot_hypergeom(hypergeometric,
                           meta_df = sample_info,
                           max_val = 10)


# Module analysis --------------------------------------------------------------
test_module <- function(seurat_object, test_var, group_var, features,
                        type = "two_tailed",
                        pairwise = TRUE){
  if(group_var == "all"){
    seurat_object$all <- "all_cells"
  }
  test_info <- lapply(features, function(x){
    module <- FetchData(seurat_object, vars = c(x, test_var, group_var)) %>%
      dplyr::mutate(module_test = get(x), group_var = get(group_var),
                    test_var = get(test_var))
    test_info_one <- lapply(unique(module$group_var), function(y){
      module_short <- module %>%
        dplyr::filter(group_var == y)
      if (length(unique(module_short$test_var)) > 2 & pairwise &
          type == "two_tailed"){
        res <- pairwise.wilcox.test(module_short$module_test,
                                    module_short$test_var,
                                    p.adjust.method = "none")
        res_dat <- res$p.value
        all_cols <- colnames(res_dat)
        res_dat <- res_dat %>%
          data.frame() %>%
          tibble::rownames_to_column("group1") %>%
          tidyr::pivot_longer(all_of(all_cols), names_to = "group2",
                              values_to = "p_value",
                              values_drop_na = TRUE) %>%
          dplyr::mutate(module = x, cluster = y) %>%
          data.frame
        return(res_dat)
        
      } else if(length(unique(module_short$test_var)) > 2 & pairwise &
                type == "one_tailed"){
        res_g <- pairwise.wilcox.test(module_short$module_test,
                                      module_short$test_var,
                                      p.adjust.method = "none",
                                      alternative = "g")
        res_dat_g <- res_g$p.value
        all_cols <- colnames(res_dat_g)
        res_dat_g <- res_dat_g %>%
          data.frame() %>%
          tibble::rownames_to_column("group1") %>%
          tidyr::pivot_longer(all_of(all_cols), names_to = "group2",
                              values_to = "p_value_higher",
                              values_drop_na = TRUE)
        
        res_l <- pairwise.wilcox.test(module_short$module_test,
                                      module_short$test_var,
                                      p.adjust.method = "none",
                                      alternative = "l")
        res_dat_l <- res_l$p.value
        all_cols <- colnames(res_dat_l)
        res_dat_l <- res_dat_l %>%
          data.frame() %>%
          tibble::rownames_to_column("group1") %>%
          tidyr::pivot_longer(all_of(all_cols), names_to = "group2",
                              values_to = "p_value_lower",
                              values_drop_na = TRUE)
        
        res_dat <- left_join(res_dat_l, res_dat_g,
                             by = c("group1", "group2")) %>%
          dplyr::mutate(module = x, cluster = y) %>%
          data.frame
        
        return(res_dat)
      } else if(length(unique(module_short$test_var)) > 2){
        res <- kruskal.test(module_test ~ test_var, data = module_short)
        res_dat <- data.frame(module = x, p_value = res$p.value,
                              cluster = y)
        return(res_dat)
      } else if(type == "two_tailed"){
        res <- wilcox.test(module_test ~ test_var, data = module_short)
        res_dat <- data.frame(module = x, p_value = res$p.value, cluster = y)
        return(res_dat)
      } else if(type == "one_tailed"){
        module_short$test_var <- factor(module_short$test_var)
        res_g <- wilcox.test(module_test ~ test_var, data = module_short,
                             alternative = "g")
        res_l <- wilcox.test(module_test ~ test_var, data = module_short,
                             alternative = "l")
        res_dat <- data.frame(module = x, p_value_higher = res_g$p.value,
                              p_value_lower = res_l$p.value,
                              reference = levels(module_short$test_var)[2],
                              cluster = y)
        return(res_dat)
      }
    }) 
    test_info <- do.call(rbind, test_info_one)
    return(test_info)
    
  })
  test_info <- do.call(rbind, test_info)
  
  if(type == "two_tailed" | !(pairwise)){
    p_vals <- p.adjust(test_info$p_value, method = "BH")
    test_info$p_adj <- p_vals
  } else if(type == "one_tailed"){
    p_vals <- test_info %>%
      pivot_longer(cols = c(p_value_higher, p_value_lower), names_to = "p_val_type",
                   values_to = "p_val")
    
    p_vals$p_adj <- p.adjust(p_vals$p_val, method = "BH")
    
    test_info <- p_vals %>%
      mutate(p_val_type = sub("p_value_", "", p_val_type)) %>%
      pivot_wider(names_from = p_val_type, values_from = c(p_adj, p_val)) %>%
      data.frame()
  }
  return(test_info)
}

for(i in names(gene_lists_use)){
  print(i)
  seurat_data <- AddModuleScore(seurat_data, features = gene_lists_use[i],
                                name = paste0(i, "_"))
  
}


features <- colnames(seurat_data[[]])[grepl("_1$", colnames(seurat_data[[]]))]

# Two tailed
seurat_data$Status <- gsub(" ", "_", seurat_data$Status)

seurat_data$Status <- factor(seurat_data$Status,
                             levels = c("nd", "aab_stage_1", "aab_stage_2", "no"))

test_info_two_tailed <- test_module(seurat_object = seurat_data,
                                    test_var = "Status",
                                    group_var = "all",
                                    features = features,
                                    type = "two_tailed")

violin_plots <- featDistPlot(seurat_data,
                             geneset = paste0(names(gene_lists_use), "_1"),
                             combine = FALSE,
                             col_by = "Status",
                             sep = "tet_hash_id",
                             color = status_colors)

p_val_plot <- lapply(names(violin_plots), function(module_x){
  p_values <- test_info_two_tailed %>%
    dplyr::filter(module == module_x) %>%
    dplyr::select(p_adj, group1, group2) %>%
    dplyr::mutate(max_val = max(seurat_data[[module_x]]),
                  sep = c(0.2, 0.5, 0.35, 0.8, 0.65, 0.2)) %>%
    dplyr::mutate(y.position = max_val + max_val * sep) %>%
    dplyr::mutate(p_adj = signif(p_adj, 3))
  
  
  violin_plot <- violin_plots[[module_x]] +
    ggpubr::stat_pvalue_manual(p_values,
                               label = "p_adj")
  
  return(violin_plot)
})
