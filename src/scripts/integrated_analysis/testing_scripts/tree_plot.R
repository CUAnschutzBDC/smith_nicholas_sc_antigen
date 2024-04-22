library(treemapify)
library(here)
library(scAnalysisR)
library(tidyverse)
library(Seurat)


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

# Set directories
save_dir <- file.path(results_dir, "R_analysis", sample)

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))

sep_columns <- c("chains", "cdr3", "cdr3_length",
                 "cdr3_nt_length", "v_gene", "d_gene", "j_gene", "c_gene",
                 "reads", "umis", "productive", "full_length",
                 "v_ins", "v_del", "v_mis", "d_ins", "d_del",
                 "d_mis", "j_ins", "j_del", "j_mis", "c_ins", "c_del", "c_mis",
                 "all_ins", "all_del", "all_mis", "vd_ins", "vd_del", "dj_ins",
                 "dj_del", "v_mis_freq", "d_mis_freq", "j_mis_freq",
                 "c_mis_freq", "all_mis_freq")
keep_columns <- c("isotype", "final_celltype", "sample", "paired",
                  "clonotype_id", "Status", 
                  "all_chains")

all_info <- seurat_data[[]] %>%
  dplyr::mutate(all_chains = chains) %>%
  dplyr::select(dplyr::all_of(c(sep_columns, keep_columns)), n_chains) %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::filter(!is.na(n_chains))


all_info_split <- cSplit(all_info, sep_columns, sep = ";", direction = "long") %>%
  dplyr::filter(!is.na(chains))


all_data <- all_info_split[all_info_split$chains == "IGH", c("v_gene", "Status")]

T1D <- all_data %>%
  dplyr::filter(Status == "T1D") %>%
  dplyr::group_by(v_gene) %>%
  dplyr::add_count(name = "v_count") %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  dplyr::mutate(v_percent = v_count / sum(v_count) * 100) %>%
  dplyr::mutate(plot_v = ifelse(v_percent > 0.25, v_gene, ""))


ggplot(T1D, aes(area = v_percent, fill = v_percent, label = plot_v)) +
  geom_treemap() +
  geom_treemap_text(colour = "white", place = "center", reflow = T)

