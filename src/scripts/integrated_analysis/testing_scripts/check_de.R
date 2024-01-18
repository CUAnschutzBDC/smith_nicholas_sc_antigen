library(tidyverse)
library(here)


source(here("src/scripts/muscat_plotting_functions.R"))

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "merged"

results_dir <- here("results")

# Set directories
save_dir <- file.path(results_dir, "R_analysis", sample)

processed_date <- read.csv(file.path(save_dir, "files",
                             "pseudobulk_all_bcell_de_processed_date.csv"))

# Save the results
collected_date_bin <- read.csv(file.path(save_dir, "files",
                             "pseudobulk_all_bcell_de.csv"))

capture_date <- read.csv(file.path(save_dir, "files",
                                   "pseudobulk_all_bcell_de_capture_correct.csv"))

no_ts_capture_date <- read.csv(file.path(save_dir, "files",
                                   "no_ts_pseudobulk_all_bcell_de_capture_correct.csv"))

table(no_ts_capture_date$contrast)

previous_res <- read.csv("/beevol/home/wellskri/Analysis/Mia_Smith/Catherine_Nicolas/20230224_BND_HC_NPOD_pLN_spl_3_NO_T1D_1_FDR/results/Catherine_combined/R_analysis/files/pseudobulk_all_bcell_de.csv")

processed_date_nd_no <- processed_date %>%
  dplyr::filter(contrast == "NO - ND")


collection_date_nd_no <- collected_date_bin %>%
  dplyr::filter(contrast == "NO - ND")

previous_res_nd_no <- previous_res %>%
  dplyr::filter(contrast == "NO - ND")

capture_res_nd_no <- capture_date %>%
  dplyr::filter(contrast == "NO - ND")

no_ts_res_nd_no <- no_ts_capture_date %>%
  dplyr::filter(contrast == "NO - ND")

dim(processed_date_nd_no)
dim(collection_date_nd_no)
dim(previous_res_nd_no)
dim(capture_res_nd_no)
dim(no_ts_res_nd_no)
length(intersect(processed_date_nd_no$gene, previous_res_nd_no$gene)) / nrow(processed_date_nd_no)
length(intersect(collection_date_nd_no$gene, previous_res_nd_no$gene)) / nrow(collection_date_nd_no)
length(intersect(capture_res_nd_no$gene, previous_res_nd_no$gene)) / nrow(capture_res_nd_no)
length(intersect(no_ts_res_nd_no$gene, previous_res_nd_no$gene)) / nrow(no_ts_res_nd_no)

# When doing the correction:
# The ND/NO is almost identical (less than 400/2500 are different)
# BUT the bigger problem is that we go from almost no changes between
# the subsets to many changes.

# SO, I fully trust the NO - ND comparison, this seems very good.
# I am more worried about the subsets.

# The two groups are:
# MV, NG, FG, and TS