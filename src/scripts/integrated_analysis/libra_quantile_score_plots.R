library(tidyverse)
library(pheatmap)
library(here)
library(viridis)


sample <- "merged"
results_dir <- here("results")
save_dir <- file.path(results_dir, "R_analysis", sample)
fig_dir <- file.path(save_dir, "images", "second_revision")

ifelse(!dir.exists(fig_dir), dir.create(fig_dir, recursive = TRUE), FALSE)

# Create custom color mapping function
create_custom_colors <- function(breaks, viridis_option = "viridis") {
  # Get number of breaks above 1
  n_colors_viridis <- sum(breaks > 1)
  
  # Create viridis colors for values > 1
  viridis_colors <- viridis(n_colors_viridis, option = viridis_option)
  
  # Create black colors for values <= 1
  black_colors <- rep("black", sum(breaks <= 1))
  
  # Combine colors
  custom_colors <- c(black_colors, viridis_colors)
  
  return(custom_colors)
}

clone_data <- file.path(save_dir, "files", "v_gene_counting", 
                        "clone_expansion.xlsx")

picked_clones <- openxlsx::readWorkbook(here("files/rAb_table_for_paper.xlsx"))

sheet_names <- openxlsx::getSheetNames(clone_data)

all_clones <- lapply(sheet_names, function(x){
  return(openxlsx::readWorkbook(clone_data, sheet = x))
  
})

all_clones <- do.call(rbind, all_clones)

keep_ids <- c(34204, 45125, 65093, 9070, 45089, 107462, 121787, 23895,
              100034, 29729, 57497, 111982, 31016, 4978, 66402, 5283, 
              9018, 33726, 7604, 111982, 16965, 110322, 29397, 40307, 
              56460, 7631, 111021)


keep_clones <- all_clones[all_clones$clone_id %in% keep_ids,]

short_clone_info <- lapply(1:nrow(picked_clones), function(x){
  print(x)
  line <- picked_clones[x,]
  clone_id <- line$`IgG.name.-.optional`
  clone_id <- gsub(".*_", "", clone_id)
  sequence_one <- line$Heavy.Chain.untrimmed
  sequence_two <- line$Light.Chain.untrimmed
  if(clone_id == "9018"){
    sequence_one <- gsub("TTTCTTATATGGGG", "", sequence_one)
  }
  new_sequence <- paste(sequence_one, sequence_two, sep = ";")
  return_lines <- keep_clones[keep_clones$clone_id == clone_id &
                              keep_clones$sequences == new_sequence,] %>%
    dplyr::distinct()
  
  if(nrow(return_lines > 0)){
    return_lines$CN_call <- line$CN.call
    
  }
  
  return(return_lines)
})

short_clone_info <- do.call(rbind, short_clone_info)

keep_barcodes <- c(
  "JH313-15_FD_CGGACACTCTGCGGCA-1",
  "JH313-15_DB_TATTACCGTCCTAGCG-1",
  "JH313-15_DB_CAGCTGGTCTCCTATA-1",
  "JH313-15_AG_TTCTTAGCAGATGAGC-1",
  "JH310-12_AP_CTTAGGAGTAGCCTAT-1", # OR JH310-12_AP_TTTCCTCGTTCGTTGA-1
  "JH313-15_FD_GCAATCAGTATGAAAC-1",
  "JH313-15_DB_CGAGAAGAGAAGCCCA-1",
  "JH313-15_AG_AAACGGGCATGGTAGG-1",
  "JH313-15_DB_AACCATGGTTCCACGG-1", # OR JH313-15_DB_GGACGTCCACGCGAAA-1
  "JH313-15_DB_GATCGTAAGCAGCGTA-1",
  "JH313-15_JB_CGCTGGATCCGTCAAA-1",
  "JH313-15_JB_GCTGCGAGTCTAGCCG-1",
  "JH310-12_BH_AACTCAGTCAACTCTT-1",
  "JH313-15_NF_AGGGTGAGTCTACCTC-1",
  "JH313-15_AG_CAAGGCCGTCGGCATC-1",
  "JH310-12_EB_CAGAATCAGGGATACC-1",
  "JH313-15_JB_CATCAAGCATTCGACA-1",
  "JH313-15_PH_CATATGGAGTATGACA-1",
  "JH310-12_BH_CTCTGGTAGGCTACGA-1", # OR JH310-12_BH_GGACATTCATACCATG-1
  "JH313-15_AG_CTCCTAGAGGTGCAAC-1", # OR JH313-15_AG_GTCTCGTCAAACAACA-1
  "JH313-15_AG_GGACGTCAGGCTAGCA-1", # OR JH313-15_AG_TGTATTCGTTCCAACA-1
  "JH310-12_CP_ACGGGCTTCAAGAAGT-1", # OR JH310-12_CP_GCACTCTCAGGCAGTA-1 OR JH310-12_CP_TGACTAGGTGAAAGAG-1
  "JH313-15_DB_CAGCTGGCACGCCAGT-1",
  "JH310-12_CP_TACTTACCAGGAACGT-1",
  "JH313-15_PH_ACGATGTAGTACGACG-1" # OR JH313-15_PH_GAAACTCGTGTTGGGA-1 OR JH313-15_PH_TTGCCGTAGATATGCA-1
)

# keep_barcodes <- c(
#   "JH313-15_FD_CGGACACTCTGCGGCA-1",
#   "JH313-15_DB_TATTACCGTCCTAGCG-1",
#   "JH313-15_DB_CAGCTGGTCTCCTATA-1",
#   "JH313-15_AG_TTCTTAGCAGATGAGC-1",
#   "JH310-12_AP_CTTAGGAGTAGCCTAT-1", "JH310-12_AP_TTTCCTCGTTCGTTGA-1",
#   "JH313-15_FD_GCAATCAGTATGAAAC-1",
#   "JH313-15_DB_CGAGAAGAGAAGCCCA-1",
#   "JH313-15_AG_AAACGGGCATGGTAGG-1",
#   "JH313-15_DB_AACCATGGTTCCACGG-1", "JH313-15_DB_GGACGTCCACGCGAAA-1",
#   "JH313-15_DB_GATCGTAAGCAGCGTA-1",
#   "JH313-15_JB_CGCTGGATCCGTCAAA-1",
#   "JH313-15_JB_GCTGCGAGTCTAGCCG-1",
#   "JH310-12_BH_AACTCAGTCAACTCTT-1",
#   "JH313-15_NF_AGGGTGAGTCTACCTC-1",
#   "JH313-15_AG_CAAGGCCGTCGGCATC-1",
#   "JH310-12_EB_CAGAATCAGGGATACC-1",
#   "JH313-15_JB_CATCAAGCATTCGACA-1",
#   "JH313-15_PH_CATATGGAGTATGACA-1",
#   "JH310-12_BH_CTCTGGTAGGCTACGA-1", "JH310-12_BH_GGACATTCATACCATG-1",
#   "JH313-15_AG_CTCCTAGAGGTGCAAC-1", "JH313-15_AG_GTCTCGTCAAACAACA-1",
#   "JH313-15_AG_GGACGTCAGGCTAGCA-1", "JH313-15_AG_TGTATTCGTTCCAACA-1",
#   "JH310-12_CP_ACGGGCTTCAAGAAGT-1", "JH310-12_CP_GCACTCTCAGGCAGTA-1", "JH310-12_CP_TGACTAGGTGAAAGAG-1",
#   "JH313-15_DB_CAGCTGGCACGCCAGT-1",
#   "JH310-12_CP_TACTTACCAGGAACGT-1",
#   "JH313-15_PH_ACGATGTAGTACGACG-1", "JH313-15_PH_GAAACTCGTGTTGGGA-1", "JH313-15_PH_TTGCCGTAGATATGCA-1"
# )


# Notes
# 31016 doesn't agree
# 5283 doesn't agree

final_clone_info <- short_clone_info %>%
  dplyr::filter(barcode %in% keep_barcodes) %>%
  dplyr::select(clone_id, full_tet_name_cutoff, scar_libra_full_hash_id,
                CN_call, dplyr::contains("myeloid_t"),
                dplyr::contains("libra_cutoff"))

final_clone_info$clone_id <- make.unique(as.character(final_clone_info$clone_id))

quantile_normalized <- final_clone_info %>%
  tibble::column_to_rownames("clone_id") %>%
  dplyr::select(dplyr::contains("myeloid_t"))

colnames(quantile_normalized) <- gsub("myeloid_t_cutoff_", "",
                                      colnames(quantile_normalized))

column_order <- c("INS.tet", "IA2.tet", "GAD.tet", "TET.tet", "DNA.tet")
row_order <- c(45089, 45125, 65093, 9070,  100034, 107462, 121787, 
               23895, 29729, 57497, 111982, 31016, 4978, 66402,
               9018, 5283, 106234, 33726, 7604, 111982, 16965, 110322,
               29397, 40307, 56450, 7631, 111026)

quantile_normalized <- quantile_normalized[order(match(rownames(quantile_normalized),
                                                       row_order)),
                                           order(match(colnames(quantile_normalized),
                                                       column_order))]

breaks <- seq(0, 3, length.out = 100)  # Adjust number of breaks as needed
colors <- create_custom_colors(breaks)

graphics.off()
pdf(file.path(fig_dir, "quantile_normalized_scores.pdf"))
pheatmap(quantile_normalized, cluster_rows = FALSE, cluster_cols = FALSE,
         color = colors, breaks = breaks)
dev.off()



# Libra normalized
libra_normalized <- final_clone_info %>%
  tibble::column_to_rownames("clone_id") %>%
  dplyr::select(dplyr::contains("libra_cutoff"))

colnames(libra_normalized) <- gsub("libra_cutoff_", "",
                                      colnames(libra_normalized))

column_order <- c("INS.tet", "IA2.tet", "GAD.tet", "TET.tet", "DNA.tet")
row_order <- c(45089, 45125, 65093, 9070,  100034, 107462, 121787, 
               23895, 29729, 57497, 111982, 31016, 4978, 66402,
               9018, 5283, 106234, 33726, 7604, 111982, 16965, 110322,
               29397, 40307, 56450, 7631, 111026)

libra_normalized <- libra_normalized[order(match(rownames(libra_normalized),
                                                       row_order)),
                                           order(match(colnames(libra_normalized),
                                                       column_order))]

breaks <- seq(0, 3, length.out = 100)  # Adjust number of breaks as needed
colors <- create_custom_colors(breaks)


graphics.off()
pdf(file.path(fig_dir, "libra_normalized_scores.pdf"))

pheatmap(libra_normalized, cluster_rows = FALSE, cluster_cols = FALSE,
         color = colors, breaks = breaks)

dev.off()
