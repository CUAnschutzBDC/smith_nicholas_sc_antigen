library(tidyverse)
library(openxlsx)
library(here)
library(splitstackshape)

sample <- "merged"

results_dir <- here("results")
save_dir <- file.path(results_dir, "R_analysis", sample)

v_counting_dir <- file.path(save_dir, "files", "v_gene_counting")


clone_list <- read.csv(file.path(v_counting_dir, 
                                 "all_clone_info.csv"))

clone_info <- read.table(file.path(save_dir, "define_clones",
                                   "immcantation_combined_clone-pass.tsv"),
                         sep = "\t", header = TRUE)

# add c_call to the so

# These are found within samples individually
sample_clone_info <- read.table(file.path(save_dir, "define_clones",
                                          "immcantation_combined_clone-pass_sample.tsv"),
                                sep = "\t", header = TRUE)

sample_clone_info <- sample_clone_info %>%
  dplyr::select(sequence_id, clone_id, locus, sample) %>%
  dplyr::rename(sample_clone_id = clone_id,
                sample_locus = locus) %>%
  dplyr::mutate(sequence_id_sample = paste(sequence_id,
                                           sample, sep = "_")) %>%
  dplyr::select(-sequence_id, -sample)

clone_info <- clone_info %>%
  dplyr::mutate(sequence_id_sample = paste(sequence_id,
                                           sample, sep = "_")) %>%
  merge(sample_clone_info, by = "sequence_id_sample", all.x = TRUE) %>%
  dplyr::mutate(cell_sample = paste(sample, cell_id, sep = "_")) %>%
  dplyr::select(-sample)


catherine_data_path <- here("files/GenScript_rAB_sequences.xlsx")
catherine_data <- openxlsx::read.xlsx(catherine_data_path) 

catherine_data_path_two <- here("files/20240314_additional_13_rAb_sequences_to_trim.xlsx")
catherine_data_two <- openxlsx::read.xlsx(catherine_data_path_two)


sep_columns <- c("v_gene", "j_gene", "chains",
                 "sequences")
keep_columns <- c("clone_id")

all_info <- catherine_data_two %>%
  dplyr::select(dplyr::all_of(c(sep_columns, keep_columns))) 


all_info_split <- cSplit(all_info, sep_columns, sep = ";", direction = "long")

heavy_info <- all_info_split[all_info_split$chains == "IGH",]
light_info <- all_info_split[all_info_split$chains != "IGH", ]

# Remove any overlaps
`%notin%` <- Negate(`%in%`)
heavy_info <- heavy_info[heavy_info$sequences %notin% catherine_data$Heavy.Chain,]

clone_info_heavy <- clone_info[clone_info$sequence %in% 
                                 c(catherine_data$Heavy.Chain,
                                   heavy_info$sequences), ] %>%
  dplyr::select(-sequence_id_sample, -sequence_id, -cell_id, -consensus_count,
                -umi_count, -cell_sample) %>%
  distinct()

clone_info_heavy <- clone_info_heavy[order(match(clone_info_heavy$sequence,
                                                 c(catherine_data$Heavy.Chain,
                                                   heavy_info$sequences))),]

clone_info_heavy$clone_id <- c(catherine_data$clone_id, heavy_info$clone_id)

# Remove any overlaps
`%notin%` <- Negate(`%in%`)
light_info <- light_info[light_info$sequences %notin% catherine_data$Light.Chain,]

clone_info_light <- clone_info[clone_info$sequence %in% 
                                 c(catherine_data$Light.Chain,
                                   light_info$sequences), ] %>%
  dplyr::select(-sequence_id_sample, -sequence_id, -cell_id, -consensus_count,
                -umi_count, -cell_sample, -sample_clone_id, -clone_id) %>%
  distinct()


clone_info_light <- clone_info_light[order(match(clone_info_light$sequence,
                                                 c(catherine_data$Light.Chain,
                                                   light_info$sequences))),]


clone_info_light$clone_id <- c(catherine_data$clone_id, light_info$clone_id)

clone_info_light$sample_clone_id <- "NA"

full_clone_info <- rbind(clone_info_heavy, clone_info_light)

write.table(full_clone_info, file.path(v_counting_dir, "all_v_info.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
