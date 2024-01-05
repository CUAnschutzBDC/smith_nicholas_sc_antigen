library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(scuttle)
library(here)
library(openxlsx)
library(djvdj)


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "JH310-12_EB"

sample_info <- args[[4]]
#sample_info <- here("files/sample_info.tsv")

results_dir <- args[[2]]
#results_dir <- here("results")

sample_metadata <- args[[3]]
#sample_metadata <- here("files/Deidentified_donor_metadata.xlsx")
sample_metadata <- openxlsx::readWorkbook(sample_metadata, detectDates = TRUE)
colnames(sample_metadata) <- make.names(colnames(sample_metadata))
sample_metadata <- sample_metadata[!(is.na(sample_metadata$Sample.Name)),]
sample_metadata <- sample_metadata[sample_metadata$Sample.Name == sample,] %>%
  dplyr::select(!X.Key)

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

vars.to.regress <- NULL

# Set directories
save_dir <- file.path(results_dir, "R_analysis", sample)

# Make directories
ifelse(!dir.exists(save_dir), dir.create(save_dir), FALSE)

ifelse(!dir.exists(file.path(save_dir, "images")),
       dir.create(file.path(save_dir, "images")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "files")),
       dir.create(file.path(save_dir, "files")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
       dir.create(file.path(save_dir, "rda_obj")), FALSE)

mt_pattern <- "^MT-" # "^MT-" for human, "^mt-" for mice

# Create seurat object ---------------------------------------------------------
seurat_object <- create_seurat_object(sample = sample,
                                      count_path = results_dir,
                                      ADT = ADT, hashtag = HTO,
                                      tenx_structure = "multi7",
                                      min_features = 200
                                      )

# Add mitochondrial percent
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,
                                                      pattern = mt_pattern)

sample_vect <- sample_metadata[1,,drop = TRUE]

seurat_object <- AddMetaData(seurat_object, metadata = sample_vect)

# Use scuttle for cutoffs ------------------------------------------------------
se <- as.SingleCellExperiment(seurat_object)
is.mito <- grep(mt_pattern, rownames(se))
per.cell <- perCellQCMetrics(se, subsets=list(Mito=is.mito))

qc.stats <- perCellQCFilters(per.cell,
                             sum.field = "sum",
                             detected.field = "detected",
                             sub.fields = "subsets_Mito_percent")

all_info <- cbind(per.cell, qc.stats)
colnames(all_info) <- paste0("cell_qc_", colnames(all_info))
all_info <- data.frame(all_info)
seurat_object <- AddMetaData(seurat_object, metadata = all_info)

rna_qual <- featDistPlot(seurat_object,
                          geneset = c("nFeature_RNA", "nCount_RNA",
                                      "percent.mt"),
                          plot_type = "violin", col_by = "cell_qc_discard")

if(ADT){
  adt_qual <- featDistPlot(seurat_object,
                            geneset = c("nFeature_ADT", "nCount_ADT",
                                        "percent.mt"),
                            plot_type = "violin", col_by = "cell_qc_discard")
}

pdf(file.path(save_dir, "images", "quality_plots.pdf"))
plot(rna_qual)
if(ADT){
  plot(adt_qual)
}
dev.off()

# Save before moving on
saveRDS(seurat_object, file = file.path(save_dir, "rda_obj",
                                        "seurat_unfilt.rds"))


seurat_object <- subset(seurat_object, subset = cell_qc_discard == FALSE)

# Normalization
# Single Cell Transform normalization
seurat_object <- SCTransform(seurat_object, vars.to.regress = vars.to.regress,
                             verbose = FALSE)

# Default normalization
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object) %>% 
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = vars.to.regress)

# Add in cell cycle
seurat_object <- CellCycleScoring(seurat_object, 
                                  g2m.features = cc.genes$g2m.genes,
                                  s.features = cc.genes$s.genes)


# Check against dropkick cells -------------------------------------------------
unfilt_path <- file.path(results_dir, sample, "outs", "multi", "count",
                         "raw_feature_bc_matrix")

full_seurat <- Read10X(unfilt_path)

full_seurat <- CreateSeuratObject(counts = full_seurat[["Gene Expression"]], 
                                  project = sample, min.cells = 0, 
                                  min.features = 0)

# Add mitochondrial percent
full_seurat[["percent.mt"]] <- PercentageFeatureSet(full_seurat,
                                                    pattern = mt_pattern)


dropkick_cells <- read.csv(file.path(save_dir, "files", "dropkick_cells.csv"),
                           row.names = 1)

full_seurat <- AddMetaData(full_seurat, dropkick_cells)

scuttle_meta <- seurat_object[[]] %>%
  dplyr::select(cell_qc_low_n_features, cell_qc_high_subsets_Mito_percent,
                cell_qc_low_lib_size) %>%
  dplyr::mutate(cell_qc_cell = "cell_qc_true")

full_seurat <- AddMetaData(full_seurat, scuttle_meta)

full_seurat$dropkick_label[is.na(full_seurat$dropkick_label)] <- "dropkick_false"
full_seurat$dropkick_label[full_seurat$dropkick_label == "True"] <- "dropkick_true"


full_seurat$cell_qc_cell[is.na(full_seurat$cell_qc_cell)] <- "cell_qc_false"

cm <- confusionMatrix(full_seurat$dropkick_label, full_seurat$cell_qc_cell)
cm <- cm / rowSums(cm)

comparison_heatmap <- pheatmap::pheatmap(cm)

full_seurat$cell_labels <- ifelse(full_seurat$dropkick_label == "dropkick_true" 
                                  & full_seurat$cell_qc_cell == "cell_qc_true",
                                  "both_cell", 
                                  ifelse(full_seurat$dropkick_label == "dropkick_true",
                                         "dropkick_only",
                                         ifelse(full_seurat$cell_qc_cell == "cell_qc_true",
                                                "cell_qc_cell", "no_cell")))


barplot_one <- scAnalysisR::stacked_barplots(full_seurat, 
                                             meta_col = "cell_labels")

rna_qual_two <- featDistPlot(full_seurat,
                             geneset = c("nFeature_RNA", "nCount_RNA",
                                         "percent.mt"),
                             plot_type = "violin", col_by = "cell_labels")

full_seurat <- subset(full_seurat, subset = cell_labels != "no_cell")

barplot_two <- scAnalysisR::stacked_barplots(full_seurat, 
                                             meta_col = "cell_labels")

graphics.off()

pdf(file.path(save_dir, "images", "dropkick_vs_cellqc.pdf"))

print(barplot_one)
plot(rna_qual_two)
print(barplot_two)
grid::grid.newpage()
print(comparison_heatmap)

dev.off()

rm(full_seurat)

# VDJ --------------------------------------------------------------------------
# Read in VDJ data
vdj_dir <- file.path(results_dir, sample, "outs/per_sample_outs",
                     sample, "vdj_b")

seurat_object <- import_vdj(input = seurat_object,
                            vdj_dir = vdj_dir,
                            filter_paired = FALSE,
                            include_mutations = TRUE)

# Tetramers --------------------------------------------------------------------
# First pull out the tetramers
all_adts <- GetAssayData(seurat_object, assay = "ADT",
                         slot = "counts")

# First remove the incorrect insulin and other unnecessary ADTs
all_sums <- rowSums(all_adts)

# 1000 is arbitrary but should work to identify the actual ADTs. This should
# be checked and adjusted.
keep_adts <- all_sums[all_sums > 1000] 

all_adts <- all_adts[names(keep_adts),]

rownames(all_adts) <- gsub("INS-tet[a|b]", "INS-tet", rownames(all_adts))

seurat_object[["ADT"]] <- CreateAssayObject(counts = all_adts)

seurat_object <- NormalizeData(seurat_object, assay = "ADT",
                               normalization.method = "CLR", margin = 2)

all_tetramers <- all_adts[grepl("tet", rownames(all_adts)),]

seurat_object[["TET"]] <- CreateAssayObject(counts = all_tetramers)

seurat_object <- NormalizeData(seurat_object, assay = "TET",
                               normalization.method = "CLR", margin = 2)

# Add in scar counts -----------------------------------------------------------
scar_counts <- read.csv(file.path(save_dir,
                                  "files", "scar_denoised.csv"),
                        row.names = 1) %>%
  t()

scar_counts <- scar_counts[ , colnames(scar_counts) %in%
                              colnames(seurat_object)]

keep_scar <- gsub("-", "_", names(keep_adts))
scar_counts <- scar_counts[keep_scar,]

rownames(scar_counts) <- gsub("INS_tet[a|b]", "INS_tet", rownames(scar_counts))

seurat_object[["SCAR_ADT_LOG"]] <- CreateAssayObject(counts = scar_counts)

norm_data <- log1p(scar_counts)

seurat_object <- SetAssayData(object = seurat_object, slot = "data", 
                              assay = "SCAR_ADT_LOG", new.data = norm_data)

seurat_object[["SCAR_ADT"]] <- CreateAssayObject(counts = scar_counts)

seurat_object <- NormalizeData(seurat_object, assay = "SCAR_ADT",
                               normalization.method = "CLR",
                               margin = 2)


all_tetramers <- scar_counts[grepl("tet", rownames(scar_counts)),]

seurat_object[["SCAR_TET"]] <- CreateAssayObject(counts = all_tetramers)

seurat_object <- NormalizeData(seurat_object, assay = "SCAR_TET",
                               normalization.method = "CLR",
                               margin = 2)

seurat_object[["SCAR_TET_LOG"]] <- CreateAssayObject(counts = all_tetramers)

norm_data <- log1p(all_tetramers)

seurat_object <- SetAssayData(object = seurat_object, slot = "data", 
                              assay = "SCAR_TET_LOG", new.data = norm_data)


# HTO demux --------------------------------------------------------------------
HTODemuxUpdate <- function(object, assay = "HTO", positive.quantile = 0.99,
                           init = NULL, nstarts = 100, kfunc = "clara", 
                           nsamples = 100, seed = 42, verbose = TRUE){
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay)
  counts <- GetAssayData(object = object, assay = assay, slot = "counts")[, 
                                                                          colnames(x = object)]
  counts <- as.matrix(x = counts)
  ncenters <- init %||% (nrow(x = data) + 1)
  switch(EXPR = kfunc, kmeans = {
    init.clusters <- kmeans(x = t(x = GetAssayData(object = object, 
                                                   assay = assay)), centers = ncenters, nstart = nstarts)
    Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
  }, clara = {
    init.clusters <- clara(x = t(x = GetAssayData(object = object, 
                                                  assay = assay)), k = ncenters, samples = nsamples)
    Idents(object = object, cells = names(x = init.clusters$clustering), 
           drop = TRUE) <- init.clusters$clustering
  }, stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'"))
  average.expression <- AverageExpression(object = object, 
                                          assays = assay, verbose = FALSE)[[assay]]
  if (sum(average.expression == 0) > 0) {
    stop("Cells with zero counts exist as a cluster.")
  }
  discrete <- GetAssayData(object = object, assay = assay)
  discrete[discrete > 0] <- 0
  
  proportion <- discrete
  for (iter in rownames(x = data)) {
    values <- counts[iter, colnames(object)]
    values.use <- values[WhichCells(object = object, idents = levels(x = Idents(object = object))[[which.min(x = average.expression[iter, 
    ])]])]
    fit <- suppressWarnings(expr = fitdistrplus::fitdist(data = values.use, 
                                                         distr = "nbinom"))
    cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
    discrete[iter, names(x = which(x = values > cutoff))] <- 1
    if (verbose) {
      message(paste0("Cutoff for ", iter, " : ", cutoff, 
                     " reads"))
    }
    if(cutoff == 0){
      div_by <- 0.1
    } else {
      div_by <- cutoff
    }
    
    proportion[iter,] <- values / div_by
  }
  npositive <- colSums(x = discrete)
  classification.global <- npositive
  classification.global[npositive == 0] <- "Negative"
  classification.global[npositive == 1] <- "Singlet"
  classification.global[npositive > 1] <- "Doublet"
  
  # Here I do it based on the proportion (count / cutoff). This way, if
  # the value was higher than the cutoff it will be greater than 1 and if it
  # is lower it will be less than 1.
  donor.id <- rownames(x = proportion)
  hash.max <- apply(X = proportion, MARGIN = 2, FUN = max)
  hash.maxID <- sapply(X = 1:ncol(x = proportion), FUN = function(x) {
    max_indices <- which(proportion[, x] == hash.max[x])
    return(donor.id[min(max_indices)])
  })
  
  hash.second <- apply(X = proportion, MARGIN = 2, FUN = function(x) {
    sorted_values <- sort(x, decreasing = TRUE)
    if (length(sorted_values) > 1) {
      second_max <- sorted_values[2] # Get the second highest value
      return(second_max)
    } else {
      return(NA) # Handle cases where there is no second highest value
    }
  })
  
  hash.secondID <- sapply(X = 1:ncol(x = proportion), FUN = function(x) {
    idx <- which(proportion[, x] == hash.second[x])
    if (length(idx) > 1) {
      # If there are multiple indices for the second highest value, choose the one not equal to hash.maxID
      idx <- idx[which(names(idx) != hash.maxID[x])]
    }
    return(donor.id[idx[1]])
  })
  
  
  hash.margin <- hash.max - hash.second
  doublet_id <- sapply(X = 1:length(x = hash.maxID), FUN = function(x) {
    return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), 
                 collapse = "_"))
  })
  classification <- classification.global
  classification[classification.global == "Negative"] <- "Negative"
  classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == 
                                                                           "Singlet")]
  classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == 
                                                                           "Doublet")]
  classification.metadata <- data.frame(hash.maxID, hash.secondID, 
                                        hash.margin, classification, classification.global)
  colnames(x = classification.metadata) <- paste(assay, c("maxID", 
                                                          "secondID", "margin", "classification", "classification.global"), 
                                                 sep = "_")
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste0(assay, "_classification")
  doublets <- rownames(x = object[[]])[which(object[[paste0(assay, 
                                                            "_classification.global")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- "Doublet"
  object$hash.ID <- Idents(object = object)
  
  
  return(list(object = object, proportions = proportion))
}

find_hash_id <- function(proportions){
  new_hash_id <- apply(proportions, MARGIN = 2, FUN = function(x){
    positive_hits <- x[x > 1]
    if(length(positive_hits) == 0){
      return(data.frame("new_hash_id" = "Negative", "full_hash_id" = "Negative"))
    } else if (length(positive_hits) == 1){
      return(data.frame("new_hash_id" = names(positive_hits),
                        "full_hash_id" = names(positive_hits)))
    } else if(any(grepl("DNA-|TET-", names(positive_hits)))) {
      return(data.frame("new_hash_id" = "other_doublet",
                        "full_hash_id" = paste(names(positive_hits),
                                               collapse = "_")))
    } else {
      return(data.frame("new_hash_id" = "islet_reactive_doublet",
                        "full_hash_id" = paste(names(positive_hits),
                                               collapse = "_")))
    }
  })
  
  new_hash_id <- do.call(rbind, new_hash_id)
  
  return(new_hash_id)
}

positive_quantile <- 0.90

# Run with the old demux and save the output
test_object <- seurat_object

test_object <- HTODemux(test_object, assay = "TET", 
                        positive.quantile = positive_quantile,
                        kfunc = "kmeans")

test_object$tet_hash_id <- test_object$hash.ID

test_object <- HTODemux(test_object, assay = "SCAR_TET", 
                        positive.quantile = positive_quantile,
                        kfunc = "kmeans")

test_object$scar_hash_id <- test_object$hash.ID

test_meta <- test_object[[]] %>%
  dplyr::select(dplyr::contains("TET_"), dplyr::contains("SCAR_TET_"),
                dplyr::contains("hash_id"))

write.csv(test_meta, file.path(save_dir, "files", "previous_htodemux_res.csv"))


rm(test_object)
rm(test_meta)

# Run with the new demux and save the object
return_res <- HTODemuxUpdate(seurat_object, assay = "TET", 
                             positive.quantile = positive_quantile,
                             kfunc = "kmeans")

seurat_object <- return_res$object
proportions <- return_res$proportions

# Based on these proporitions, add a new column that is
# "islet doublet" if there is more than 2 > 1 and all > 1 are islet reactive
# "other dobulet" if there is more than 2 > 1 and not all > 1 are islet reactive
# singlet based on name if there is 1 > 1
# negative based on name if there is 0 > 1
new_hash_id <- find_hash_id(proportions)

seurat_object <- AddMetaData(seurat_object, metadata = new_hash_id)


# Make an assay for the proportions
seurat_object[["TET_PROPORTIONS"]] <- CreateAssayObject(data = proportions)

seurat_object$tet_hash_id <- seurat_object$hash.ID

return_res <- HTODemuxUpdate(seurat_object, assay = "SCAR_TET", 
                                positive.quantile = positive_quantile,
                                kfunc = "kmeans")

seurat_object <- return_res$object
proportions <- return_res$proportions

# Based on these proporitions, add a new column that is
# "islet doublet" if there is more than 2 > 1 and all > 1 are islet reactive
# "other dobulet" if there is more than 2 > 1 and not all > 1 are islet reactive
# singlet based on name if there is 1 > 1
# negative based on name if there is 0 > 1
new_hash_id <- find_hash_id(proportions)

colnames(new_hash_id) = paste0("scar_", colnames(new_hash_id))

seurat_object <- AddMetaData(seurat_object, metadata = new_hash_id)

# Make an assay for the proportions
seurat_object[["SCAR_TET_PROPORTIONS"]] <- CreateAssayObject(data = proportions)


seurat_object$scar_hash_id <- seurat_object$hash.ID



# seurat_object <- HTODemux(seurat_object, assay = "SCAR_TET",
#                           positive.quantile = 0.90,
#                           kfunc = "kmeans")
# 
# table(seurat_object$scar_hash_id, seurat_object$hash.ID)

cm <- confusionMatrix(seurat_object$tet_hash_id,
                      seurat_object$scar_hash_id)

cm <- cm / rowSums(cm)

graphics.off()

pdf(file.path(save_dir, "images", "hto_demux_heatmap.pdf"))

print(pheatmap::pheatmap(cm, main = "raw_vs_scar"))

dev.off()

 if(HTO){
  # Demultiplex
  seurat_object <- HTODemux(seurat_object, assay = "HTO",
                            positive.quantile = 0.90)


  # Plots
  ridge_p <- RidgePlot(seurat_object, assay = "HTO",
                       features = rownames(stoc_data[["HTO"]])[1:2], ncol = 2)
  scatter_p <- FeatureScatter(seurat_object,
                              feature1 = "hto_Hashtag-STOC-86-008",
                              feature2 = "hto_Hashtag-STOC-86-009")
  Idents(seurat_object) <- "HTO_classification.global"
  vln_p <- VlnPlot(seurat_object, features = "nCount_RNA",
                   pt.size = 0.1, log = TRUE)

  pdf(file.path(save_dir, "images", "HTO_plots.pdf"))
  ridge_p
  scatter_p
  vln_p
  dev.off()
}

saveRDS(seurat_object, file = file.path(save_dir, "rda_obj",
                                        "seurat_start.rds"))
