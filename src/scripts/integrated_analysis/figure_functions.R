#' Generates an hclust object based on the comparison (query) and a status.
#' 
#' This assumes that gene ontology was run using `gprofiler2::gost`. The
#' final output will be a dotplot where the size is the number of
#' intersectig genes and the color is the log of the p value. These
#' dots will be on the end of the dendrogram built on the similarity 
#' between kett terms. Currently this function only works for KEGG because
#' the genes in each term need to be downloaded and only KEGG is currently
#' supported
#' 
#' @param query The query from the gse_res, the gse_res will be subset
#' by this query
#' @param gse_res the results table returned by `gprofiler2::gost`
#' @param p_val_cutoff what p_val to use to subset terms
make_histogram <- function(query, gse_res, p_val_cutoff = 0.01){
  
  # Subset to only kegg terms and the query of interest
  gse_res <- gse_res[gse_res$source == "KEGG" & 
                       gse_res$query == query,]
  
  # Make a list of all terms
  gse_list <- gse_res$term_id
  names(gse_list) <- gse_res$term_name
  
  # Find all genes associated with the kegg term
  # This could probably go faster if all genes from all kegg terms are
  # first pulled and keggConv and mapIds are only run once.
  all_genes <- lapply(gse_list, function(x){
    genes <- get_gene_list(x)
    
    return(genes)
  })
  
  # Make a matrix of the correct size, it should be all gene sets by all gene
  # sets
  total_values <- length(all_genes) ^ 2
  final_jaccard_matrix <- Matrix::Matrix(data = rep(0, total_values),
                                         nrow = length(all_genes),
                                         ncol = length(all_genes),
                                         sparse = TRUE)
  
  # Find the jaccard distance for all pairwise gene sets
  for (x in 1:length(all_genes)){
    for (y in 1:length(all_genes)){
      
      save_val <- jaccard(a = all_genes[[x]],
                          b = all_genes[[y]])
      
      # Save to both relevant positions in the matrix so you don't
      # run twice
      final_jaccard_matrix[x, y] <- save_val
      final_jaccard_matrix[y, x] <- save_val
    }
  }
  
  colnames(final_jaccard_matrix) <- names(all_genes)
  rownames(final_jaccard_matrix) <- names(all_genes)
  
  # Make into a distance matrix and run hclust
  hc_full <- stats::hclust(stats::as.dist(1- final_jaccard_matrix),
                           method = "ward.D2")
  
  # Subset based on the desired p-value
  gse_small <- gse_res[gse_res$p_value < p_val_cutoff,]
  
  jaccard_small <- final_jaccard_matrix[rownames(final_jaccard_matrix)
                                        %in% gse_small$term_name,
                                        colnames(final_jaccard_matrix)
                                        %in% gse_small$term_name]
  
  # Repeat clustering on the smaller data set
  hc_small <- stats::hclust(stats::as.dist(1- jaccard_small),
                            method = "ward.D2")
  
  return(list(full = hc_full, small = hc_small, gse_res = gse_res))
  
}

#' Finds all gene IDs from a a kegg term
#' 
#' This takes a kegg term and a path of all kegg terms and maps
#' to ncbi gene ids
#' 
#' @param kegg_term The list of kegg terms to query. If multiple are provided
#' a list of genes from all terms combined will be returned.
get_gene_list <- function(kegg_term){
  path <- keggLink("pathway", "hsa")
  
  search_term <- gsub("KEGG:", "", kegg_term)
  
  tryCatch({
    # Find genes in term
    all_genes <- lapply(search_term, function(term){
      genes_in_term <- path[grepl(term, path)]
      return(genes_in_term)
    })
    all_genes <- unlist(all_genes)
    
    # Translate to gene id
    all_gene_id <- keggConv("ncbi-geneid", unique(names(all_genes)))
    
    all_gene_id <- gsub("ncbi-geneid:", "", all_gene_id)
    
    # Map to gene symbols
    gene_symbols <- mapIds(org.Hs.eg.db, keys = all_gene_id, 
                           keytype = "ENTREZID", column = "SYMBOL")  
    
  }, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
    gene_symbols <- "NA"
  })
  
  return(as.character(gene_symbols))
}

#' Makes a histogram based on the output from `make_histogram`.
#' 
#' This assumes that `make_histogram` has already been run. Will take
#' that output and make a dendrogram
#' 
#' @param hc The hclust object output from `make_h`
#' @param gse_res the results table returned by `gprofiler2::gost`
#' @param y_lim how to set the y-lim, this is important to modify if
#' you save the figure in different sizes, it moves where the key is 
#' relative to the text
#' @param title The title of the plot
plot_histogram <- function(hc, gse_res, ylim = c(2, -2), title = NULL){
  # Convert hclust to dendrogram
  dend <- as.dendrogram(hc)
  
  # Make dendrogram into an object that plays well with ggplot
  ddata <- ggdendro::dendro_data(dend, type = "rectangle")
  
  # Add gse information
  dot_data <- merge(ddata$labels, gse_res, by.x = "label", by.y = "term_name",
                    all.x = TRUE, all.y = FALSE)
  
  
  # Make the plot with the dendrogram passed to geom segment and the 
  # gse res as the geom_point
  p <- ggplot(ggdendro::segment(ddata)) + 
    ggplot2::geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    ggplot2::geom_text(data = ddata$labels, ggplot2::aes(x = x, y = y - 0.06,
                                                         label = label),
                       hjust = 0, angle = 0, size = 3) +
    ggplot2::geom_point(data = dot_data, ggplot2::aes(x = x, y = y,
                                                      color = -log(p_value),
                                                      size = intersection_size)) +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    ggplot2::coord_flip() + 
    ggplot2::scale_y_reverse(expand = c(0.2, 0)) +
    ggplot2::ylim(ylim) +
    ggplot2::theme(axis.line = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title = element_blank())
  
  if(!is.null(title)){
    p <- p +
      ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  
  return(p)
  
}

#' Makes a volcano
#' 
#' Makes a volcano plot based on a previously run DE analysis. Assumes 
#' columns named `p_val_adj`, `avg_log2FC` and `gene` are present. Also
#' takes a list of kegg pathways. Any gene in any of those pathways will be
#' colored on top of the volcano plot. Logic to fix the p-values was borrowed
#' from `EnhancedVolcano`
#' 
#' @param de_genes The de genes to plot
#' @param kegg_pathways A list of kegg pathways. Genes from these will
#' be highlighted
#' @param group_name Name for the kegg pathways. This will be the name
#' of the plot
#' @param color What color to color the kegg pathway genes
make_volcano <- function(de_genes, kegg_pathways, group_name,
                         color = NULL, plotly = FALSE){
  # Pull out all genes to highlight on the volcano plot
  highlight_genes <- get_gene_list(kegg_pathways)
  
  if(is.null(color)){
    color <- "#C41E3A"
  }
  
  # Fix y values (taken from enhanced volcano)
  if (min(de_genes$p_val_adj, na.rm=TRUE) == 0) {
    # <= v1.2
    #warning(paste("One or more P values is 0.",
    #  "Converting to minimum possible value..."),
    #  call. = FALSE)
    #toptable[which(toptable[[y]] == 0), y] <- .Machine$double.xmin
    warning(paste('One or more p-values is 0.',
                  'Converting to 10^-1 * current',
                  'lowest non-zero p-value...'),
            call. = FALSE)
    
    de_genes[which(de_genes[["p_val_adj"]] == 0), "p_val_adj"] <- min(
      de_genes[which(de_genes[["p_val_adj"]] != 0), "p_val_adj"],
      na.rm = TRUE) * 10^-1
  }
  
  # Add in colors
  de_genes$color <- ifelse(de_genes$gene %in% highlight_genes  & 
                             de_genes$p_val_adj < 0.05 &
                             abs(de_genes$avg_log2FC) > 0.25,
                           "kegg_gene", "other_gene")
  
  de_genes$color <- factor(de_genes$color, levels = c("other_gene",
                                                      "kegg_gene"))
  
  de_genes <- de_genes[order(de_genes$color),]
  
  color_mapping <- c("kegg_gene" = color, "other_gene" = "#999999")
  
  
  # Make a volcano plot
  if(plotly){
    volcano_plot <- ggplot2::ggplot(de_genes, ggplot2::aes(x = avg_log2FC,
                                                           y = -log(p_val_adj),
                                                           color = color,
                                                           text = gene))  
  } else {
    volcano_plot <- ggplot2::ggplot(de_genes, ggplot2::aes(x = avg_log2FC,
                                                           y = -log(p_val_adj),
                                                           color = color))
  }
  volcano_plot <- volcano_plot +
    ggplot2::geom_point() + 
    ggplot2::xlim(-2, 2) +
    ggplot2::scale_color_manual(values = color_mapping) +
    ggplot2::ggtitle(group_name) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  return(volcano_plot)
}

make_circos_plot <- function(circos_df, color = NULL,
                             grid_color = NULL){
  if(nrow(circos_df) > 0){
    
    chordDiagram(circos_df, annotationTrack = "grid", 
                 preAllocateTracks = 1, col = color,
                 grid.col = grid_color)
    
    # we go back to the first track and customize sector labels
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                  facing = "clockwise", niceFacing = TRUE, adj = c(-0.1, 0.5))
    }, bg.border = NA) # here set bg.border to NA is important
    
  }
  
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
    select_cols <- c(group_by)
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
      dplyr::select(dplyr::all_of(c(select_cols, "IGH", chain))) %>%
      dplyr::filter(!is.na(!!as.name(chain)))
    test_chain <- chain
  }
  
  return_df <- return_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(select_cols, test_chain, "IGH")))) %>% 
    dplyr::add_count(name = "value") %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    dplyr::filter(value >= subset_counts)
  
  colnames(return_df) <- c(select_cols, "heavy", "light", "value")
  
  # Order the data frame
  heavy_order <- strcapture("([A-Za-z]+)([0-9]+)-([0-9]+)", 
                            unique(return_df$heavy),
                            list(prefix = "", numeric1 = 0, numeric2 = 0))
  
  heavy_order$original <- unique(return_df$heavy)
  
  heavy_order <- data.frame(heavy_order)
  
  heavy_order <- heavy_order[order(heavy_order$numeric1, heavy_order$numeric2),]
  
  light_order <- strcapture("([A-Za-z]+)([0-9]+)-([0-9]+)", 
                            unique(return_df$light),
                            list(prefix = "", numeric1 = 0, numeric2 = 0))
  
  light_order$original <- unique(return_df$light)
  
  light_order <- data.frame(light_order)
  
  light_order <- light_order[order(light_order$prefix,
                                   light_order$numeric1, 
                                   light_order$numeric2),]  
  
  return_df$heavy <- factor(return_df$heavy, levels = heavy_order$original)
  return_df$light <- factor(return_df$light, levels = light_order$original)
  
  return_df <- return_df[order(return_df$heavy,
                               return_df$light), ]
  
  if(group_by != "all"){
    
    return_df$color <- color_list[as.character(return_df[[group_by]])]
    color_values <- return_df$color
    full_df <- return_df
    return_df <- return_df %>%
      dplyr::select(-dplyr::all_of(group_by), -color)
  } else{
    color_values <- NULL
    full_df <- return_df
  }
  
  return(list("df" = return_df, "color_list" = color_values,
              "full_df" = full_df))
}
