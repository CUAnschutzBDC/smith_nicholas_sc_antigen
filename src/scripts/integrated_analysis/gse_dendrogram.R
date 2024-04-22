#library(enrichplot)
library(here)
library(KEGGREST)
library(org.Hs.eg.db)
library(ggraph)
library(ggdendro)

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))


save_dir <- here("results/R_analysis/merged")
gse_res <- read.csv(file.path(save_dir, "files/de/all_GSE_results.csv"))

name_mapping <- c("T1D" = "T1D_Islet_Reactive_ND_Islet_Reactive",
                  "AAB" = "AAB_Islet_Reactive_ND_Islet_Reactive")

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


make_histogram <- function(query, status, gse_res, p_val_cutoff = 0.01){

  gse_res <- gse_res[gse_res$source == "KEGG" & 
                           gse_res$query == query,]

  gse_list <- gse_res$term_id
  names(gse_list) <- gse_res$term_name
  
  path <- keggLink("pathway", "hsa")
  
  
  all_genes <- lapply(gse_list, function(x){
    search_term <- gsub("KEGG:", "", x)
    
    tryCatch({
      genes_in_term <- path[grepl(search_term, path)]
      
      all_gene_id <- keggConv("ncbi-geneid", names(genes_in_term))
      
      all_gene_id <- gsub("ncbi-geneid:", "", all_gene_id)
      
      gene_symbols <- mapIds(org.Hs.eg.db, keys = all_gene_id, 
                             keytype = "ENTREZID", column = "SYMBOL")  
      
    }, error = function(e) {
      cat("Error:", conditionMessage(e), "\n")
      gene_symbols <- "NA"
    })
    
    return(as.character(gene_symbols))
  })
  
  total_values <- length(all_genes) ^ 2
  final_jaccard_matrix <- Matrix::Matrix(data = rep(0, total_values),
                                         nrow = length(all_genes),
                                         ncol = length(all_genes),
                                         sparse = TRUE)
  
  for (x in 1:length(all_genes)){
    for (y in 1:length(all_genes)){
      
      save_val <- jaccard(a = all_genes[[x]],
                          b = all_genes[[y]])
      final_jaccard_matrix[x, y] <- save_val
      final_jaccard_matrix[y, x] <- save_val
    }
  }
  
  colnames(final_jaccard_matrix) <- names(all_genes)
  rownames(final_jaccard_matrix) <- names(all_genes)
  
  hc_full <- stats::hclust(stats::as.dist(1- final_jaccard_matrix),
                      method = "ward.D2")
  
  gse_small <- gse_res[gse_res$p_value < p_val_cutoff,]
  
  jaccard_small <- final_jaccard_matrix[rownames(final_jaccard_matrix)
                                        %in% gse_small$term_name,
                                        colnames(final_jaccard_matrix)
                                        %in% gse_small$term_name]
  
  hc_small <- stats::hclust(stats::as.dist(1- jaccard_small),
                            method = "ward.D2")
  
  return(list(full = hc_full, small = hc_small, gse_res = gse_res))
  
}

plot_histogram <- function(hc, gse_res, ylim = c(2, -2), title = NULL){
  # Convert hclust to dendrogram
  dend <- as.dendrogram(hc)
  
  
  ddata <- dendro_data(dend, type = "rectangle")
  
  dot_data <- merge(ddata$labels, gse_res, by.x = "label", by.y = "term_name",
                    all.x = TRUE, all.y = FALSE)
  
  
  p <- ggplot(segment(ddata)) + 
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

hc_res_t1d <- make_histogram(status = names(name_mapping)[[1]],
                             query = name_mapping[[1]], 
                             gse_res = gse_res)


hc_res_aab <- make_histogram(status = names(name_mapping)[[2]],
                             query = name_mapping[[2]], 
                             gse_res = gse_res)

pdf(file.path(save_dir, "images", "final_figures",
              "kegg_clustering_t1d_all.pdf"), 
    height = 12, width = 9) 

plot_histogram(hc = hc_res_t1d$full, ylim = c(2, -1),
               title = "T1D vs ND pathways",
               gse_res = hc_res_t1d$gse_res)

dev.off()


pdf(file.path(save_dir, "images", "final_figures",
              "kegg_clustering_t1d_pval_01.pdf"), 
    height = 9, width = 9) 

plot_histogram(hc = hc_res_t1d$small, ylim = c(2, -1),
               title = "T1D vs ND pathways",
               gse_res = hc_res_t1d$gse_res)

dev.off()


pdf(file.path(save_dir, "images", "final_figures",
              "kegg_clustering_aab_all.pdf"), 
    height = 12, width = 9) 

plot_histogram(hc = hc_res_aab$full, ylim = c(2.1, -1),
               title = "AAB vs ND pathways",
               gse_res = hc_res_aab$gse_res)

dev.off()


pdf(file.path(save_dir, "images", "final_figures",
              "kegg_clustering_aab_pval_01.pdf"), 
    height =6, width = 9) 

plot_histogram(hc = hc_res_aab$small, ylim = c(2, -1),
               title = "AAB vs ND pathways",
               gse_res = hc_res_aab$gse_res)

dev.off()
