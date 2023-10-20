make_plots <- function(clone_info_expanded, name, percent_keep = 25,
                       locus_use = "IGH", clone_min = 3){
  expanded <- clone_info_expanded %>%
    dplyr::filter(locus == locus_use, clone_count > clone_min) 
  
  # Number of clones across samples
  tetramer_group_across <- expanded %>%
    dplyr::group_by(clone_id, tet_hash_id) %>%
    dplyr::add_count(name = "tetramer_count") %>%
    dplyr::mutate(tetramer_fraction = tetramer_count / clone_count * 100) %>%
    dplyr::mutate(clone_name = paste0("clone_", clone_id)) %>%
    dplyr::select(clone_id, clone_name, tet_hash_id, tetramer_count,
                  tetramer_fraction, clone_count) %>%
    dplyr::distinct()
  
  plot_one <- single_plot(tetramer_group = tetramer_group_across,
                          name = name, type = "across",
                          percent_keep = percent_keep)
  
  # Number of clones within samples
  expanded <- clone_info_expanded %>%
    dplyr::filter(locus == locus_use, sample_clone_count > clone_min) 
  
  tetramer_group_within <- expanded %>%
    dplyr::group_by(sample_clone_id, tet_hash_id) %>%
    dplyr::add_count(name = "tetramer_count") %>%
    dplyr::mutate(tetramer_fraction = tetramer_count / sample_clone_count * 100) %>%
    dplyr::mutate(sample_clone_name = paste0("clone_", sample_clone_id)) %>%
    dplyr::select(sample_clone_id, sample_clone_name, tet_hash_id, tetramer_count,
                  tetramer_fraction, sample_clone_count, sample) %>%
    dplyr::distinct()
  
  colnames(tetramer_group_within) <- gsub("sample_", "", 
                                          colnames(tetramer_group_within))
  
  plot_two <- single_plot(tetramer_group = tetramer_group_within,
                          name = name, type = "within",
                          percent_keep = percent_keep)
  
  # Clones shared between samples
  expanded <- clone_info_expanded %>%
    dplyr::filter(locus == locus_use, clone_count > clone_min) 
  
  tetramer_group_sample <- expanded %>%
    dplyr::group_by(clone_id, tet_hash_id) %>%
    dplyr::add_count(name = "tetramer_count") %>%
    dplyr::mutate(tetramer_fraction = tetramer_count / clone_count * 100) %>%
    dplyr::mutate(clone_name = paste0("clone_", clone_id)) %>%
    dplyr::select(clone_id, clone_name, tet_hash_id, tetramer_count,
                  tetramer_fraction, clone_count, sample)
  
  keep_clones <- tetramer_group_sample %>%
    dplyr::ungroup() %>%
    dplyr::select(clone_name, sample) %>%
    dplyr::distinct() %>%
    dplyr::group_by(clone_name) %>%
    dplyr::add_count(name = "number_samples") %>%
    dplyr::filter(number_samples > 1)
  
  tetramer_group_sample <- tetramer_group_sample %>%
    dplyr::filter(clone_name %in% unique(keep_clones$clone_name)) %>%
    dplyr::select(clone_id, clone_name, tet_hash_id, tetramer_count,
                  tetramer_fraction, clone_count) %>%
    dplyr::distinct()
  
  plot_three <- single_plot(tetramer_group = tetramer_group_sample,
                            name = name, type = "between",
                            percent_keep = percent_keep)
  
  return(list("across" = plot_one,
              "within" = plot_two,
              "between" = plot_three))
  
  
}

single_plot <- function(tetramer_group, name,
                        type, percent_keep = 25){
  
  if(type == "across"){
    save_name <- paste0("All clones across ", name, " samples")
  } else if(type == "within"){
    save_name <- paste0("All clones within ", name, " samples")
  } else if(type == "between"){
    save_name <- paste0("All clones shared between ", name, " samples")
  } else {
    stop("Type can only be across, within, or between")
  }
  
  interesting_clones <- tetramer_group %>%
    dplyr::filter(tet_hash_id %in% c("GAD-tet", "IA2-tet",
                                     "INS-tet")) %>%
    dplyr::group_by(clone_name) %>%
    dplyr::mutate(percent_diabetes = sum(tetramer_fraction)) %>%
    dplyr::filter(percent_diabetes > percent_keep)
  
  plot_data <- tetramer_group %>%
    dplyr::filter(clone_name %in% unique(interesting_clones$clone_name))
  
  count_data <- plot_data %>%
    dplyr::select(clone_name, clone_count) %>%
    dplyr::mutate(clone_count = paste0(clone_count, " cells")) %>%
    dplyr::distinct()
  
  return_plot1 <-ggplot2::ggplot(plot_data, ggplot2::aes(x = clone_name,
                                                        y = tetramer_fraction,
                                                        fill = tet_hash_id)) +
    ggplot2::geom_bar(position = "stack", stat = "identity") +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggplot2::scale_fill_manual(values = tetramer_colors) +
    ggplot2::ggtitle(save_name)
  
  return_plot2 <- return_plot1 + 
    ggplot2::geom_text(
      data = count_data, 
      ggplot2::aes(clone_name, 50, label = clone_count),
      show.legend = F,
      inherit.aes = F,
      color = "white", 
      angle = 90
    )
  
  if(nrow(plot_data) > 0){
    return(list(return_plot1, return_plot2))
  } else {
    return(NULL)
  }
}
