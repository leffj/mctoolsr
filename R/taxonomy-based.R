# mctoolsr

######################
# TAXONOMY FUNCTIONS #
######################

#' @keywords internal
# generates a data frame with all levels of taxanomic info as columns
.compile_taxonomy = function(biom_dat){
  # get only taxonomy observation metadata from a biom file
  obs_md = biom::observation_metadata(biom_dat)
  if(class(obs_md) == 'list'){
    # replace label for otus with only 1 taxonomy level
    obs_md = sapply(obs_md, function(x) {
      names(x) = gsub('taxonomy$', 'taxonomy1', names(x))
      list(x)})
    # get the taxonomy levels
    otu_full_md = names(which.max(sapply(obs_md, length))) # need to get the otu with all metadata levels
    taxa_levels = names(obs_md[otu_full_md][[1]])[
      grepl('taxonomy', names(obs_md[[1]]), ignore.case = TRUE)]
    # compile taxonomy for each level
    tax_comp = data.frame(row.names = names(obs_md))
    for(l in 1:length(taxa_levels)){
      level_tax_tmp = sapply(obs_md, function(x) x[taxa_levels[l]])
      names(level_tax_tmp) = names(obs_md)
      tax_comp[, taxa_levels[l]] = level_tax_tmp
    }
    tax_comp[is.na(tax_comp)] = 'unclassified'
    tax_comp
  } else if(class(obs_md) == 'data.frame'){
    if(ncol(obs_md) < 2) { 
      #can happen if taxonomy loaded as another type of metadata in biom
      warning(paste0('Error compiling taxonomy. Check that taxonomy is ', 
                    'formatted correctly in biom file. Proceeding without ', 
                    'taxonomy.'))
      NULL
    } else obs_md
  } else {
    warning(paste0('Error compiling taxonomy. Check that taxonomy is ', 
                   'formatted correctly in biom file. Proceeding without ', 
                   'taxonomy.'))
    NULL
  }
}

#' @keywords internal
# Convert a taxonomy character vector pulled from a text OTU 
# table to a data frame
.parse_taxonomy = function(taxonomy_vec){
  tmp = sapply(taxonomy_vec, function(x) strsplit(as.character(x), 
                                                  split = '; *')[[1]])
  # check for long taxonomy rows, which may indicate a problem
  lengths = sapply(tmp, length)
  if(max(lengths) > 15) {
    prob_line = which(lengths == max(lengths))
    warning(paste0('Issue with provided taxonomy. Check for quotes near line ', 
                   prob_line, '. Proceeding without taxonomy.'))
    return(NULL)
  }
  # if not all taxonomic levels present for each OTU, this is necessary to 
  # convert to data frame
  if(class(tmp) == 'list'){
    n_obs = sapply(tmp, length)
    seq_max = seq_len(max(n_obs))
    tmp = (sapply(tmp, "[", i = seq_max))
  }
  tmp = as.data.frame(t(tmp))
  potential_names = c('taxonomy1', 'taxonomy2', 'taxonomy3', 'taxonomy4', 
                      'taxonomy5', 'taxonomy6', 'taxonomy7', 'taxonomy8', 
                      'taxonomy9', 'taxonomy10', 'taxonomy11', 'taxonomy12')
  names(tmp) = potential_names[1:ncol(tmp)]
  tmp
}

#' @title Calculate Values for Coarser Taxonomic Groupings
#' @description Given input as generated from \code{load_taxon_table()}, 
#'   calculate relative abundances (or absolute abundances) of taxa at a coarser
#'   taxonomic threshold. Taxonomy must have been loaded for this function to 
#'   work.
#' @param input The input dataset
#' @param level The taxonomy level to summarize by [integer]
#' @param relative Return relative abundances or not (default = TRUE)
#' @param report_higher_tax Whether or not to return taxonomic strings higher
#'   than the indecated level (default = TRUE)
#' @examples 
#' ## Return relative abundances of phyla for each sample
#' \dontrun{
#' summarize_taxonomy(input = example_input, level = 2)
#' }
summarize_taxonomy = function(input, level, relative = TRUE, 
                              report_higher_tax = TRUE){
  if(report_higher_tax) taxa_strings = apply(input$taxonomy_loaded[1:level], 1, 
                                             paste0, collapse = '; ')
  else taxa_strings = input$taxonomy_loaded[, level]
  no_taxa = length(unique(taxa_strings))
  tax_sum = apply(input$data_loaded, 2, function(x) by(x, taxa_strings, sum))
  if(no_taxa == 1) tax_sum = data.frame(t(tax_sum), 
                                        row.names = unique(taxa_strings))
  else tax_sum = as.data.frame(tax_sum)
  if(relative){
    convert_to_relative_abundances(tax_sum)
  } else tax_sum
}

#' @title Plot stacked bar plots to represent taxa compompositions
#' @description Stacked bar plots will be generated for each factor level 
#'   indicated by \code{group_factor} to display their taxonomic compositions. 
#'   Only the top few taxa will be displayed as indicated by \code{num_taxa}. 
#'   Mean values are calculated within the factor levels.
#' @param taxa_summary_df The taxon summary dataframe (the output of \code{
#'   summarize_taxonomy()})
#' @param metadata_map The mapping file dataframe
#' @param group_factor The factor (metadata header label) used to create the 
#'   bars. Means will be taken for each factor level.
#' @param num_taxa The number of top most abundant taxa to display. Additional 
#'   will be grouped into "Other".  
plot_taxa_bars = function(taxa_summary_df, metadata_map, group_factor, num_taxa){
  taxa_summary_df$taxon = row.names(taxa_summary_df)
  taxa_summary_df_melted = reshape2::melt(taxa_summary_df, variable.name = 'Sample_ID', 
                                          id.vars = 'taxon')
  group_by_levels = metadata_map[match(taxa_summary_df_melted$Sample_ID, 
                                       row.names(metadata_map)), group_factor]
  taxa_summary_df_melted$group_by = group_by_levels
  mean_tax_vals = dplyr::summarise(dplyr::group_by(taxa_summary_df_melted, 
                                                   group_by, taxon), 
                                   mean_value = mean(value))
  # get top taxa and convert other to 'other'
  mean_tax_vals_sorted = mean_tax_vals[order(mean_tax_vals$mean_value, 
                                             decreasing = TRUE), ]
  top_taxa = unique(mean_tax_vals_sorted$taxon)[1:num_taxa]
  mean_tax_vals_sorted$taxon[!mean_tax_vals_sorted$taxon %in% top_taxa] = 'Other'
  to_plot = dplyr::summarise(dplyr::group_by(mean_tax_vals_sorted, group_by, 
                                             taxon), 
                             mean_value = sum(mean_value))
  # plot
  ggplot2::ggplot(to_plot, ggplot2::aes(group_by, mean_value, fill = taxon)) +
    ggplot2::geom_bar(stat = 'identity') + 
    ggplot2::ylab('') + ggplot2::xlab('') +
    ggplot2::theme(legend.title = ggplot2::element_blank())
}

#' @title Filter taxa from an individual taxa summary table
#' @details Can use one or more of the parameters to do filtering. Threshold 
#'  filtering takes precidence over taxa filtering. If taxa to keep and 
#'  taxa to remove are both included, taxa to remove will be 
#'  removed if the parameter entries conflict.
#' @param tax_table Input taxa summary table from \code{summarize_taxonomy()} 
#'  (dataframe).
#' @param filter_thresh Filter taxa less than this number based on mean 
#'  table values.
#' @param taxa_to_keep Keep only taxa that contain these names. Vector or string.
#' @param taxa_to_remove Remove taxa that contain these names. Vector or string.
filter_taxa_from_table = function(tax_table, filter_thresh = 0, taxa_to_keep, 
                                  taxa_to_remove){
  means = apply(tax_table[, 1:ncol(tax_table)], 1, 
                function(x){mean(x,na.rm=TRUE)})
  number_retained = sum((means >= filter_thresh) *1)
  taxa_keep = names(means[means >= filter_thresh])
  if(!missing(taxa_to_keep)) {taxa_keep = taxa_keep[taxa_keep %in% taxa_to_keep]}
  if(!missing(taxa_to_remove)) {taxa_keep = taxa_keep[!taxa_keep %in% taxa_to_remove]}
  tax_table[row.names(tax_table) %in% taxa_keep, ]
}

#' @title Filter taxa from a loaded dataset
#' @details Can use one or more of the parameters to do filtering. Threshold 
#'          filtering takes precidence over taxa filtering. If taxa to keep and 
#'          taxa to remove are both included, taxa to remove will be 
#'          removed if the parameter entries conflict. Taxa are found using 
#'          grep, so it is not necessary to use the entire names of taxa as long
#'          as they are not ambiguous. Can filter using taxa IDs (i.e. row 
#'          names in taxa table).
#' @param input Input data (a list variable) from \code{load_taxa_table()}.
#' @param filter_thresh Filter OTUs less than this number based on mean OTU 
#'        table values.
#' @param taxa_to_keep Keep only taxa that contain these names. Vector or string.
#' @param taxa_to_remove Remove taxa that contain these names. Vector or string.
#' @param at_spec_level If included, only keep/remove matches at this specific 
#'        taxonomy level(s) (a number/numbers referring to the taxonomy 
#'        column(s)).
#' @param taxa_IDs_to_keep Keep only taxa with these IDs (row names).
#' @param taxa_IDs_to_remove Remove taxa with these IDs (row names).
filter_taxa_from_input = function(input, filter_thresh, taxa_to_keep, 
                                 taxa_to_remove, at_spec_level, 
                                 taxa_IDs_to_keep, taxa_IDs_to_remove){
  rows_keep = seq(1, nrow(input$data_loaded))
  if(!missing(filter_thresh)){
    means = apply(input$data_loaded[, 1:ncol(input$data_loaded)], 1, 
                  function(x){mean(x, na.rm=TRUE)})
    number_retained = sum((means >= filter_thresh) *1)
    rows_keep = rows_keep[means >= filter_thresh]
  }
  # specify taxa levels to search
  if(missing(at_spec_level)){
    tax_levels = 1:ncol(input$taxonomy_loaded)
  } else tax_levels = at_spec_level
  # if particular taxa specified to keep, identify those rows.
  if(!missing(taxa_to_keep)){
    rows_keep_tmp = sapply(taxa_to_keep, FUN = function(x){
      grep(x, apply(as.data.frame(input$taxonomy_loaded[, tax_levels]), 1, 
                    paste0, collapse = ''))
    })
    if(length(rows_keep_tmp[[1]]) == 0){
      stop('Taxa not found.')
    }
    rows_keep = intersect(rows_keep, rows_keep_tmp)
  }
  # if particular taxa IDs specified to keep, identify those rows.
  if(!missing(taxa_IDs_to_keep)){
    rows_keep_tmp = sapply(taxa_IDs_to_keep, FUN = function(x){
      match(x, row.names(input$taxonomy_loaded))
    })
    if(length(rows_keep_tmp[[1]]) == 0){
      stop('Taxa IDs not found.')
    }
    rows_keep = intersect(rows_keep, rows_keep_tmp)
  }
  # if particular taxa to remove, identify those rows
  if(!missing(taxa_to_remove)){
    rows_remove = sapply(taxa_to_remove, FUN = function(x){
      grep(x, apply(as.data.frame(input$taxonomy_loaded[, tax_levels]), 1, 
                    paste0, collapse = ''))
    })
    if(length(rows_remove[[1]]) == 0){
      stop('Taxa not found.')
    }
    rows_keep = rows_keep[! rows_keep %in% unlist(rows_remove)]
  }
  # if particular taxa IDs to remove, identify those rows
  if(!missing(taxa_IDs_to_remove)){
    rows_remove = sapply(taxa_IDs_to_remove, FUN = function(x){
      match(x, row.names(input$taxonomy_loaded))
    })
    if(length(rows_remove[[1]]) == 0){
      stop('Taxa IDs not found.')
    }
    rows_keep = rows_keep[! rows_keep %in% unlist(rows_remove)]
  }
  output = list(data_loaded = input$data_loaded[rows_keep, ],
                map_loaded = input$map_loaded, 
                taxonomy_loaded = droplevels(input$taxonomy_loaded[rows_keep, 
                                                                   ]))
  removed = nrow(input$data_loaded) - nrow(output$data_loaded)
  message(paste0(removed, ' taxa removed'))
  output
}


#' @title Plot Taxa Summary Heatmap
#' @description A quick way to create a heatmap from a taxa summary dataframe.
#'  Samples are grouped by a category that is specified in 'type_header'.
#' @param tax_table A taxa table dataframe. The output of 
#'  \code{summarize_taxonomy()}.
#' @param metadata_map A metadata mapping dataframe. Usually, 
#'  \code{input$map_loaded}
#' @param min_rel_abund The minimum mean relative abundance for a taxon to not
#'  be grouped into 'Other'. Between 0 and 1.
#' @param type_header The metadata_map header label used to group samples.
#' @param scale_by Whether to scale colors by (a) 'sample_types', (b) 'taxa', or
#'  (c) 'all'.
plot_ts_heatmap = function(tax_table, metadata_map, min_rel_abund, type_header, 
                           scale_by) {
  # group all taxa lower than threshold into other
  lt_thresh = tax_table[rowMeans(tax_table) < min_rel_abund, ]
  gt_thresh = tax_table[rowMeans(tax_table) >= min_rel_abund, ]
  Other = colSums(lt_thresh)
  sumtax_mod = rbind(gt_thresh, Other = Other)
  # get means
  sumtax_smry = taxa_summary_by_sample_type(sumtax_mod, metadata_map, 
                                            type_header, smry_fun = mean)
  sumtax_smry = round(sumtax_smry*100, 1)
  melted = reshape2::melt(sumtax_smry)
  if (scale_by == 'sample_types') {
    to_plot = dplyr::mutate(dplyr::group_by(melted, Var2), 
                            scaled = scales::rescale(value))
  } else if (scale_by == 'taxa') {
    to_plot = dplyr::mutate(dplyr::group_by(melted, Var1), 
                            scaled = scales::rescale(value))
  } else if (scale_by == 'all') {
    to_plot = dplyr::mutate(melted, scaled = scales::rescale(value))
  } else stop("scale_by must be 'sample_types' or 'taxa' or 'all'.")
  # https://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/
  p = ggplot2::ggplot(to_plot, ggplot2::aes(Var1, Var2, fill = scaled)) +
    ggplot2::geom_tile(color = 'black', size = 0.25) + 
    ggplot2::scale_fill_gradientn(colours = c('blue', 'white', 'red'), 
                                  values = c(min(to_plot$scaled), 
                                             mean(to_plot$scaled), 
                                             max(to_plot$scaled))) +
    ggplot2::xlab('') + ggplot2::ylab('') +
    ggplot2::theme(legend.position = 'none', 
                   axis.ticks = ggplot2::element_blank(), 
                   axis.text.x = ggplot2::element_text(angle = 90, 
                                                       hjust = 1, vjust = 0.25), 
                   axis.text = ggplot2::element_text(color = 'gray20'))+
    ggplot2::geom_text(data = to_plot, ggplot2::aes(label = value), size = 3) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) + 
    ggplot2::scale_y_discrete(expand = c(0, 0))
  p
}