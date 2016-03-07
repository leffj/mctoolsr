# mctoolsr

#################
# NORMALIZATION #
#################

#' @title Rarefy samples in a taxa table
#' @description Rarefying means to randomly select a specified number of 
#'  observations from each sample. This approach is used to normalize for 
#'  differences in the number of observations (i.e. sequencing depth) across
#'  samples. Samples with fewer observations (sequences) than the specified
#'  depth will be removed.
#' @param input The input dataset as loaded by \code{load_taxa_table()}.
#' @param depth The number of observations to randomly select for each sample.
#' @return Will return a dataset object similar to the input.
single_rarefy = function(input, depth) {
  data_filt_samples = input$data_loaded[, colSums(input$data_loaded) >= depth]
  data_rar = as.data.frame(t(vegan::rrarefy(t(data_filt_samples), depth)))
  matched_data = .match_data_components(data_rar, input$map_loaded, 
                                        input$taxonomy_loaded)
  message(paste0(nrow(matched_data$map_loaded), ' samples remaining'))
  matched_data
}


#' @title Calculate mean taxa values across a specified factor
#' @param input The input dataset as loaded by \code{load_taxa_table()} or
#'  an otu table of class \code{data.frame}.
#' @param metadata_map [Optional]. The metadata mapping data frame. Required if
#'  input is a \code{data.frame}.
#' @param summarize_by_factor Category in mapping file to summarize by.
#' @return If input is a list, returns a list with a taxon table (data_loaded) 
#'  and a mapping data frame (map_loaded). It will automatically return 
#'  taxonomy in the list if provided in the input.
calc_taxa_means = function(input, summarize_by_factor, metadata_map) {
  .calc_tt_means = function(table, metadata_map, summarize_by_factor){
    as.data.frame(t(apply(table, 1, function(x) {
      tapply(x, metadata_map[, summarize_by_factor], mean)
    })))
  }
  if(class(input) == 'list') {
    tt_means = .calc_tt_means(input$data_loaded, input$map_loaded, 
                              summarize_by_factor)
    mean_map = .summarize_map(input$map_loaded, summarize_by_factor)
    map_loaded = mean_map[match(colnames(tt_means), row.names(mean_map)), ]
    output = list(data_loaded = tt_means, map_loaded = map_loaded)
    if(!is.null(input$taxonomy_loaded)) {
      c(output, list(taxonomy_loaded = input$taxonomy_loaded))
    } else output
  } else if(class(input) == 'data.frame') {
    taxa_table_means = .calc_tt_means(input, metadata_map, summarize_by_factor)
    mean_map = .summarize_map(metadata_map, summarize_by_factor)
    .match_data_components(taxa_table_means, mean_map, NULL)
  } else stop('input is of an incorrect variable type.')
}