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
