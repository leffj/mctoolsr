# mctoolsr

#############
# UTILITIES #
#############

#' @keywords internal
.filt_map = function(map, filter_cat, filter_vals, keep_vals){
  if(!missing(filter_vals) & !missing(keep_vals)){
    stop('Can only handle filter_vals or keep_vals, not both.')
  }
  if(!filter_cat %in% names(map)){
    stop('filter_cat not found in mapping file headers. Check spelling.')
  }
  # filter out values within mapping category
  else if(!missing(filter_cat) & !missing(filter_vals)){
    map_f = map[!map[, filter_cat] %in% filter_vals, , drop = FALSE]
    map_f = droplevels(map_f)
    if(nrow(map_f) == 0){
      stop('All rows filtered out. Check spelling of filter parameters.')
    }
  }
  # keep certain values with mapping category
  else if(!missing(filter_cat) & !missing(keep_vals)){
    map_f = map[map[,filter_cat] %in% keep_vals, , drop = FALSE]
    map_f = droplevels(map_f)
    if(nrow(map_f) == 0){
      stop('All rows filtered out. Check spelling of filter parameters.')
    }
  }
  map_f
}

#' @keywords internal
.summarize_map = function(metadata_map, summarize_by_factor) {
  .smry_fun = function(x) {
    if (is.numeric(x)) {
      mean(x)
    } else {
      if (length(unique(x)) == 1) {
        unique(x)
      } else
        NA
    }
  }
  # change row names for NA values in summarize_by_factor with warning
  na_idxs = is.na(metadata_map[, summarize_by_factor])
  if (sum(na_idxs) > 0) {
    warning(
      paste0(
        'NA values present in "summarize_by_factor". NAs will be ',
        'referred to as "NO_VALUE".'
      )
    )
    vec = as.character(metadata_map[, summarize_by_factor])
    vec[na_idxs] = 'NO_VALUE'
    metadata_map[, summarize_by_factor] = factor(vec)
  }
  mean_map = NULL
  for (i in seq_along(metadata_map)) {
    name = colnames(metadata_map)[i]
    if (class(unlist(metadata_map[i])) == 'factor') {
      x = as.character(unlist(metadata_map[i]))
    } else {
      x = unlist(metadata_map[i])
    }
    result = tapply(x, metadata_map[, summarize_by_factor], .smry_fun)
    newnames = c(colnames(mean_map), name)
    mean_map = cbind(mean_map, result)
    colnames(mean_map) = newnames
  }
  # mean_map = dplyr::summarise_each(dplyr::group_by_(metadata_map,
  #                                                   summarize_by_factor),
  #                                  dplyr::funs(.smry_fun))
  mean_map = as.data.frame(as.matrix(mean_map))
  # row.names(mean_map) = mean_map[, summarize_by_factor]
  mean_map
}

#' @keywords internal
.match_data_components = function(tax_table, metadata_map, taxonomy){
  samplesToUse = intersect(names(tax_table), row.names(metadata_map))
  tax_table.use = tax_table[, match(samplesToUse, names(tax_table)), 
                            drop = FALSE]
  tax_table.use = tax_table.use[rowSums(tax_table.use) != 0, , drop = FALSE]
  map.use = metadata_map[match(samplesToUse, row.names(metadata_map)), , 
                         drop = FALSE]
  map.use = droplevels(map.use)
  if(!missing('taxonomy') & !is.null(taxonomy)) {
    taxonomy.use = taxonomy[match(row.names(tax_table.use), 
                                  row.names(taxonomy)), ]
    taxonomy.use = droplevels(taxonomy.use)
    list(data_loaded = tax_table.use, map_loaded = map.use, 
         taxonomy_loaded = taxonomy.use)
  } else {
    list(data_loaded = tax_table.use, map_loaded = map.use)
  }
}

#' @title Rename samples in an mctoolsr dataset
#' @description Rename the samples by substituting column names in the taxa 
#'  table and row names in the metadata map with values from a metadata map
#'  column that you specify. Note that all values in the metadata map column 
#'  must be unique.
#' @param input The input dataset as loaded by \code{\link{load_taxa_table}}.
#' @param name_header The header value in the metadata map that will be used
#'  to rename the samples.
#' @concept Taxa table manipulation
#' @examples 
#' fruits_veggies$map_loaded$alt_id =
#' paste0('alt', 1:nrow(fruits_veggies$map_loaded))
#' rename_samples(fruits_veggies, 'alt_id')
rename_samples = function(input, name_header) {
  colnames(input$data_loaded) = input$map_loaded[, name_header]
  row.names(input$map_loaded) = input$map_loaded[, name_header]
  input
}
