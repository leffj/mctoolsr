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
    map.f = map[!map[, filter_cat] %in% filter_vals, ]
    map.f = droplevels(map.f)
    if(nrow(map.f) == 0){
      stop('All rows filtered out. Check spelling of filter parameters.')
    }
  }
  # keep certain values with mapping category
  else if(!missing(filter_cat) & !missing(keep_vals)){
    map.f = map[map[,filter_cat] %in% keep_vals, ]
    map.f = droplevels(map.f)
    if(nrow(map.f) == 0){
      stop('All rows filtered out. Check spelling of filter parameters.')
    }
  }
  map.f
}

#' @keywords internal
.match_data_components = function(tax_table, metadata_map, taxonomy){
  samplesToUse = intersect(names(tax_table), row.names(map))
  tax_table.use = tax_table[, match(samplesToUse, names(tax_table))]
  tax_table.use = tax_table.use[rowSums(tax_table.use) != 0, ]
  map.use = map[match(samplesToUse, row.names(map)),]
  map.use = droplevels(map.use)
  if(!missing('taxonomy')) {
    taxonomy.use = taxonomy[match(row.names(tax_table.use), row.names(taxonomy)), ]
    taxonomy.use = droplevels(taxonomy.use)
    list(data_loaded = tax_table.use, map_loaded = map.use, taxonomy_loaded = taxonomy.use)
  } else {
    list(data_loaded = tax_table.use, map_loaded = map.use)
  }
}

#' @title Convert Taxon Table to Relative Abundances
#' @description Convert taxon table or taxon table in data input to relative 
#' abundances
#' @param input Either a dataset as generated using 'load_taxon_table' which
#' includes a mapping file or an individual taxon table
convert_to_relative_abundances = function(input){
  if('map_loaded' %in% names(input)){
    seq_cts = colSums(input$data_loaded)
    rel_abund_table = as.data.frame(t(apply(input$data_loaded, 1, 
                                            function(x) x / seq_cts)))
    list(data_loaded = rel_abund_table, map_loaded = input$map_loaded, 
         taxonomy_loaded = input$taxonomy_loaded)
  } else {
    seq_cts = colSums(input, na.rm = TRUE)
    rel_abund_table = as.data.frame(t(apply(input, 1, function(x) x / seq_cts)))
    rel_abund_table
  }
  
}