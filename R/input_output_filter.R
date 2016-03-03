# mctoolsr

###################################
# LOADING/FILTER SAMPLE FUNCTIONS #
###################################

#' @title Load a taxa table for use with mctoolsr
#' @description Load in a taxa table (aka. an OTU table) and a corresponding 
#'  mapping file with metadata values. The samples in the loaded taxa table 
#'  and mapping file will be in the same order and only samples in both will 
#'  be loaded. The function can optionally filter samples of a specific type 
#'  based on the mapping file. This can also be done later via the filter_data() 
#'  function.
#' @param tab_fp Taxa table filepath.
#' @param map_fp Metadata mapping filepath.
#' @param filter_cat The map_fp header string for the factor you would like 
#'  to use to filter samples.
#' @param filter_vals The values within the filter category (vector or single 
#'  value) you would like to use to remove samples from the imported data.
#' @param keep_vals Alternatively, keep only samples represented by these 
#'  values.
#' @return A list variable with (1) the loaded taxa table, and (2) the loaded 
#'  mapping file.
#' @examples 
#' \dontrun{
#' load_taxa_table("filepath_to_OTU_table.txt", "filepath_to_mapping_file.txt",
#'   "sample_type", filter_vals = "blank")
#' }
load_taxa_table = function(tab_fp, map_fp, filter_cat, filter_vals, keep_vals){
  # load data
  if(tools::file_ext(tab_fp) == 'biom'){
    data_b = biom::read_biom(tab_fp)
    data = as.data.frame(as.matrix(biom::biom_data(data_b)))
    data_taxonomy = .compile_taxonomy(data_b)
  }
  else if(tools::file_ext(tab_fp) == 'txt'){
    if(readChar(tab_fp, nchars = 4) == "#OTU"){
      data = read.table(tab_fp, sep = '\t', comment.char = '', header = TRUE, 
                        check.names = FALSE, row.names = 1, quote = "\"")
    } else {
      data = read.table(tab_fp, sep='\t', skip=1, comment.char = '', 
                        header = TRUE, check.names = FALSE, row.names = 1, 
                        quote = "\"")
    }
    if(names(data)[ncol(data)] == 'taxonomy'){
      data_taxonomy = .parse_taxonomy(data$taxonomy)
      if(!is.null(data_taxonomy)) {
        row.names(data_taxonomy) = row.names(data)
      }
      data$taxonomy = NULL
    } else data_taxonomy = NULL
  }
  else stop('Input file must be either biom (.biom) or tab-delimited (.txt) format.')
  # import mapping file
  map = tryCatch(read.table(map_fp, sep = '\t', comment.char = '', header = TRUE, 
                            check.names = FALSE, row.names = 1, quote = "\""),
                 error = function(c) conditionMessage(c),
                 warning = function(c) {
                   c$message = paste0("Error loading mapping file. Please ",
                            "check that there are no quotation marks.")
                   stop(c)
                   },
                 message = function(c) conditionMessage(c)
                 )
  if(class(map) != 'data.frame') {
    warning(paste0('Problem loading mapping file. Note that the mapping file ',
                   'should have more than one metadata column. Check that ', 
                   'All rows have the same number of columns and that there ', 
                   'are no duplicate sample IDs. It could also be that the ', 
                   'specified filepath is not correct.'))}
  # optionally, subset data
    # cant subset if trying to filter out certain values and keep certain values
    # use one or the other
  if(!missing(filter_cat)){
    map_f = .filt_map(map, filter_cat, filter_vals, keep_vals)
  } else {
    map_f = map
  }
  # match up data from dissimilarity matrix with mapping file
  matched_data = .match_data_components(data, map_f, data_taxonomy)
  message(paste0(nrow(matched_data$map_loaded), ' samples loaded'))
  matched_data
}

#' @title Load a dissimilarity matrix for use with mctoolsr
#' @description Load in a dissimilarity matrix and a corresponding metadata
#'  mapping file
#' @param dm_fp Dissimilarity matrix filepath (tab-delimited text).
#' @param map_fp Metadata mapping filepath.
#' @param filter_cat The map_fp header string for the factor you would like 
#'  to use to filter samples.
#' @param filter_vals The values within the filter category (vector or single 
#'  value) you would like to use to remove samples from the imported data.
#' @param keep_vals Alternatively, keep only samples represented by these 
#'  values.
#' @return A list variable with (1) the loaded dissimilarity matrix, and (2) 
#'  the loaded mapping file.
#' @examples 
#' \dontrun{
#' load_dm("filepath_to_dissim_matrix.txt", "filepath_to_mapping_file.txt",
#'   "sample_type", filter_vals = "blank")
#' }
load_dm = function(dm_fp, map_fp, filter_cat, filter_vals, keep_vals){
  dm = read.table(dm_fp, sep='\t', comment.char='', header=T, check.names=F, 
                  row.names=1)
  map = read.table(map_fp, sep='\t', comment.char='', header=T, check.names=F, 
                   row.names=1)
  # optionally, subset data
  # cant subset if trying to filter out certain values and keep certain values
  # use one or the other
  if(!missing(filter_cat)){
    map_f = .filt_map(map, filter_cat, filter_vals, keep_vals)
  } else map_f = map
  # match up data from dissimilarity matrix with mapping file
  samplesToUse = intersect(names(dm), row.names(map_f))
  dm.use = as.dist(dm[match(samplesToUse,names(dm)), 
                      match(samplesToUse, names(dm))])
  map.use = map_f[match(samplesToUse, row.names(map_f)), ]
  # output
  list(dm_loaded = dm.use, map_loaded = map.use)
}

#' @title Load two dissimilarity matrices for use with mctoolsr
#' @description Load in two dissimilarity matrices and a corresponding metadata
#'  mapping file. Useful for Mantel tests when dms are already generated.
#' @param dm1_fp Dissimilarity matrix filepath (tab-delimited text).
#' @param dm2_fp Dissimilarity matrix filepath (tab-delimited text).
#' @param map_fp Metadata mapping filepath.
#' @param filter_cat The map_fp header string for the factor you would like 
#'  to use to filter samples.
#' @param filter_vals The values within the filter category (vector or single 
#'  value) you would like to use to remove samples from the imported data.
#' @param keep_vals Alternatively, keep only samples represented by these 
#'  values.
#' @return A list variable with (1) the loaded dissimilarity matrix 1, (2) the
#'  loaded dissimilarity matrix 2, and (3) the loaded mapping file.
#' @examples 
#' \dontrun{
#' load_dm("filepath_to_dissim_matrix1.txt", "filepath_to_dissim_matrix1.txt", 
#'   "filepath_to_mapping_file.txt", "sample_type", filter_vals = "blank")
#' }
load_2_dms = function(dm1_fp, dm2_fp, map_fp, filter_cat, filter_vals, keep_vals){
  dm1 = read.table(dm1_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  dm2 = read.table(dm2_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  map = read.table(map_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  # optionally, subset data
  # cant subset if trying to filter out certain values and keep certain values
  # use one or the other
  if(!missing(filter_cat)){
    map_f = .filt_map(map, filter_cat, filter_vals, keep_vals)
  } else map_f = map
  # match up data from dissimilarity matrix with mapping file
  samplesToUse = intersect(intersect(names(dm1), row.names(map_f)), names(dm2))
  dm1.use = as.dist(dm1[match(samplesToUse,names(dm1)), match(samplesToUse,names(dm1))])
  dm2.use = as.dist(dm2[match(samplesToUse,names(dm2)), match(samplesToUse,names(dm2))])
  map.use = map_f[match(samplesToUse,row.names(map_f)), ]
  # output
  list(dm1_loaded = dm1.use, dm2_loaded = dm2.use, map_loaded = map.use)
}

#' @title Filter Samples from Dataset
#' @description Filter out or keep particular samples in a dataset based on 
#'  contextual metadata.
#' @param input The input dataset as loaded by \code{load_taxa_table()}.
#' @param filter_cat The map_fp header string for the factor you would like 
#'  to use to filter samples.
#' @param filter_vals The values within the filter category (vector or single 
#'  value) you would like to use to remove samples from the imported data.
#' @param keep_vals Alternatively, keep only samples represented by these 
#'  values.
#' @return A list variable with (1) the loaded taxa table, (2) 
#'  the loaded mapping file, and optionally (3) the loaded taxonomy information.
#' @examples 
#' \dontrun{
#' ex_in_filt = filter_data(input = "example_input", filter_cat = "Sample_type", 
#'                          filter_vals = c("mushrooms", "strawberries"))
#' }
filter_data = function(input, filter_cat, filter_vals, keep_vals){
  if(!missing(filter_cat)){
    map_f = .filt_map(input$map_loaded, filter_cat, filter_vals, keep_vals)
  } else map_f = input$map_loaded
  # match up data from dissimilarity matrix with mapping file
  if('taxonomy_loaded' %in% names(input)){
    matched_data = .match_data_components(input$data_loaded, map_f, 
                                          input$taxonomy_loaded)
    message(paste0(nrow(matched_data$map_loaded), ' samples remaining'))
    matched_data
  } else {
    matched_data = .match_data_components(input$data_loaded, map_f, NULL)
    message(paste0(nrow(matched_data$map_loaded), ' samples remaining'))
    matched_data
  }
}

#' @title Filter Samples from a Dataset based on number of sequences
#' @param input The input dataset as loaded by \code{load_taxa_table()}.
#' @param min_seqs A sample must have at minimum this number of sequences 
#'  in order to be retained.
#' @return A list variable with (1) the loaded taxa table, (2) 
#'  the loaded mapping file, and optionally (3) the loaded taxonomy information.
filter_samples_by_counts = function(input, min_seqs) {
  data_filt = input$data_loaded[, colSums(input$data_loaded) >= min_seqs]
  # match up data from dissimilarity matrix with mapping file
  if('taxonomy_loaded' %in% names(input)){
    matched_data = .match_data_components(data_filt, input$map_loaded, 
                                          input$taxonomy_loaded)
    message(paste0(nrow(matched_data$map_loaded), ' samples remaining'))
    matched_data
  } else {
    matched_data = .match_data_components(data_filt, input$map_loaded, NULL)
    message(paste0(nrow(matched_data$map_loaded), ' samples remaining'))
    matched_data
  }
}

#' @title Filter Samples from Dissimilarity Matrix
#' @description Filter out or keep particular samples in a dissimilarity matrix 
#'  based on contextual metadata.
#' @param input_dm The input dissimilarity matrix with corresponding mapping 
#'  file as generated by \code{load_dm()} or \code{calc_mean_dissimilarities()} 
#'  with the option to produce a resulting metadata map.
#' @param filter_cat The metadata map header string for the factor you would 
#'  like to use to filter samples.
#' @param filter_vals The values within the filter category (vector or single 
#'  value) you would like to use to remove samples from the dissimilarity 
#'  matrix and metadata map.
#' @param keep_vals Alternatively, keep only samples represented by these 
#'  values.
#' @return A list variable with (1) the loaded dissimilarity matrix and (2) 
#'  the loaded mapping file.
filter_dm = function(input_dm, filter_cat, filter_vals, keep_vals){
#   if(!missing(filter_vals)) {
#     if(!missing(keep_vals)) {
#       stop('Can only use one of "filter_vals" or "keep_vals" at a time.')
#     }
#     toKeep = !input_dm$map_loaded[, filter_cat] %in% keep_vals
#   }
#   else if(!missing(keep_vals)) {
#     toKeep = input_dm$map_loaded[, filter_cat] %in% keep_vals
#   }
  map_filt = .filt_map(input_dm$map_loaded, filter_cat, filter_vals, keep_vals)
  # match up data from dissimilarity matrix with mapping file
  dm = as.matrix(input_dm$dm_loaded)
  samplesToUse = intersect(colnames(dm), row.names(map_filt))
  dm_use = as.dist(dm[match(samplesToUse, colnames(dm)), 
                      match(samplesToUse, colnames(dm))])
  map_use = map_filt[match(samplesToUse, row.names(map_filt)), ]
  # output
  list(dm_loaded = dm_use, map_loaded = map_use)
}

#' @title Match up samples from two datasets
#' @description Function to match up sample order from two datasets that contain
#'  some overlapping sample IDs. Sample IDs that are not present in both
#'  datasets will be dropped. The output is a list containing the two filtered
#'  datasets in the same order as they were input.
#' @param ds1 The first dataset as loaded by \code{load_taxa_table()}
#' @param ds2 The second dataset.
#' @return A list variable with the matched ds1 as the first element and ds2
#'  as the second element
match_datasets = function(ds1, ds2){
  common_samples = intersect(names(ds1$data_loaded), names(ds2$data_loaded))
  ds1$map_loaded$common_sample = row.names(ds1$map_loaded) %in% common_samples
  ds1_filt = filter_data(ds1, filter_cat = 'common_sample', filter_vals = FALSE)
  ds1_filt$data_loaded = ds1_filt$data_loaded[, 
                                              match(common_samples, 
                                                    names(ds1_filt$data_loaded))]
  ds1_filt$map_loaded = ds1_filt$map_loaded[match(common_samples, 
                                                  row.names(ds1_filt$map_loaded)), 
                                            ]
  ds2$map_loaded$common_sample = row.names(ds2$map_loaded) %in% common_samples
  ds2_filt = filter_data(ds2, 'common_sample', FALSE)
  ds2_filt$data_loaded = ds2_filt$data_loaded[, 
                                              match(common_samples, 
                                                    names(ds2_filt$data_loaded))]
  ds2_filt$map_loaded = ds2_filt$map_loaded[match(common_samples, 
                                                  row.names(ds2_filt$map_loaded)), 
                                            ]
  list(ds1 = ds1_filt, ds2 = ds2_filt)
}

#' @title Export an OTU table as a text file
#' @description A convenient way to export a loaded OTU table as a text file. 
#'  Taxonomy strings will appear in the right most column. This is also a good
#'  way to save an OTU table to be loaded later. Output format is tab-delimited.
#' @param input The input dataset as loaded by \code{load_taxa_table()} or
#'  an otu table of class \code{data.frame}.
#' @param out_fp The output filepath.
#' @param map_fp (OPTIONAL) The metadata map output filepath if you want to 
#'  write it to file.
export_otu_table = function(input, out_fp, map_fp){
  if(class(input) == 'list') {
    table = input$data_loaded
    if(!is.null(input$taxonomy_loaded)) {
      taxonomy = apply(input$taxonomy_loaded, 1, paste, collapse = '; ')
      out_tab = data.frame(OTU_ID = row.names(table), table, taxonomy, 
                           check.names = FALSE)
    } else {
      out_tab = data.frame(OTU_ID = row.names(table), table, 
                           check.names = FALSE)
    }
  } else if(class(input) == 'data.frame') {
    table = input
    if(!missing(map_fp)){
      stop(paste0('Cannot write a metadata map unless input is a list that 
                  contains the metadata.'))
    }
    out_tab = data.frame(OTU_ID = row.names(table), table, 
                         check.names = FALSE)
  }
  names(out_tab)[1] = '#OTU ID'
  write('#Exported from mctoolsr', out_fp)
  suppressWarnings(write.table(out_tab, out_fp, sep = '\t', row.names = FALSE, 
                               append = TRUE))
  if(!missing(map_fp)){
    if(class(map_fp) != 'character') stop('map_fp must be a valid filepath.')
    map_file = data.frame(row.names(input$map_loaded), input$map_loaded, 
                          check.names = FALSE)
    names(map_file)[1] = '#Sample ID'
    write.table(map_file, map_fp, sep = '\t', row.names = FALSE)
  }
}



