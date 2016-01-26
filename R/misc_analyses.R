# mctoolsr

#################
# MISC ANALYSES #
#################

#' @title Determine taxa that are common (core) across sample types
#' @description This function returns a data frame with taxa that meet 
#'  proportion user thresholds and are determined to be common across multiple 
#'  sample types.
#' @param input The input dataset as loaded by \code{load_taxa_table()}.
#' @param type_header Mapping file header (in quotation marks) of factor for 
#'  which you are testing for differences.
#' @param prop_types The proportion of sample types in which a taxon must
#'  be present to be considered core.
#' @param prop_reps The proportion of replicate samples within a sample type
#'  in which a taxon must be present to be considered core.
#' @return A data frame with each core taxon as a row. The proportion of sample
#'  types where the taxon was present, the mean count of each 
#'  sample type, and the corresponding taxonomy (if provided), are also 
#'  returned.
core_taxa = function(input, type_header, prop_types = 1, prop_reps = 0.5) {
  # indicate taxa that are present within ≥ <prop_reps> of replicates
  df = as.data.frame(t(input$data_loaded))
  df$type = input$map_loaded[, type_header]
  df_m = reshape2::melt(df, id.vars = 'type', variable.name = 'OTU_ID')
  gt_0 = function(x) sum(x > 0)
  by_type = dplyr::summarise(dplyr::group_by(df_m, type, OTU_ID), 
                             reps = length(value), reps_obs = gt_0(value), 
                             mean_val = mean(value))
  by_type = dplyr::mutate(by_type, prop_reps_obs = reps_obs / reps, 
                          obs = prop_reps_obs >= prop_reps)
  # indicate taxa that are core across sample types
  # calc proportion of sample types that observed each OTU
  by_OTU = dplyr::summarise(dplyr::group_by(by_type, OTU_ID), 
                            prop_sample_types = sum(obs) / length(obs))
  # add on mean values
  means = reshape2::dcast(by_type, OTU_ID ~ type, value.var = 'mean_val')
  if(!identical(by_OTU$OTU_ID, means$OTU_ID)) stop('Unknown error.')
  by_OTU = cbind(by_OTU, means[, 2:ncol(means)])
  # return core with taxonomy (if provided)
  core = as.data.frame(dplyr::filter(by_OTU, prop_sample_types >= prop_types))
  if(is.null(input$taxonomy_loaded)) core else {
    tax = input$taxonomy_loaded[match(core$OTU_ID, 
                                      row.names(input$taxonomy_loaded)), ]
    core$taxonomy = apply(tax, 1, paste0, collapse = '; ')
    core
  }
}

#' @title Calculate mean proportions of shared taxa between paired samples
#' @description This function will calculate the mean proportion of taxa 
#'  (i.e. OTUs) that are shared between pairs of samples in a given data set.
#'  It has the option to calculate this value for only pairs where each sample
#'  is of a different, specified type. For example, this is useful if you want
#'  to calculate proportion of taxa that are shared between red and blue 
#'  samples only and not red and red or red and green, etc.  
#' @param input Data from load_taxa_table()
#' @param type_header (Optional) The header label from the metadata map that
#'  has the labels which indicate the types of samples you want to calculate 
#'  proportions for.
#' @param sample_types (Optional) Use with type_header to indicate which two 
#'  sample types for which you want to calculate proportions. Supply a character
#'  vector of length 2.
#' @param within_cat (Optional) Specify a header label from the metadata map 
#'  that tells the function to restrict pairs to within levels of this factor.
calc_prop_shared_taxa = function(input, type_header, sample_types, within_cat) {
  calc_prop_shared_main = function(input, type_header, sample_types) {
    # get all pairs of samples
    pairs = as.data.frame(t(combn(colnames(input$data_loaded), 2)))
    # pairs of certain types
    if(!missing(type_header) & !missing(sample_types)) {
      head(pairs)
      pairs$S1_type = 
        input$map_loaded[match(pairs[, 1], row.names(input$map_loaded)), 
                         type_header]
      pairs$S2_type = 
        input$map_loaded[match(pairs[, 2], row.names(input$map_loaded)), 
                         type_header]
      pairs = filter(pairs, S1_type != S2_type, S1_type %in% sample_types, 
                     S2_type %in% sample_types)
      # check if missing sample type(s)
      mis_sts = sample_types[! sample_types %in% c(levels(pairs$S1_type), levels(pairs$S2_type))]
      if(length(mis_sts) > 0) {
        stop(paste("Missing sample type: ", mis_sts, " "))
      }
    } else if(missing(type_header) & missing(sample_types)) {
      # all pairs in data set 
    } else {
      stop("Must use 'cat_header' and 'sample_types' together.")
    }
    # function to calc proportion shared for each pair
    .prop_shared_pair = function(pair, tax_table) {
      pair_data = tax_table[, c(pair[1], pair[2])]
      pair_data_f = pair_data[rowSums(pair_data) > 0, ]  # filt 0s
      no_total = nrow(pair_data_f)
      no_shared = sum(pair_data_f[, 1] * pair_data_f[, 2] > 0)
      prop_shared = no_shared / no_total
      prop_shared
    }
    pairs$prop_shared = apply(pairs, 1, .prop_shared_pair, input$data_loaded)
    mean(pairs$prop_shared)
  }
  if(missing(within_cat)) {
    calc_prop_shared_main(input, type_header, sample_types)
  } else {
    props = c()
    for(cat_lev in unique(input$map_loaded[, within_cat])) {
      tmp_filt = suppressMessages(filter_data(input, within_cat, 
                                              keep_vals = cat_lev))
      tmp_prop = calc_prop_shared_main(tmp_filt, type_header, sample_types)
      props = c(props, tmp_prop)
    }
    mean(props)
  }
}