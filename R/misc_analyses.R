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