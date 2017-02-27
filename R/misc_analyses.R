# mctoolsr

#################
# MISC ANALYSES #
#################

#' @title Determine taxa that are common (core) across sample types
#' @description This function returns a data frame with taxa that are 
#'  determined to be common across multiple sample types based on default or
#'  user defined thresholds.
#' @param input The input dataset as loaded by \code{\link{load_taxa_table}}.
#' @param type_header Mapping file header (in quotation marks) of factor for
#'  which you are testing for differences.
#' @param prop_types [OPTIONAL] The proportion of sample types in which a taxon must
#'  be present to be considered core.
#' @param prop_reps [OPTIONAL] The proportion of replicate samples within a sample type
#'  in which a taxon must be present to be considered core.
#' @return A data frame with each core taxon as a row. The proportion of sample
#'  types where the taxon was present, the mean count of each
#'  sample type, and the corresponding taxonomy (if provided), are also
#'  returned.
#' @concept Misc analyses
#' @examples 
#' # taxa core to fruits and vegetables
#' core_taxa(fruits_veggies, type_header = 'Sample_type')
#' # relax the proportion of produce/sample types
#' core_taxa(fruits_veggies, type_header = 'Sample_type', prop_types = 0.6)
core_taxa = function(input, type_header, prop_types = 1, prop_reps = 0.5) {
  # indicate taxa that are present within ≥ <prop_reps> of replicates
  df = as.data.frame(t(input$data_loaded))
  df$type = input$map_loaded[, type_header]
  df_m = reshape2::melt(df, id.vars = 'type', variable.name = 'OTU_ID')
  gt_0 = function(x) sum(x > 0)
  by_type = dplyr::summarise_(
    dplyr::group_by_(df_m, "type", "OTU_ID"),
    reps = ~ length(value), reps_obs = ~ gt_0(value),
    mean_val = ~ mean(value)
  )
  by_type = dplyr::mutate_(by_type, prop_reps_obs = ~ reps_obs / reps,
                          obs = ~ prop_reps_obs >= prop_reps)
  # indicate taxa that are core across sample types
  # calc proportion of sample types that observed each OTU
  by_OTU = dplyr::summarise_(dplyr::group_by_(by_type, "OTU_ID"),
                            prop_sample_types = ~ sum(obs) / length(obs))
  # add on mean values
  means = reshape2::dcast(by_type, OTU_ID ~ type, value.var = 'mean_val')
  if (!identical(by_OTU$OTU_ID, means$OTU_ID))
    stop('Unknown error.')
  by_OTU = cbind(by_OTU, means[, 2:ncol(means)])
  # return core with taxonomy (if provided)
  core = as.data.frame(dplyr::filter_(by_OTU, ~ prop_sample_types >= prop_types))
  if (is.null(input$taxonomy_loaded))
    core
  else {
    tax = input$taxonomy_loaded[match(core$OTU_ID,
                                      row.names(input$taxonomy_loaded)),]
    core$taxonomy = apply(tax, 1, paste0, collapse = '; ')
    core
  }
}

#' @title Calculate mean proportions of shared taxa between paired samples
#' @description This function will calculate the mean proportion of taxa (i.e.
#'   OTUs) that are shared between pairs of samples in a given data set. It has
#'   the option to calculate this value for only pairs where each sample is of a
#'   different, specified type. For example, this is useful if you want to
#'   calculate proportion of taxa that are shared between red and blue samples
#'   only and not red and red or red and green, etc.
#' @param input The input dataset as loaded by \code{\link{load_taxa_table}}.
#' @param type_header [OPTIONAL] The header label from the metadata map that has
#'   the labels which indicate the types of samples you want to calculate 
#'   proportions for.
#' @param sample_types [OPTIONAL] Use with \code{type_header} to indicate which
#'   two sample types for which you want to calculate proportions. Supply a
#'   character vector of length 2.
#' @param within_cat [OPTIONAL] Specify a header label from the metadata map 
#'   that tells the function to restrict pairs to within levels of this factor. 
#'   If included, will return a data frame instead of a vector.
#' @details If only \code{input} is provided, the proportion is calculated
#'   across all samples.
#' @seealso \link[mctoolsr]{calc_prop_taxa_from_sample_type}
#' @concept Misc analyses
#' @examples 
#' # proportion of taxa shared between spinach and strawberries
#' calc_prop_shared_taxa(
#'   fruits_veggies, type_header = 'Sample_type',
#'   sample_types = c('Spinach', 'Strawberries')
#' )
#' # proportion of taxa shared between spinach and mushrooms
#' calc_prop_shared_taxa(
#'   fruits_veggies, type_header = 'Sample_type', 
#'   sample_types = c('Strawberries', 'Mushrooms')
#' )
#' # proportion of taxa shared between spinach and strawberries within each
#' # farming type (greater overlap within conventional)
#' calc_prop_shared_taxa(
#'   fruits_veggies, type_header = 'Sample_type', 
#'   sample_types = c('Spinach', 'Strawberries'), within_cat = 'Farm_type'
#' )
calc_prop_shared_taxa = function(input, type_header, sample_types, within_cat) {
  calc_prop_shared_main = function(input, type_header, sample_types) {
    # get all pairs of samples
    pairs = as.data.frame(t(combn(colnames(
      input$data_loaded
    ), 2)))
    # pairs of certain types
    if (!missing(type_header) & !missing(sample_types)) {
      pairs$S1_type =
        input$map_loaded[match(pairs[, 1], row.names(input$map_loaded)),
                         type_header]
      pairs$S2_type =
        input$map_loaded[match(pairs[, 2], row.names(input$map_loaded)),
                         type_header]
      conds = pairs$S1_type != pairs$S2_type &
        pairs$S1_type %in% sample_types & pairs$S2_type %in% sample_types
      pairs = pairs[conds,]
      # check if missing sample type(s)
      mis_sts = sample_types[!sample_types %in% c(levels(pairs$S1_type),
                                                  levels(pairs$S2_type))]
      if (length(mis_sts) > 0) {
        stop(paste("Missing sample type: ", mis_sts, " "))
      }
    } else if (missing(type_header) & missing(sample_types)) {
      # all pairs in data set
    } else {
      stop("Must use 'cat_header' and 'sample_types' together.")
    }
    # function to calc proportion shared for each pair
    .prop_shared_pair = function(pair, tax_table) {
      pair_data = tax_table[, c(pair[1], pair[2])]
      pair_data_f = pair_data[rowSums(pair_data) > 0,]  # filt 0s
      no_total = nrow(pair_data_f)
      no_shared = sum(pair_data_f[, 1] * pair_data_f[, 2] > 0)
      prop_shared = no_shared / no_total
      prop_shared
    }
    pairs$prop_shared = apply(pairs, 1, .prop_shared_pair, input$data_loaded)
    mean(pairs$prop_shared)
  }
  if (missing(within_cat)) {
    calc_prop_shared_main(input, type_header, sample_types)
  } else {
    props = c()
    cat_levs = unique(input$map_loaded[, within_cat])
    for (cat_lev in cat_levs) {
      tmp_filt = suppressMessages(filter_data(input, within_cat,
                                              keep_vals = cat_lev))
      tmp_prop = calc_prop_shared_main(tmp_filt, type_header, sample_types)
      props = c(props, tmp_prop)
    }
    data.frame(cat = cat_levs, prop = props)
  }
}

#' @title Calculate the proportion of taxa in a set of samples that are
#'  also observed in another sample type
#' @description Calculate the mean proportion of taxa (i.e. OTUs) in a set of
#'   samples that are also observed in another sample type. One use of this
#'   function is to get a rough sense of how likely a potential source is
#'   contributing to the compostion of another sample type.
#' @param input The input dataset as loaded by \code{\link{load_taxa_table}}.
#' @param type_header The header label from the metadata map that has the
#'  labels which indicate the types of samples you want to use.
#' @param primary_type The sample type for which you would like to calculate
#'  the proportion of taxa that are also observed in the \code{source_type}.
#' @param source_type The sample type to use as the source type.
#' @param within_cat [OPTIONAL] Specify a header label from the metadata map
#'  that tells the function to restrict pairs to within levels of this factor. 
#'  If included, will return a data frame instead of a vector.
#' @seealso \link[mctoolsr]{calc_prop_shared_taxa}
#' @concept Misc analyses
#' @examples 
#' # proportion of taxa on strawberries that were also found on spinach
#' calc_prop_taxa_from_sample_type(
#'   fruits_veggies, type_header = 'Sample_type', primary_type = 'Strawberries',
#'   source_type = 'Spinach'
#' )
#' # proportion of taxa on strawberries that were also found on spinach within
#' # each farming type
#' calc_prop_taxa_from_sample_type(
#'   fruits_veggies, type_header = 'Sample_type', primary_type = 'Strawberries',
#'   source_type = 'Spinach', within_cat = 'Farm_type'
#' )
calc_prop_taxa_from_sample_type = function(input, type_header, primary_type,
                                           source_type, within_cat) {
  .calc_prop_main = function(input, type_header, primary_type,
                             source_type) {
    pairs = as.data.frame(t(combn(colnames(
      input$data_loaded
    ), 2)))
    # pairs of a certain type
    pairs$S1_type =
      input$map_loaded[match(pairs[, 1], row.names(input$map_loaded)),
                       type_header]
    pairs$S2_type =
      input$map_loaded[match(pairs[, 2], row.names(input$map_loaded)),
                       type_header]
    types = c(primary_type, source_type)
    conds = pairs$S1_type != pairs$S2_type &
      pairs$S1_type %in% types & pairs$S2_type %in% types
    pairs = pairs[conds,]
    # check if missing sample type(s)
    mis_sts = types[!types %in% c(levels(pairs$S1_type), levels(pairs$S2_type))]
    if (length(mis_sts) > 0)
      stop(paste("Missing sample type: ", mis_sts, " "))
    # function to calc proportion shared for each pair
    .prop_source_in_primary = function(pair, tax_table, primary_type) {
      pair_data = tax_table[, c(pair[1], pair[2])]
      primary_idx = which(c(pair[3], pair[4]) == primary_type)
      no_primary = sum(pair_data[, primary_idx] > 0)
      no_shared = sum(pair_data[, 1] * pair_data[, 2] > 0)
      prop_overlap = no_shared / no_primary
      prop_overlap
    }
    pairs$prop_overlap = apply(pairs, 1, .prop_source_in_primary,
                               input$data_loaded, primary_type)
    mean(pairs$prop_overlap)
  }
  if (missing(within_cat)) {
    .calc_prop_main(input, type_header, primary_type, source_type)
  } else {
    props = c()
    cat_levs = unique(input$map_loaded[, within_cat])
    for (cat_lev in cat_levs) {
      tmp_filt = suppressMessages(filter_data(input, within_cat,
                                              keep_vals = cat_lev))
      n_samples = nrow(tmp_filt$map_loaded)
      if (n_samples < 2)
        stop(paste(cat_lev, 'has less than 2 samples.'))
      tmp_prop = .calc_prop_main(tmp_filt, type_header, primary_type,
                                 source_type)
      props = c(props, tmp_prop)
    }
    data.frame(cat = cat_levs, prop = props)
  }
}

#' @title Return the most abundant taxa in a dataset
#' @description Return the n most abundant taxa as calculated by the mean 
#' sequence counts across all samples.
#' @param input The input dataset as loaded by \code{\link{load_taxa_table}}.
#' @param number_taxa The number of top taxa to display.
#' @concept Misc analyses
return_top_taxa = function(input, number_taxa){
  taxa_ordered = input$taxonomy_loaded[order(rowMeans(input$data_loaded), 
                                             decreasing = T), ]
  head(taxa_ordered, n = number_taxa)
}

#' @title Calculate the changes in taxonomic relative abundances compared to 
#'   controls
#' @description Calculate the change relative to controls (percent) within 
#'   sample groups (blocks) on a per taxon basis.
#' @param ts A taxa summary dataframe as generated by 
#'   \code{\link{summarize_taxonomy}}.
#' @param metadata_map A metadata mapping dataframe. Usually, 
#'   \code{input$map_loaded}.
#' @param block_header The header name from the metadata map that contains the 
#'   block (i.e. sample group) designations.
#' @param treatment_header The header name from the metadata map that contains 
#'   the treatment designations including which samples are controls.
#' @param control_label The designation within the \code{treatment_header} 
#'   column that corresponds to the control samples.
#' @param zero_substitute [OPTIONAL] A value to substitute in for zeros to avoid
#'   infinite percent changes. It can be appropriate to use a lower detection
#'   limit for this value.
#' @concept Misc analyses
#' @examples 
#' # calculate the percent difference in organic from conventional for each 
#' # fruit/veggie
#' ts = summarize_taxonomy(fruits_veggies, level = 2)
#' ts_filt = filter_taxa_from_table(ts, filter_thresh = 0.01)
#' calc_taxa_changes(ts_filt, fruits_veggies$map_loaded, 'Sample_type',
#'                   'Farm_type', 'Conventional')
calc_taxa_changes = function(ts, metadata_map, block_header,
                             treatment_header, control_label, zero_substitute) {
  if(!missing(zero_substitute)) ts[ts == 0] = zero_substitute
  # block control means
  ts2 = dplyr::mutate(ts, Taxon = row.names(ts))
  ts_melt = reshape2::melt(
    ts2, id.vars = 'Taxon',
    value.name = 'Relative_abundance',
    variable.name = 'Sample_ID'
  )
  ts_melt$block = metadata_map[match(ts_melt$Sample_ID, row.names(metadata_map)),
                               block_header]
  ts_melt$Tx = metadata_map[match(ts_melt$Sample_ID, row.names(metadata_map)),
                            treatment_header]
  ts_melt_filt = ts_melt[ts_melt$Tx == control_label, ]
  # ts_melt_filt = dplyr::filter_(ts_melt, lazyeval::interp("Tx == control_label"), control_label = control_label)
  control_means = dplyr::summarise_(dplyr::group_by_(ts_melt_filt, "Taxon", "block"),
                                    mean_RA = "mean(Relative_abundance)")
  # calc percent change
  ts_melt$block_control_mean =
    control_means$mean_RA[match(
      paste0(ts_melt$block, ts_melt$Taxon),
      paste0(control_means$block, control_means$Taxon)
    )]
  ts_melt =
    dplyr::mutate_(
      ts_melt,
      pct_change =
        paste0(
          "(Relative_abundance - block_control_mean) /",
          "block_control_mean * 100"
        ), 
      log2ratio = "log(Relative_abundance / block_control_mean,
                      base = 2)"
    )
  ts_melt
}
