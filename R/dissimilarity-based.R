# mctoolsr

#######################
# DISSIMILARITY-BASED #
#######################

#' @title Calculate a dissimilarity matrix from a taxa table
#' @description Currently calculates Bray-Curtis dissimilarities after 
#'  performing square-root transformations on the data. Soon, other metrics
#'  will be available.
#' @param tax_table The taxa table.
#' @return A variable of class 'dist'.
calc_dm = function(tax_table){
  # check for and warn about samples with no sequences
  if(min(colSums(tax_table)) == 0){
    warning('Some samples have no sequences. Samples with low sequence counts
            should be filtered out by rarefying or another acceptable method.')
  }
  # transform otu table (square root transformation)
  otuTable.xform = t(sqrt(tax_table))
  # create dissimilarity matrix from otu table
  otuTable.dist = vegan::vegdist(otuTable.xform, method='bray')
  otuTable.dist
}

#' @title Calculate Point Coordinates in an Ordination
#' @description Use before plotting ordination.
#' @param dm Dissimilarity matrix.
#' @param ord_type The type of ordination. 'NMDS' or 'constrained' are the 
#'  current accepted values.
#' @param metadata_map Required if 'constrained' ord_type.
#' @param constrain_factor Required if 'constrained' ord_type.
calc_ordination = function(dm, ord_type, metadata_map, constrain_factor){
  dm = as.dist(dm)
  if(ord_type == 'NMDS' | ord_type == 'nmds'){
    dm_mds = vegan::metaMDS(dm, k=2)
    data.frame(dm_mds$points)
  }
  else if(ord_type == 'constrained'){
    cap = vegan::capscale(formula = dm ~ metadata_map[, constrain_factor])
    data.frame(vegan::scores(cap)$sites)
  }
  else stop('Only NMDS implementd so far.')
  
}

#' @title Generate an Ordination Plot
#' @description Used to generate a plot from the output of \code{
#'  calc_ordination()}
#' @param input The input dataset as loaded by \code{load_taxa_table()}.
#' @param ordination_axes The output of \code{calc_ordination()}.
#' @param color_cat The metadata map header used to color points.
#' @param shape_cat The metadata map header used for points' shapes (optional).
#' @param hulls Whether or not to include an outline around sample categories.
plot_ordination = function(input, ordination_axes, color_cat, shape_cat, 
                           hulls = FALSE){
  if(missing(color_cat)){
    warning('No mapping category to color by.')
    color_vec = rep('none', length(labels(dm)))
  } else color_vec = input$map_loaded[, color_cat]
  to_plot = data.frame(ordination_axes, cat = color_vec)
  names(to_plot)[3] = 'cat'
  headers = colnames(to_plot)
  # hulls prep
  if(hulls){
    .find_hulls = function(df) {df[chull(df), ]}
    hull_vals = dplyr::do(dplyr::group_by(to_plot, cat), .find_hulls(.))
  }
  # plot w/ shape
  if(!missing(shape_cat)){
    to_plot = data.frame(to_plot, cat2 = input$map_loaded[,shape_cat])
    p = ggplot2::ggplot(to_plot, ggplot2::aes_string(headers[1], headers[2]))
    p = p + ggplot2::geom_point(size = 3, alpha = 0.8, 
                                ggplot2::aes(color = cat, shape = cat2))
    p = p + ggplot2::theme_bw()
    p = p + ggplot2::xlab(colnames(to_plot)[1]) + 
      ggplot2::ylab(colnames(to_plot)[2])
  }
  # plot without shape
  else{
    p = ggplot2::ggplot(to_plot, ggplot2::aes_string(headers[1], headers[2]))
    p = p + ggplot2::geom_point(size = 3, alpha = 0.8, ggplot2::aes(color=cat))
    p = p + ggplot2::theme_bw()
    p = p + ggplot2::xlab(colnames(to_plot)[1]) + 
      ggplot2::ylab(colnames(to_plot)[2])
  }
  if(hulls){
    p = p + ggplot2::geom_polygon(data = hull_vals, 
                                  ggplot2::aes(fill = cat, color = cat), 
                                  alpha = 0.1)
  }
  p
}

#' @title Generate an NMDS Plot Quickly
#' @description Used to generate a quick ordination.
#' @param dm Dissimilarity matrix.
#' @param metadata_map The metadata mapping dataframe.
#' @param color_cat The metadata map header used to color points.
#' @param shape_cat The metadata map header used for points' shapes (optional).
plot_nmds = function(dm, metadata_map = NULL, color_cat, shape_cat){
  if(missing(color_cat)){
    warning('No mapping category to color by.')
    color_vec = rep('none', length(labels(dm)))
  } else color_vec = metadata_map[, color_cat]
  # format data and do NMDS
  dm = as.dist(dm)
  dm.mds = vegan::metaMDS(dm, k=2)
  # plot w shape
  if(!missing(shape_cat)){
    points = data.frame(dm.mds$points, cat = color_vec, 
                        cat2 = metadata_map[, shape_cat])
    ggplot2::ggplot(points, 
                    ggplot2::aes(MDS1, MDS2, color = cat, shape = cat2)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) + 
      ggplot2::theme_bw() +
      ggplot2::scale_color_discrete('') + 
      ggplot2::scale_shape_discrete('') 
  }
  # plot without shape
  else{
    points = data.frame(dm.mds$points, cat = color_vec)
    ggplot2::ggplot(points, ggplot2::aes(MDS1, MDS2, color = cat)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) + 
      ggplot2::theme_bw() +
      ggplot2::scale_color_discrete('') + 
      ggplot2::scale_shape_discrete('') 
  }
}

# Interactive plots NOT WORKING
# # @title Generate interactive plot
# # @description Create an interactive plot with the 'clickme' package
# # @param color_by The mapping file header for the factor to color by
# plot_interactive = function(ordination_axes, data, color_by, name_by){
#   to_plot = data.frame(ordination_axes, point_colors = data$map_loaded[, color_by], 
#                        point_names = data$map_loaded[, name_by])
#   clickme('points', to_plot[, 1], to_plot[, 2],
#           color_groups = factor(to_plot$point_colors),
#           names = to_plot$point_names)
# }
# # data = otu_data
# # color_by = 'group'
# # output_dir = outdir


# secondary_permanova = function(data, split_cat, cat_level, test_factor){
#   data.tmp = filter_data(data, filter_cat = split_cat, keep_vals = cat_level)
#   data.tmp = filter_data(data.tmp, filter_cat = test_factor, filter_vals = '#N/A')
#   dm.tmp = calc_dm(data.tmp$data_loaded)
#   results = adonis(dm.tmp ~ data.tmp$map_loaded[, test_factor])
#   results
# }
# 
# assess_within_category_effects = function(data, category, test_factor){
#   for(i in 1:length(unique(data$map_loaded[, category]))){
#     level = unique(data$map_loaded[, category])[i]
#     print(paste(cat('\n'), level, cat('\n')))
#     results = secondary_permanova(data, category, level, test_factor)
#     print(results)
#   }
# }

#' @title Convert dissimilarity matrix to 3 column format
#' @description This is useful for performing analyses on dissimilarity values
#' @param dm Dissimilarity matrix of either class 'dist' or class 'data.frame'
convert_dm_to_3_column = function(dm){
  if(class(dm) == 'data.frame'){
    dmat = as.dist(dm)
  }
  else{
    dmat = dm
  }
  dmat.clmns = data.frame(t(combn(unlist(labels(dmat)),2)),as.numeric(dmat))
  names(dmat.clmns) = c('x1','x2','dist')
  dmat.clmns
}

.add_metadata_to_df = function(x, y, z){
  stop('Deprecated - Please use ".add_metadata_to_dm_clmns".')
}

.add_metadata_to_dm_clmns = function(dmat_clmns, map, cat){
  cat1 = map[match(dmat_clmns$x1, row.names(map)), cat]
  cat2 = map[match(dmat_clmns$x2, row.names(map)), cat]
  dmat_clmns_wCat = cbind(dmat_clmns, cat1, cat2)
  names(dmat_clmns_wCat) = c(names(dmat_clmns), paste(cat, "_1", sep=''), 
                             paste(cat, "_2", sep=''))
  dmat_clmns_wCat
}

# cats_equal = function(x, col1, col2){
#   if(x[col1] == x[col2]){"same"} else{"different"}
# }

#' @keywords internal
# Test which order of two paired strings is the recognized order by comparing to a vector of accepted
# categories
.get_combination_category = function(x, accepted_categories){
  if(paste(x, collapse='__') %in% accepted_categories) {return(paste(x, collapse='__'))}
  else if(paste(rev(x), collapse='__') %in% accepted_categories) {return(paste(rev(x), collapse='__'))}
  else {return("Not accepted category")}
}

#' @keywords internal
# Get the combination category of two vectors containing strings. 
# This command dereplicates reverse order combinations.
.id_treatment_combination = function(col1, col2){
  # get the list of unique categories
  unique_levels = unique(c(as.character(col1), as.character(col2)))
  # get all possible combinations
  combinations = rbind(t(combn(unique_levels, 2)), 
                       t(as.data.frame(lapply(unique_levels, FUN=rep, times=2))))
  combinations = paste(combinations[,1], combinations[,2], sep='__')
  # identify the combination for each pair of categories testing each order
  comparison_types = apply(data.frame(col1, col2), 1, .get_combination_category, 
                           accepted_categories = combinations)
  comparison_types
}

#' @keywords internal
# convert 3 column dissimilarities back to matrix format
.convert_one_column_to_matrix = function(df){
  # initialize matrix with dimensions equal to number of unique categories
  uNames = unique(c(as.character(df[, 1]), as.character(df[, 2])))
  mean_dists_mat = data.frame(matrix(ncol = length(uNames), nrow = length(uNames)))
  names(mean_dists_mat) = uNames
  row.names(mean_dists_mat) = uNames
  for(i in 1:nrow(df)){
    mean_dists_mat[as.character(df[i, 1]), as.character(df[i, 2])] = as.character(df[i, 3])
  }
  for(i in 1:nrow(mean_dists_mat)){
    for(k in 1:ncol(mean_dists_mat)){
      if(is.na(mean_dists_mat[i,k])){
        mean_dists_mat[i,k] = mean_dists_mat[k,i]
      }
    }
  }
  mean_dists_mat[is.na(mean_dists_mat)] = 0
  mean_dists_mat
}

#' @title Calculate mean dissimilarities using a metadata factor
#' @description Calculate mean dissimilarities across all levels of a given factor
#' 
#' @param dissim_mat Dissimilarity matrix - typically created using 'calc_dm()'.
#' @param metadata_map The metadata mapping dataframe.
#' @param summarize_by_factor Category in mapping file to summarize by.
#' @param return_map Whether or not to return summarized mapping files. If true,
#'  will return a list (default: FALSE).
#' @return Mean dissimilarities.
calc_mean_dissimilarities = function(dissim_mat, metadata_map, summarize_by_factor, 
                                     return_map = FALSE){
  .sumry_fun = function(x){
    if(is.numeric(x)){
      mean(x)
    } else {
      if(length(unique(x)) == 1){
        unique(x)
      } else NA
    }
  }
  dm_clmns = convert_dm_to_3_column(dissim_mat)
  # list sample 1 and sample 2 factor categories in new clmns
  dm_clmns_wCat = .add_metadata_to_dm_clmns(dm_clmns, metadata_map, summarize_by_factor)
  # only take samples in mapping file
  dm_clmns_wCat = dm_clmns_wCat[!is.na(dm_clmns_wCat[, 4]) & !is.na(dm_clmns_wCat[, 5]), ]
  # remove rows where distances are comparing samples from the same cat
  dm_clmns_wCat_reduced = dm_clmns_wCat[dm_clmns_wCat[, 4] != dm_clmns_wCat[, 5], ]
  # get pairwise comparison while accounting for differences in category order
  tx_combo = .id_treatment_combination(col1 = dm_clmns_wCat_reduced[, 4], 
                                       col2 = dm_clmns_wCat_reduced[, 5])
  dm_clmns_wCat_reduced = cbind(dm_clmns_wCat_reduced, tx_combo)
  # calc mean dissimilarities
  means = dplyr::summarize(dplyr::group_by(dm_clmns_wCat_reduced, tx_combo), 
                           mean_dist = mean(dist))
  # convert back to matrix format
  means2 = data.frame(do.call(rbind, strsplit(as.character(means$tx_combo), 
                                              split = '__')), 
                      mean_dist = means$mean_dist)
  if(return_map){
    mean_map = dplyr::summarise_each(dplyr::group_by_(metadata_map, summarize_by_factor), 
                                     dplyr::funs(.sumry_fun))
    list(dm = as.dist(.convert_one_column_to_matrix(means2)), map_loaded = mean_map)
  } else as.dist(.convert_one_column_to_matrix(means2))
}