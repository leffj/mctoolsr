# mctoolsr

#######################
# DISSIMILARITY-BASED #
#######################

#' @title Calculate a dissimilarity matrix from a taxa table
#' @description The default is to calculate Bray-Curtis dissimilarities after 
#'  performing square-root transformations on the data. Send requests to add
#'  additional metrics.
#' @param tax_table$data_loaded The observations loaded into the taxa table.
#' @param method The method to use to calculate the dissimilarity metric.
#'  Available methods include: 'bray_sq_trans', 'bray', 'jaccard', '.
#' @return A variable of class 'dist'.
#' @examples 
#' \dontrun{
#' DM = calc_dm(input$data_loaded)
#' }
calc_dm = function(tax_table, method = 'bray_sq_trans'){
  # check for and warn about samples with no sequences
  if(min(colSums(tax_table)) == 0){
    warning('Some samples have no sequences. Samples with low sequence counts
            should be filtered out by rarefying or another acceptable method.')
  }
  if(method == 'bray_sq_trans') {
    # transform otu table (square root transformation)
    tax_table_xform = t(sqrt(tax_table))
    # create dissimilarity matrix from otu table
    vegan::vegdist(tax_table_xform, method = 'bray')
  } else if(method == 'bray') {
    vegan::vegdist(t(tax_table), method = 'bray')
  } else if(method == 'jaccard') {
    vegan::vegdist(t(tax_table), method = 'jaccard')
  } else if(method == 'jaccard') {
    vegan::vegdist(t(tax_table), method = 'morisita')
  } else {
    stop(paste0('Invalid dissimilarity metric specified. See documentation ', 
                'for allowable metrics.'))
  }
}

#' @title Calculate Point Coordinates in an Ordination
#' @description Use before plotting ordination.
#' @param dm Dissimilarity matrix.
#' @param ord_type The type of ordination. 'NMDS' or 'constrained' are the 
#'  current accepted values.
#' @param metadata_map Required if 'constrained' ord_type.
#' @param constrain_factor Required if 'constrained' ord_type.
#' @return A data frame consisting of the coordinates.
#' @examples 
#' \dontrun{
#' ORD = calc_ordination(DM, 'nmds')
#' }
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
#' @param ... Additional arguments passed to ggplot2's \code{geom_point()}.
#' @examples 
#' \dontrun{
#' plot_ordination(input, ORD, "Sample_type", "Farm_type", hull=TRUE)
#' }
plot_ordination = function(input, ordination_axes, color_cat, shape_cat, 
                           hulls = FALSE, ...){
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
    to_plot = data.frame(to_plot, cat2 = input$map_loaded[, shape_cat])
    p = ggplot2::ggplot(to_plot, ggplot2::aes_string(headers[1], headers[2]))
    p = p + ggplot2::geom_point(ggplot2::aes(color = cat, shape = cat2),
                                ..., size = 3, alpha = 0.8)
    p = p + ggplot2::theme_bw()
    p = p + ggplot2::xlab(colnames(to_plot)[1]) + 
      ggplot2::ylab(colnames(to_plot)[2])
    p = p + ggplot2::theme(legend.title = ggplot2::element_blank())
  }
  # plot without shape
  else{
    p = ggplot2::ggplot(to_plot, ggplot2::aes_string(headers[1], headers[2]))
    p = p + ggplot2::geom_point(ggplot2::aes(color=cat),
                                ..., size = 3, alpha = 0.8)
    p = p + ggplot2::theme_bw()
    p = p + ggplot2::xlab(colnames(to_plot)[1]) + 
      ggplot2::ylab(colnames(to_plot)[2])
    p = p + ggplot2::theme(legend.title = ggplot2::element_blank())
  }
  if(hulls){
    p = p + ggplot2::geom_polygon(data = hull_vals, 
                                  ggplot2::aes(fill = cat, color = cat), 
                                  alpha = 0.1)
    p = p + ggplot2::theme(legend.title = ggplot2::element_blank())
  }
  p
}

#' @title Generate an NMDS Plot Quickly
#' @description Used to generate a quick ordination.
#' @param dm Dissimilarity matrix.
#' @param metadata_map The metadata mapping data frame.
#' @param color_cat The metadata map header used to color points.
#' @param shape_cat The metadata map header used for points' shapes (optional).
#' @examples 
#' \dontrun{
#' plot_nmds(DM, metadata_map = input$map_loaded, 'Sample_type', 'Farm_type')
#' }
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

#' @title Generate a dendrogram based on a dissimilarity matrix
#' @param dm Dissimilarity matrix.
#' @param metadata_map The metadata mapping dataframe. 
#'  Typically, input$map_loaded.
#' @param labels The metadata mapping dataframe column name representing the 
#'  intended leaf labels.
#' @param color_by The metadata mapping dataframe column name representing the 
#'  intended leaf label colors.
#' @param method The clustering method to use when creating the dendrogram.
#' @param ... Other parameters passed on to geom_text
#' @examples 
#' \dontrun{
#' plot_dendrogram(DM, input$map_loaded, 'Farm_type', 
#'                 'Sample_type')
#' }
plot_dendrogram = function(dm, metadata_map, labels, color_by, 
                           method = 'complete', ...) {
  if (!requireNamespace("ggdendro", quietly = TRUE)) {
    stop(paste0("'ggdendro' package needed for this function ", 
                "to work. Please install it."), call. = FALSE)
  }
  hc = hclust(dm, method)
  ddata = ggdendro::dendro_data(hc)
  map_rows = match(ddata$labels$label, row.names(metadata_map))
  leaf_labels = as.character(metadata_map[map_rows, labels])
  if(!missing(color_by)) {
    sample_categories = metadata_map[map_rows, color_by]
    ddata$leaf_labels = data.frame(ddata$labels[, 1:2], leaf_labels, 
                                   sample_categories)
  } else {
    ddata$leaf_labels = data.frame(ddata$labels[, 1:2], leaf_labels)
  }
  
  p = ggplot2::ggplot()
  p = p + ggplot2::geom_segment(data = ddata$segments, 
                                ggplot2::aes(x, y, xend = xend, yend = yend))
  if(!missing(color_by)) {
    p = p + ggplot2::geom_text(data = ddata$leaf_labels, 
                               ggplot2::aes(x, y, label = leaf_labels, 
                                            color = sample_categories), 
                               ...,
                               hjust = 0, angle = -90)
  } else {
    p = p + ggplot2::geom_text(data = ddata$leaf_labels, 
                               ggplot2::aes(x, y, label = leaf_labels), 
                               ...,
                               color = 'black', hjust = 0, angle = -90)
  }
  p = p + ggdendro::theme_dendro()
  p = p + ggplot2::scale_y_continuous(lim = c(-0.5, max(ddata$segments$y)*1.05))
  p
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
#' @return The dissimilarity matrix as a 3 column data frame with columns:
#' c('x1', 'x2', 'dist')
#' @examples 
#' \dontrun{
#' DM_cols = convert_dm_to_3_column(DM)
#' }
convert_dm_to_3_column = function(dm){
  if(class(dm) == 'data.frame'){
    dmat = as.dist(dm)
  }
  else{
    dmat = dm
  }
  dmat.clmns = data.frame(t(combn(unlist(labels(dmat)), 2)), as.numeric(dmat))
  names(dmat.clmns) = c('x1', 'x2', 'dist')
  dmat.clmns
}

.add_metadata_to_df = function(x, y, z){
  stop('Deprecated - Please use "add_metadata_to_dm_clmns".')
}

#' @title Add metadata to an additional column in column formatted 
#'  dissimilarities dataframe
#' @param dm_clmns The dissimilarities dataframe produced using \code{
#'  convert_dm_to_3_column()}.
#' @param metadata_map The metadata dataframe.
#' @param cat The header string from the metadata map corresponding to the data
#'  you would like to add to the dissimilarities dataframe.
#' @return A 5 column data frame with the metadata added, the two additional 
#' columnds corresponds to x1 and x2's metadata column equivalent
#' @examples 
#' \dontrun{
#' add_metadata_to_dm_clmns(DM_cols, input$map_loaded, "Sample_type")
#' }
add_metadata_to_dm_clmns = function(dm_clmns, metadata_map, cat){
  cat1 = metadata_map[match(dm_clmns$x1, row.names(metadata_map)), cat]
  cat2 = metadata_map[match(dm_clmns$x2, row.names(metadata_map)), cat]
  dm_clmns_wCat = cbind(dm_clmns, cat1, cat2)
  names(dm_clmns_wCat) = c(names(dm_clmns), paste(cat, "_1", sep=''), 
                           paste(cat, "_2", sep=''))
  dm_clmns_wCat
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
#' @param dm Dissimilarity matrix - typically created using 'calc_dm()'.
#' @param metadata_map The metadata mapping dataframe.
#' @param summarize_by_factor Category in mapping file to summarize by.
#' @param return_map Whether or not to return summarized mapping files. If true,
#'  will return a list (default: FALSE).
#' @return Mean dissimilarities.
#' @examples 
#' \dontrun{
#' calc_mean_dissimilarities(DM, input$map_loaded, 'Sample_Farming')
#' }
calc_mean_dissimilarities = function(dm, metadata_map, summarize_by_factor, 
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
  # check that dm labels match metadata sample IDs
  if(!identical(labels(dm), row.names(metadata_map))) {
    warning('Dissimilarity matrix labels and metadata sample IDs do not match.')
  }
  dm_clmns = convert_dm_to_3_column(dm)
  # list sample 1 and sample 2 factor categories in new clmns
  dm_clmns_wCat = add_metadata_to_dm_clmns(dm_clmns, metadata_map, summarize_by_factor)
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
    mean_map = dplyr::summarise_each(dplyr::group_by_(metadata_map, 
                                                      summarize_by_factor), 
                                     dplyr::funs(.sumry_fun))
    mean_map = as.data.frame(as.matrix(mean_map))
    row.names(mean_map) = mean_map[, summarize_by_factor]
    dm_loaded = as.dist(.convert_one_column_to_matrix(means2))
    map_loaded = mean_map[match(labels(dm_loaded), row.names(mean_map)), ]
    list(dm_loaded = dm_loaded, map_loaded = map_loaded) 
  } else as.dist(.convert_one_column_to_matrix(means2))
}


#' @title Calculate pairwise PERMANOVA results
#' @description The \code{adonis()} function in the \code{vegan} package does
#'  not provide a way to calculate pairwise post-hoc comparisons between
#'  factor levels. This function calculates PERMANOVA results using the 
#'  \code{adonis()} function for all factor level pairs. Raw p values are 
#'  returned, but it is recommended to use the provided FDR corrected p values
#'  since the multiple comparisons can raise your likelihood of false positive
#'  differences.
#' @param dm Dissimilarity matrix of class \code{dist}.
#' @param metadata_map The metadata mapping dataframe with samples matching and
#'  in the same order as the ones in the provided dm.
#' @param compare_header The header in the metadata mapping dataframe with the
#'  factor levels to use for the pairwise comparisons.
#' @return A dataframe with R2 and P values.
#' @examples 
#' \dontrun{
#' calc_pairwise_permanovas(DM, input$map_loaded, 'Sample_Farming')
#' }
calc_pairwise_permanovas = function(dm, metadata_map, compare_header) {
  comp_var = metadata_map[, compare_header]
  comp_pairs = combn(levels(comp_var), 2)
  pval = c()
  R2 = c()
  for(i in 1:ncol(comp_pairs)) {
    pair = comp_pairs[, i]
    dm_w_map = list(dm_loaded = dm, map_loaded = metadata_map)
    dm_w_map$map_loaded$in_pair = comp_var %in% pair
    dm_w_map_filt = filter_dm(dm_w_map, filter_cat = 'in_pair', keep_vals = TRUE)
    m = vegan::adonis(dm_w_map_filt$dm_loaded ~ 
                        dm_w_map_filt$map_loaded[, compare_header])
    pval = c(pval, m$aov.tab$`Pr(>F)`[1])
    R2 = c(R2, m$aov.tab$R2[1])
  }
  results = data.frame(t(comp_pairs), R2, pval)
  results$pvalBon = pval * length(pval)
  results$pvalFDR = round(
    pval * (length(pval) / rank(pval, ties.method = "average")), 
    3)
  results
}