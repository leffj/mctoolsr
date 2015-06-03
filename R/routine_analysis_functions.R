
#############
# FUNCTIONS #
#############

#' @title Load a taxon table for use with mctoolsr
#' @description Load in a taxon table (aka. an OTU table) and a corresponding mapping file 
#' with metadata values. The samples in the loaded taxon table and mapping file 
#' will be in the same order and only samples in both will be loaded. The function
#' can optionally filter samples of a specific type based on the mapping file.
#' This can also be done later via the filter_data() function.
#' 
#' @param tab_fp Taxon table file path
#' @return A list variable with (1) the loaded taxon table, and (2) the loaded mapping file
load_taxon_table = function(tab_fp, map_fp, filter_cat, filter_vals, keep_vals){
  require(tools)
  # load data
  if(file_ext(tab_fp) == 'biom'){
    require(biom)
    data_b = read_biom(tab_fp)
    data = as.data.frame(as.matrix(biom_data(data_b)))
    data_taxonomy = compile_taxonomy(data_b)
  }
  else if(file_ext(tab_fp) == 'txt'){
    data = read.table(tab_fp,sep='\t',skip=1,comment.char='',header=T,check.names=F,row.names=1)
    if(names(data)[ncol(data)] == 'taxonomy'){
      data$taxonomy = NULL
    }
  }
  else stop('Input file must be either biom (.biom) or tab-delimited (.txt) format.')
  map = read.table(map_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  if(class(map) != 'data.frame') warning('Mapping file should have more than one metadata column.')
  # optionally, subset data
    # cant subset if trying to filter out certain values and keep certain values
    # use one or the other
  if(!missing(filter_cat)){
    map.f = .filt_map(map, filter_cat, filter_vals, keep_vals)
  } else map.f = map
  # match up data from dissimilarity matrix with mapping file
  .match_data_components(data, map.f, data_taxonomy)
}


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

.match_data_components = function(tax_table, map, taxonomy){
  samplesToUse = intersect(names(tax_table), row.names(map))
  tax_table.use = tax_table[, match(samplesToUse, names(tax_table))]
  tax_table.use = tax_table.use[rowSums(tax_table.use) != 0, ]
  map.use = map[match(samplesToUse, row.names(map)),]
  if(!missing('taxonomy')) {
    taxonomy.use = taxonomy[match(row.names(tax_table.use), row.names(taxonomy)), ]
    list(data_loaded = tax_table.use, map_loaded = map.use, taxonomy_loaded = taxonomy.use)
  } else {
    list(data_loaded = tax_table.use, map_loaded = map.use)
  }
}


# generates a data frame with all levels of taxanomic info as columns
compile_taxonomy = function(biom_data){
  # get only taxonomy observation metadata from a biom file
  obs_md = observation_metadata(biom_data)
  # replace label for otus with only 1 taxonomy level
  obs_md = sapply(obs_md, function(x) {
    names(x) = gsub('taxonomy$', 'taxonomy1', names(x))
    list(x)})
  # get the taxonomy levels
  otu_full_md = names(which.max(sapply(obs_md, length))) # need to get the otu with all metadata levels
  taxa_levels = names(obs_md[otu_full_md][[1]])[
    grepl('taxonomy', names(obs_md[[1]]), ignore.case = TRUE)]
  # compile taxonomy for each level
  tax_comp = data.frame(row.names = names(obs_md))
  for(l in 1:length(taxa_levels)){
    level_tax_tmp = sapply(obs_md, function(x) x[taxa_levels[l]])
    names(level_tax_tmp) = names(obs_md)
    tax_comp[, taxa_levels[l]] = level_tax_tmp
  }
  tax_comp[is.na(tax_comp)] = 'unclassified'
  tax_comp
}

# level is a single number referring to the taxonomic level
# relative refers to whether output should be sequence counts or relative abundances
# report_higher_tax indicates whether to display all higher taxonomic strings or just 
  # the level of interest
summarize_taxonomy = function(data, level, relative = TRUE, report_higher_tax = TRUE){
  if(report_higher_tax) taxa_strings = apply(data$taxonomy_loaded[1:level], 1, paste0, collapse = '; ')
  else taxa_strings = data$taxonomy_loaded[, level]
  tax_sum = as.data.frame(apply(data$data_loaded, 2, function(x) by(x, taxa_strings, sum)))
  if(relative){
    tax_sum/colSums(data$data_loaded)
  } else tax_sum
}


load_ts_table = function(tab_fp, map_fp, filter_cat, filter_vals, keep_vals){
  require(tools)
  # load data
  if(file_ext(tab_fp) == 'biom'){
    require(biom)
    data = read_biom(tab_fp)
    data = as.data.frame(as.matrix(biom_data(data)))
  }
  else if(file_ext(tab_fp) == 'txt'){
    data = read.table(tab_fp, header=TRUE, sep="\t", row.names=1, comment.char="", check.names=FALSE)
  }
  else stop('Input file must be either biom (.biom) or tab-delimited (.txt) format.')
  map = read.table(map_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  # optionally, subset data
  # cant subset if trying to filter out certain values and keep certain values
  # use one or the other
  if(!missing(filter_cat)){
    map.f = .filt_map(map, filter_cat, filter_vals, keep_vals)
  } else map.f = map
  # match up data from dissimilarity matrix with mapping file
  .match_data_components(data, map.f)
}


load_dm = function(dm_fp, map_fp, filter_cat, filter_vals, keep_vals){
  dm = read.table(dm_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  map = read.table(map_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  # optionally, subset data
  # cant subset if trying to filter out certain values and keep certain values
  # use one or the other
  if(!missing(filter_cat)){
    map.f = .filt_map(map, filter_cat, filter_vals, keep_vals)
  } else map.f = map
  # match up data from dissimilarity matrix with mapping file
  samplesToUse = intersect(names(dm), row.names(map.f))
  dm.use = as.dist(dm[match(samplesToUse,names(dm)), match(samplesToUse,names(dm))])
  map.use = map.f[match(samplesToUse,row.names(map.f)), ]
  # output
  list(dm_loaded = dm.use, map_loaded = map.use)
}

load_2_dms = function(dm1_fp, dm2_fp, map_fp, filter_cat, filter_vals, keep_vals){
  dm1 = read.table(dm1_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  dm2 = read.table(dm2_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  map = read.table(map_fp,sep='\t',comment.char='',header=T,check.names=F,row.names=1)
  # optionally, subset data
  # cant subset if trying to filter out certain values and keep certain values
  # use one or the other
  if(!missing(filter_cat)){
    map.f = .filt_map(map, filter_cat, filter_vals, keep_vals)
  } else map.f = map
  # match up data from dissimilarity matrix with mapping file
  samplesToUse = intersect(intersect(names(dm1), row.names(map.f)), names(dm2))
  dm1.use = as.dist(dm1[match(samplesToUse,names(dm1)), match(samplesToUse,names(dm1))])
  dm2.use = as.dist(dm2[match(samplesToUse,names(dm2)), match(samplesToUse,names(dm2))])
  map.use = map.f[match(samplesToUse,row.names(map.f)), ]
  # output
  list(dm1_loaded = dm1.use, dm2_loaded = dm2.use, map_loaded = map.use)
}

filter_data = function(data, filter_cat, filter_vals, keep_vals){
  # input is list from 'load_data' function
  # cant subset if trying to filter out certain values and keep certain values
  # use one or the other
  if(!missing(filter_cat)){
    map.f = .filt_map(data$map_loaded, filter_cat, filter_vals, keep_vals)
  } else map.f = map
  # match up data from dissimilarity matrix with mapping file
  if('taxonomy_loaded' %in% names(data)){
    .match_data_components(data$data_loaded, map.f, data$taxonomy_loaded)
  } else {
    .match_data_components(data$data_loaded, map.f)
  }
}


filter_taxa = function(table, filter_thresh, taxa_to_keep, taxa_to_remove){
  stop('Deprecated. Please use "filter_taxa_from_table"')
}

filter_taxa_from_table = function(table, filter_thresh, taxa_to_keep, taxa_to_remove){
  # filter taxa from otu table or taxa summary table based on mean abundance
  # optionally, specify additional taxa to keep
  means = apply(table[, 1:ncol(table)], 1, function(x){mean(x,na.rm=TRUE)})
  number_retained = sum((means >= filter_thresh) *1)
  taxa_keep = names(means[means >= filter_thresh])
  if(!missing(taxa_to_keep)) {taxa_keep = taxa_keep[taxa_keep %in% taxa_to_keep]}
  if(!missing(taxa_to_remove)) {taxa_keep = taxa_keep[!taxa_keep %in% taxa_to_remove]}
  table[row.names(table) %in% taxa_keep, ]
}


#' @details Can use one or more of the parameters to do filtering. Threshold 
#'          filtering takes precidence over taxa filtering. If taxa to keep and 
#'          taxa to remove are both included, taxa to remove will be 
#'          removed if the parameter entries conflict.
#' @param input Input data (a list variable) from 'load_data()' functions
#' @param filter_thresh Filter OTUs less than this number based on mean OTU 
#'        table values.
#' @param taxa_to_keep Keep only taxa that contain these names. Vector or string.
#' @param taxa_to_remove Remove taxa that contain these names. Vector or string.
#' @param at_spec_level If included, only keep/remove matches at this specific 
#'        taxonomy level(s) (a number/numbers referring to the taxonomy 
#'        column(s)).
filter_taxa_from_data = function(input, filter_thresh, taxa_to_keep, 
                                 taxa_to_remove, at_spec_level){
  rows_keep = seq(1, nrow(input$data_loaded))
  if(!missing(filter_thresh)){
    means = apply(input$data_loaded[, 1:ncol(input$data_loaded)], 1, 
                  function(x){mean(x, na.rm=TRUE)})
    number_retained = sum((means >= filter_thresh) *1)
    rows_keep = rows_keep[means >= filter_thresh]
  }
  if(missing(at_spec_level)){
    tax_levels = 1:ncol(input$taxonomy_loaded)
  } else tax_levels = at_spec_level
  if(!missing(taxa_to_keep)){
    rows_keep_tmp = sapply(taxa_to_keep, FUN = function(x){
      grep(x, apply(as.data.frame(input$taxonomy_loaded[, tax_levels]), 1, 
                    paste0, collapse = ''))
    })
    rows_keep = intersect(rows_keep, rows_keep_tmp)
    }
  if(!missing(taxa_to_remove)){
    rows_remove = sapply(taxa_to_remove, FUN = function(x){
      grep(x, apply(as.data.frame(input$taxonomy_loaded[, tax_levels]), 1, 
                    paste0, collapse = ''))
      })
    rows_keep = rows_keep[! rows_keep %in% unlist(rows_remove)]
  }
  list(data_loaded = input$data_loaded[rows_keep, ],
       map_loaded = input$map_loaded, 
       taxonomy_loaded = input$taxonomy_loaded[rows_keep, ])
}
  
  
export_otu_table = function(tab, tax_fp, seq_fp, outfp){
  tax = read.table(tax_fp,sep='\t',comment.char='',header=F,check.names=F,row.names=1)
  seqs = read.table(seq_fp,sep='\t',comment.char='',header=F,check.names=F,row.names=1)
  otus = row.names(tab)
  tab.out = data.frame(tab, taxonomy = tax[match(otus, row.names(tax)),1],
                       sequence = seqs[match(otus, row.names(seqs)),1])
  write.table(tab.out, outfp, sep='\t', col.names=NA)
}


single_rarefy = function(data, depth) {
  require(vegan)
  data_filt_samples = data$data_loaded[, colSums(data$data_loaded) >= depth]
  data_rar = as.data.frame(t(rrarefy(t(data_filt_samples), depth)))
  # match up data from dissimilarity matrix with mapping file
  samplesToUse = intersect(names(data_rar), row.names(data$map_loaded))
  data.use = data_rar[, match(samplesToUse, names(data_rar))]
  data.use = data.use[rowSums(data.use) != 0,]
  map.use = data$map_loaded[match(samplesToUse, row.names(data$map_loaded)), ]
  if('taxonomy_loaded' %in% names(data)) {
    taxonomy_loaded.use = data$taxonomy_loaded[match(row.names(data.use), row.names(data$taxonomy_loaded)), ]
    list(data_loaded = data.use, map_loaded = map.use, taxonomy_loaded = taxonomy_loaded.use)
  } else {
    list(data_loaded = data.use, map_loaded = map.use)
  }
}


calc_dm = function(tab){
  require(vegan)
  # check for and warn about samples with no sequences
  if(min(colSums(tab)) == 0){
    warning('Some samples have no sequences. Samples with low sequence counts
            should be filtered out by rarefying or another acceptable method.')
  }
  # transform otu table (square root transformation)
  otuTable.xform = t(sqrt(tab))
  # create dissimilarity matrix from otu table
  otuTable.dist = vegdist(otuTable.xform, method='bray')
  otuTable.dist
}

calc_ordination = function(dm, ord_type, map, constrain_factor){
  require(vegan)
  dm = as.dist(dm)
  if(ord_type == 'NMDS' | ord_type == 'nmds'){
    dm_mds = metaMDS(dm, k=2)
    data.frame(dm_mds$points)
  }
  else if(ord_type == 'constrained'){
    cap = capscale(formula = dm ~ map[, constrain_factor])
    data.frame(scores(cap)$sites)
  }
  else stop('Only NMDS implementd so far.')
  
}

plot_ordination = function(data, ordination_axes, color_cat, shape_cat){
  require(ggplot2)
  if(missing(color_cat)){
    warning('No mapping category to color by.')
    color_vec = rep('none', length(labels(dm)))
  } else color_vec = data$map_loaded[, color_cat]
  to_plot = data.frame(ordination_axes, cat = color_vec)
  headers = colnames(to_plot)
  # plot w shape
  if(!missing(shape_cat)){
    to_plot = data.frame(to_plot, cat2 = data$map_loaded[,shape_cat])
    ggplot(to_plot, aes_string(headers[1], headers[2])) +
      geom_point(size = 3, alpha = 0.8, aes(color = cat, shape = cat2)) + theme_bw() +
#       scale_color_discrete('') + scale_shape_discrete('') +
      xlab(colnames(to_plot)[1]) + ylab(colnames(to_plot)[2])
  }
  # plot without shape
  else{
    ggplot(to_plot, aes_string(headers[1], headers[2])) +
      geom_point(size = 3, alpha = 0.8, aes(color=cat)) + theme_bw() +
#       scale_color_discrete('') + scale_shape_discrete('') +
      xlab(colnames(to_plot)[1]) + ylab(colnames(to_plot)[2])
  }
}

plot_nmds = function(dm, map = NULL, color_cat, shape_cat){
  require(ggplot2)
  if(missing(color_cat)){
    warning('No mapping category to color by.')
    color_vec = rep('none', length(labels(dm)))
  } else color_vec = map[, color_cat]
  # format data and do NMDS
  dm = as.dist(dm)
  dm.mds = metaMDS(dm, k=2)
  # plot w shape
  if(!missing(shape_cat)){
    points = data.frame(dm.mds$points, cat = color_vec, 
                        cat2 = map[, shape_cat])
    ggplot(points, aes(MDS1, MDS2, color = cat, shape = cat2)) +
      geom_point(size = 3, alpha = 0.8) + theme_bw() +
      scale_color_discrete('') + scale_shape_discrete('') 
  }
  # plot without shape
  else{
    points = data.frame(dm.mds$points, cat = color_vec)
    ggplot(points, aes(MDS1, MDS2, color = cat)) +
      geom_point(size = 3, alpha = 0.8) + theme_bw() +
      scale_color_discrete('') + scale_shape_discrete('') 
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

# functions to convert and manipulate dissimilarity matrices in 3 column format

secondary_permanova = function(data, split_cat, cat_level, test_factor){
  data.tmp = filter_data(data, filter_cat = split_cat, keep_vals = cat_level)
  data.tmp = filter_data(data.tmp, filter_cat = test_factor, filter_vals = '#N/A')
  dm.tmp = calc_dm(data.tmp$data_loaded)
  results = adonis(dm.tmp ~ data.tmp$map_loaded[, test_factor])
  results
}

assess_within_category_effects = function(data, category, test_factor){
  for(i in 1:length(unique(data$map_loaded[, category]))){
    level = unique(data$map_loaded[, category])[i]
    print(paste(cat('\n'), level, cat('\n')))
    results = secondary_permanova(data, category, level, test_factor)
    print(results)
  }
}

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

add_metadata_to_df = function(dmat_clmns, map, cat){
  cat1 = map[match(dmat_clmns$x1,row.names(map)),cat]
  cat2 = map[match(dmat_clmns$x2,row.names(map)),cat]
  dmat_clmns_wCat = cbind(dmat_clmns, cat1, cat2)
  names(dmat_clmns_wCat) = c(names(dmat_clmns), paste(cat, "_1", sep=''), paste(cat, "_2", sep=''))
  dmat_clmns_wCat
}

cats_equal = function(x, col1, col2){
  if(x[col1] == x[col2]){"same"} else{"different"}
}

#' Test which order of two paired strings is the recognized order by comparing to a vector of accepted
#' categories
get_combination_category = function(x, accepted_categories){
  if(paste(x, collapse='__') %in% accepted_categories) {return(paste(x, collapse='__'))}
  else if(paste(rev(x), collapse='__') %in% accepted_categories) {return(paste(rev(x), collapse='__'))}
  else {return("Not accepted category")}
}

#' Get the combination category of two vectors containing strings. 
#' This command dereplicates reverse order combinations.
id_treatment_combination = function(col1, col2){
  # get the list of unique categories
  unique_levels = unique(c(as.character(col1), as.character(col2)))
  # get all possible combinations
  combinations = rbind(t(combn(unique_levels, 2)), 
                       t(as.data.frame(lapply(unique_levels, FUN=rep, times=2))))
  combinations = paste(combinations[,1], combinations[,2], sep='__')
  # identify the combination for each pair of categories testing each order
  comparison_types = apply(data.frame(col1, col2), 1, get_combination_category, 
                           accepted_categories = combinations)
  comparison_types
}

#' convert 3 column dissimilarities back to matrix format
convert_one_column_to_matrix = function(df){
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
#' @param dissim_mat Dissimilarity matrix - typically created using 'calc_dm()'
#' @param summarize_by_factor Category in mapping file to summarize by
#' @return Mean dissimilarities
calc_mean_dissimilarities = function(dissim_mat, map, summarize_by_factor){
  require(dplyr)
  dm_clmns = convert_dm_to_3_column(dissim_mat)
  # list sample 1 and sample 2 factor categories in new clmns
  dm_clmns_wCat = add_metadata_to_df(dm_clmns, map, summarize_by_factor)
  # only take samples in mapping file
  dm_clmns_wCat = dm_clmns_wCat[!is.na(dm_clmns_wCat[, 4]) & !is.na(dm_clmns_wCat[, 5]), ]
  # remove rows where distances are comparing samples from the same cat
  dm_clmns_wCat_reduced = dm_clmns_wCat[dm_clmns_wCat[, 4] != dm_clmns_wCat[, 5], ]
  # get pairwise comparison while accounting for differences in category order
  tx_combo = id_treatment_combination(col1 = dm_clmns_wCat_reduced[, 4], col2 = dm_clmns_wCat_reduced[, 5])
  dm_clmns_wCat_reduced = cbind(dm_clmns_wCat_reduced, tx_combo)
  # calc mean dissimilarities
  means = summarize(group_by(dm_clmns_wCat_reduced, tx_combo), mean_dist = mean(dist))
  # convert back to matrix format
  means2 = data.frame(do.call(rbind, strsplit(as.character(means$tx_combo), split = '__')), 
             mean_dist = means$mean_dist)
  convert_one_column_to_matrix(means2)
}


