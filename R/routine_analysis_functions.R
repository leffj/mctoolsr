# source('~/Software/mctoolsr/R/differences_in_taxa_functions.R')
# devtools::use_package('tools')
# devtools::use_package('biom')
# devtools::use_package('reshape2')
# devtools::use_package('dplyr')
# devtools::use_package('ggplot2')
# devtools::use_package('vegan')
# devtools::use_package('nlme')
# devtools::use_package('VennDiagram', 'Suggests')

.onAttach = function(libname, pkgname) {
  packageStartupMessage("You're loading mctoolsr. Direct inquiries to:
                        'https://github.com/leffj/mctoolsr'")
}

.onLoad = function(libname, pkgname) {
  op = options()
  op.devtools = list(
    mctoolsr.path = "",
    mctoolsr.install.args = "",
    mctoolsr.name = "mctoolsr",
    mctoolsr.desc.author = '"Jon Leff <jonathan.leff@colorado.edu> [aut, cre]"',
    mctoolsr.desc.license = "GPL-3",
    mctoolsr.desc.suggests = NULL,
    mctoolsr.desc = list()
  )
  toset = !(names(op.devtools) %in% names(op))
  if(any(toset)) options(op.devtools[toset])
  
  invisible()
}

#############
# FUNCTIONS #
#############

#' @title Load a taxon table for use with mctoolsr
#' @description Load in a taxon table (aka. an OTU table) and a corresponding 
#'  mapping file with metadata values. The samples in the loaded taxon table 
#'  and mapping file will be in the same order and only samples in both will 
#'  be loaded. The function can optionally filter samples of a specific type 
#'  based on the mapping file. This can also be done later via the filter_data() 
#'  function.
#' @param tab_fp Taxon table filepath.
#' @param map_fp Metadata mapping filepath.
#' @param filter_cat The map_fp header string for the factor you would like 
#'  to use to filter samples.
#' @param filter_vals The values within the filter category (vector or single 
#'  value) you would like to use to remove samples from the imported data.
#' @param keep_vals Alternatively, keep only samples represented by these 
#'  values.
#' @return A list variable with (1) the loaded taxon table, and (2) the loaded 
#'  mapping file.
#' @examples 
#' \dontrun{
#' load_taxon_table("filepath_to_OTU_table.txt", "filepath_to_mapping_file.txt",
#'   "sample_type", filter_vals = "blank")
#' }
load_taxon_table = function(tab_fp, map_fp, filter_cat, filter_vals, keep_vals){
  # load data
  if(tools::file_ext(tab_fp) == 'biom'){
    data_b = biom::read_biom(tab_fp)
    data = as.data.frame(as.matrix(biom::biom_data(data_b)))
    data_taxonomy = .compile_taxonomy(data_b)
  }
  else if(file_ext(tab_fp) == 'txt'){
    if(readChar(tab_fp, nchars = 4) == "#OTU"){
      data = read.table(tab_fp, sep='\t', comment.char='', header=T, 
                        check.names=F, row.names=1)
    } else {
      data = read.table(tab_fp, sep='\t', skip=1, comment.char='', header=T, 
                        check.names=F, row.names=1)
    }
    if(names(data)[ncol(data)] == 'taxonomy'){
      data_taxonomy = .parse_taxonomy(data$taxonomy)
      row.names(data_taxonomy) = row.names(data)
      data$taxonomy = NULL
    }
  }
  else stop('Input file must be either biom (.biom) or tab-delimited (.txt) format.')
  map = read.table(map_fp, sep = '\t', comment.char = '', header = T, 
                   check.names = F, row.names = 1)
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

.match_data_components = function(tax_table, map, taxonomy){
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

# generates a data frame with all levels of taxanomic info as columns
.compile_taxonomy = function(biom_dat){
  # get only taxonomy observation metadata from a biom file
  obs_md = biom::observation_metadata(biom_dat)
  if(class(obs_md) == 'list'){
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
  } else if(class(obs_md) == 'data.frame'){
    obs_md
  } else stop('Error compiling taxonomy.')
}

#' @description Convert a taxonomy character vector pulled from a text OTU 
#' table to a data frame
.parse_taxonomy = function(taxonomy_vec){
  tmp = sapply(taxonomy_vec, function(x) strsplit(as.character(x), split = '; *')[[1]])
  # if not all taxonomic levels present for each OTU, this is necessary to 
  # convert to data frame
  if(class(tmp) == 'list'){
    n_obs = sapply(tmp, length)
    seq_max = seq_len(max(n_obs))
    tmp = (sapply(tmp, "[", i = seq_max))
  }
  tmp = as.data.frame(t(tmp))
  potential_names = c('taxonomy1', 'taxonomy2', 'taxonomy3', 'taxonomy4', 
                      'taxonomy5', 'taxonomy6', 'taxonomy7', 'taxonomy8', 
                      'taxonomy9', 'taxonomy10', 'taxonomy11', 'taxonomy12')
  names(tmp) = potential_names[1:ncol(tmp)]
  tmp
}

#' @title convert_to_relative_abundances
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

# Means if relative abundance
# level is a single number referring to the taxonomic level
# relative refers to whether output should be sequence counts or relative abundances
# report_higher_tax indicates whether to display all higher taxonomic strings or just 
  # the level of interest
summarize_taxonomy = function(data, level, relative = TRUE, 
                              report_higher_tax = TRUE){
  if(report_higher_tax) taxa_strings = apply(data$taxonomy_loaded[1:level], 1, 
                                             paste0, collapse = '; ')
  else taxa_strings = data$taxonomy_loaded[, level]
  no_taxa = length(unique(taxa_strings))
  tax_sum = apply(data$data_loaded, 2, function(x) by(x, taxa_strings, sum))
  if(no_taxa == 1) tax_sum = data.frame(t(tax_sum), 
                                        row.names = unique(taxa_strings))
  else tax_sum = as.data.frame(tax_sum)
  if(relative){
    convert_to_relative_abundances(tax_sum)
  } else tax_sum
}

load_ts_table = function(tab_fp, map_fp, filter_cat, filter_vals, keep_vals){
  # load data
  if(file_ext(tab_fp) == 'biom'){
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

plot_taxa_bars = function(taxa_summary_df, metadata_map, factor, num_taxa){
  taxa_summary_df$taxon = row.names(taxa_summary_df)
  taxa_summary_df_melted = melt(taxa_summary_df, variable.name = 'Sample_ID', 
                                id.vars = 'taxon')
  group_by_levels = metadata_map[match(taxa_summary_df_melted$Sample_ID, 
                                       row.names(metadata_map)), factor]
  taxa_summary_df_melted$group_by = group_by_levels
  mean_tax_vals = summarise(group_by(taxa_summary_df_melted, group_by, taxon), 
                            mean_value = mean(value))
  # get top taxa and convert other to 'other'
  mean_tax_vals_sorted = mean_tax_vals[order(mean_tax_vals$mean_value, 
                                             decreasing = TRUE), ]
  top_taxa = unique(mean_tax_vals_sorted$taxon)[1:num_taxa]
  mean_tax_vals_sorted$taxon[!mean_tax_vals_sorted$taxon %in% top_taxa] = 'Other'
  to_plot = summarise(group_by(mean_tax_vals_sorted, group_by, taxon), 
                      mean_value = sum(mean_value))
  # plot
  ggplot(to_plot, aes(group_by, mean_value, fill = taxon)) +
    geom_bar(stat = 'identity') + ylab('') + xlab('') +
    theme(legend.title=element_blank())
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

#' @description Function to match up sample order from two datasets that contain
#' some overlapping sample IDs. Sample IDs that are not present in both
#' datasets will be dropped. The output is a list containing the two filtered
#' datasets in the same order as they were input.
match_datasets = function(ds1, ds2){
  common_samples = intersect(names(ds1$data_loaded), names(ds2$data_loaded))
  ds1$map_loaded$common_sample = row.names(ds1$map_loaded) %in% common_samples
  ds1_filt = filter_data(ds1, 'common_sample', FALSE)
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
    if(length(rows_keep_tmp[[1]]) == 0){
      stop('Taxon not found.')
    }
    rows_keep = intersect(rows_keep, rows_keep_tmp)
    }
  if(!missing(taxa_to_remove)){
    rows_remove = sapply(taxa_to_remove, FUN = function(x){
      grep(x, apply(as.data.frame(input$taxonomy_loaded[, tax_levels]), 1, 
                    paste0, collapse = ''))
      })
    if(length(rows_remove[[1]]) == 0){
      stop('Taxon not found.')
    }
    rows_keep = rows_keep[! rows_keep %in% unlist(rows_remove)]
  }
  list(data_loaded = input$data_loaded[rows_keep, ],
       map_loaded = input$map_loaded, 
       taxonomy_loaded = droplevels(input$taxonomy_loaded[rows_keep, ]))
}

# export_otu_table = function(tab, tax_fp, seq_fp, outfp){
#   tax = read.table(tax_fp,sep='\t',comment.char='',header=F,check.names=F,row.names=1)
#   seqs = read.table(seq_fp,sep='\t',comment.char='',header=F,check.names=F,row.names=1)
#   otus = row.names(tab)
#   tab.out = data.frame(tab, taxonomy = tax[match(otus, row.names(tax)),1],
#                        sequence = seqs[match(otus, row.names(seqs)),1])
#   write.table(tab.out, outfp, sep='\t', col.names=NA)
# }
export_otu_table = function(input, out_fp){
  table = input$data_loaded
  taxonomy = apply(input$taxonomy_loaded, 1, paste, collapse = '; ')
  out_tab = data.frame(OTU_ID = row.names(table), table, taxonomy)
  names(out_tab)[1] = '#OTU ID'
  write('#Exported from mctoolsr', out_fp)
  suppressWarnings(write.table(out_tab, out_fp, sep = '\t', row.names = FALSE, append = TRUE))
}

single_rarefy = function(data, depth) {
  data_filt_samples = data$data_loaded[, colSums(data$data_loaded) >= depth]
  data_rar = as.data.frame(t(vegan::rrarefy(t(data_filt_samples), depth)))
  .match_data_components(data_rar, data$map_loaded, data$taxonomy_loaded)
}

calc_dm = function(tab){
  # check for and warn about samples with no sequences
  if(min(colSums(tab)) == 0){
    warning('Some samples have no sequences. Samples with low sequence counts
            should be filtered out by rarefying or another acceptable method.')
  }
  # transform otu table (square root transformation)
  otuTable.xform = t(sqrt(tab))
  # create dissimilarity matrix from otu table
  otuTable.dist = vegan::vegdist(otuTable.xform, method='bray')
  otuTable.dist
}

calc_ordination = function(dm, ord_type, map, constrain_factor){
  dm = as.dist(dm)
  if(ord_type == 'NMDS' | ord_type == 'nmds'){
    dm_mds = vegan::metaMDS(dm, k=2)
    data.frame(dm_mds$points)
  }
  else if(ord_type == 'constrained'){
    cap = vegan::capscale(formula = dm ~ map[, constrain_factor])
    data.frame(vegan::scores(cap)$sites)
  }
  else stop('Only NMDS implementd so far.')
  
}

plot_ordination = function(data, ordination_axes, color_cat, shape_cat, 
                           hulls = FALSE){
  if(missing(color_cat)){
    warning('No mapping category to color by.')
    color_vec = rep('none', length(labels(dm)))
  } else color_vec = data$map_loaded[, color_cat]
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
    to_plot = data.frame(to_plot, cat2 = data$map_loaded[,shape_cat])
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

plot_nmds = function(dm, map = NULL, color_cat, shape_cat){
  if(missing(color_cat)){
    warning('No mapping category to color by.')
    color_vec = rep('none', length(labels(dm)))
  } else color_vec = map[, color_cat]
  # format data and do NMDS
  dm = as.dist(dm)
  dm.mds = vegan::metaMDS(dm, k=2)
  # plot w shape
  if(!missing(shape_cat)){
    points = data.frame(dm.mds$points, cat = color_vec, 
                        cat2 = map[, shape_cat])
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

.add_metadata_to_df = function(x, y, z){
  stop('Deprecated - Please use "add_metadata_to_dm_clmns".')
}

add_metadata_to_dm_clmns = function(dmat_clmns, map, cat){
  cat1 = map[match(dmat_clmns$x1, row.names(map)), cat]
  cat2 = map[match(dmat_clmns$x2, row.names(map)), cat]
  dmat_clmns_wCat = cbind(dmat_clmns, cat1, cat2)
  names(dmat_clmns_wCat) = c(names(dmat_clmns), paste(cat, "_1", sep=''), 
                             paste(cat, "_2", sep=''))
  dmat_clmns_wCat
}

cats_equal = function(x, col1, col2){
  if(x[col1] == x[col2]){"same"} else{"different"}
}

#' Test which order of two paired strings is the recognized order by comparing to a vector of accepted
#' categories
.get_combination_category = function(x, accepted_categories){
  if(paste(x, collapse='__') %in% accepted_categories) {return(paste(x, collapse='__'))}
  else if(paste(rev(x), collapse='__') %in% accepted_categories) {return(paste(rev(x), collapse='__'))}
  else {return("Not accepted category")}
}

#' Get the combination category of two vectors containing strings. 
#' This command dereplicates reverse order combinations.
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

#' convert 3 column dissimilarities back to matrix format
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
#' @param dissim_mat Dissimilarity matrix - typically created using 'calc_dm()'
#' @param summarize_by_factor Category in mapping file to summarize by
#' @param return_map Whether or not to return summarized mapping files. If true,
#'  will return a list (default: FALSE)
#' @return Mean dissimilarities
calc_mean_dissimilarities = function(dissim_mat, map, summarize_by_factor, 
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
  dm_clmns_wCat = add_metadata_to_dm_clmns(dm_clmns, map, summarize_by_factor)
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
    mean_map = dplyr::summarise_each(dplyr::group_by_(map, summarize_by_factor), 
                                   dplyr::funs(.sumry_fun))
    list(dm = as.dist(.convert_one_column_to_matrix(means2)), map_loaded = mean_map)
  } else as.dist(.convert_one_column_to_matrix(means2))
}


#' @details Function to show contributions of specific taxa to variation among 
#'          communities using Mann-Whitney (2 factor levels), Kruskal-Wallis 
#'          (more than 2) tests, or more complex models
#' @param taxa_smry_df Taxa summary data frame
#' @param metadata_map Mapping file
#' @param out_fp Test results output filepath (OPTIONAL)
#' @param factor Mapping file header (in quotation marks) of factor for which 
#'        you are testing for differences
#' @param filter_level The minimum mean value needed in at least one 
#'        of the factor levels for a taxon to be retained in the analysis
#' @param test_type either 'MW', 'KW', or 'custom' (i.e. Wilcoxon/Mann-Whitney U 
#'        for 2 factor levels or Kruskal-Wallis for more than two factor levels).
#'        See details for custom test/model implementation.
#' @param custom_test_function Name of custom test function
#' @param smry_fun The function to summarize values by (Default: mean)
taxa_summary_by_sample_type = function(taxa_smry_df, metadata_map, factor, 
                                       filter_level, test_type, grouping_factor, 
                                       custom_test_function, smry_fun = mean, 
                                       out_fp){
  if(!missing(filter_level)){
    # filter taxa summary table by abundance in any/either factor level
    taxa_smry_df = .filter_taxa_dit(taxa_smry_df, metadata_map, filter_level, 
                                    factor, smry_fun = smry_fun)
  }
  if(!missing(grouping_factor)){
    test_results = .run_test(taxa_smry_df, metadata_map, factor, test_type, 
                             grouping_factor, smry_fun = smry_fun)
  }
  else if(!missing(custom_test_function)){
    test_results = .run_test(taxa_smry_df, metadata_map, factor, test_type,
                             cust_test = custom_test_function, 
                             smry_fun = smry_fun)
  } 
  else{
    test_results = .run_test(taxa_smry_df, metadata_map, factor, test_type,
                             smry_fun = smry_fun)
  }
  # Sort by pvalues 
  test_results = test_results[with(test_results, order(pvals)), ]
  # output data
  if(!missing(out_fp)){
    write.table(test_results, file = out_fp, sep = ",", row.names = TRUE, 
                col.names = NA)  
  } else test_results
}

# Plot venn diagram showing OTU overlap between multiple sample categories
#' @param pres_thresh The threshold relative abundance for an OTU to be 
#' considered present for a given sample category
## testing
# infp = '/Users/leffj/Software/mctoolsr/testing/otu_tab1.biom'
# mapfp = '/Users/leffj/Software/mctoolsr/testing/map1.txt'
plot_venn_diagram = function(input_data, category, pres_thresh){
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    stop("'VennDiagram' package needed for this function to work. Please 
         install it.", call. = FALSE)
  }
  .venn_3cat = function(data, thresh){
    c1count = sum(data[, 2] >= thresh)
    c2count = sum(data[, 3] >= thresh)
    c3count = sum(data[, 4] >= thresh)
    n12ct = sum(rowSums(data[, c(2, 3)] > thresh) == 2)
    n13ct = sum(rowSums(data[, c(2, 4)] > thresh) == 2)
    n23ct = sum(rowSums(data[, c(3, 4)] > thresh) == 2)
    n123ct = sum(rowSums(data[, c(2, 3, 4)] > thresh) == 3)
    grid.newpage()
    draw.triple.venn(c1count, c2count, c3count, n12ct, n23ct, n13ct, n123ct,
                     category = names(data)[2:4],
                     fill = c('blue', 'red', 'orange'))
  }
  
  .venn_4cat = function(data, thresh){
    c1count = sum(data[, 2] > thresh)
    c2count = sum(data[, 3] > thresh)
    c3count = sum(data[, 4] > thresh)
    c4count = sum(data[, 5] > thresh)
    n12ct = sum(rowSums(data[, c(2, 3)] > thresh) == 2)
    n13ct = sum(rowSums(data[, c(2, 4)] > thresh) == 2)
    n14ct = sum(rowSums(data[, c(2, 5)] > thresh) == 2)
    n23ct = sum(rowSums(data[, c(3, 4)] > thresh) == 2)
    n24ct = sum(rowSums(data[, c(3, 5)] > thresh) == 2)
    n34ct = sum(rowSums(data[, c(4, 5)] > thresh) == 2)
    n123ct = sum(rowSums(data[, c(2, 3, 4)] > thresh) == 3)
    n124ct = sum(rowSums(data[, c(2, 3, 5)] > thresh) == 3)
    n134ct = sum(rowSums(data[, c(2, 4, 5)] > thresh) == 3)
    n234ct = sum(rowSums(data[, c(3, 4, 5)] > thresh) == 3)
    n1234ct = sum(rowSums(data[, c(2, 3, 4, 5)] > thresh) == 4)
    grid.newpage()
    draw.quad.venn(c1count, c2count, c3count, c4count, n12ct, n13ct, n14ct, 
                   n23ct, n24ct, n34ct, n123ct, n124ct, n134ct, n234ct, n1234ct,
                   category = names(data)[2:5],
                   fill = c('blue', 'red', 'green', 'orange'))
  }
  
  # get relative abundances
  otu_RAs = decostand(input_data$data_loaded, method='total', MARGIN=2)
  # get means of otus by metadata category to check for presence
  if(! category %in% names(input_data$map_loaded)){
    stop('Category header not found in mapping file.')
  }
  otu_RAs_t = as.data.frame(t(otu_RAs))
  otu_RAs_t$cat = input_data$map_loaded[, category]
  otu_RAs_melted = melt(otu_RAs_t, id.vars = 'cat')
  otu_RAs_means = dplyr::summarize(dplyr::group_by(otu_RAs_melted, variable, 
                                                   cat), 
                            mean_RA = mean(value))
  # plot diagram, discount OTUs with relative abundances lower than threshold
  otu_RAs_means_cast = reshape::dcast(otu_RAs_means, variable ~ cat, value.var = 'mean_RA')
  if(ncol(otu_RAs_means_cast) - 1 == 3){
    .venn_3cat(otu_RAs_means_cast, pres_thresh)
  } else if (ncol(otu_RAs_means_cast) - 1 == 4){
    .venn_4cat(otu_RAs_means_cast, pres_thresh)
  } else {
    stop('Can only plot Venn with 3 or 4 category levels.')
  }
}

#' @description Calculate diversity values from a taxon table
calc_diversity = function(taxon_table, metric){
  metrics = c('richness', 'shannon', 'simpson')
  if(metric == 'richness'){
    apply(taxon_table, 2, function(x) length(x[x > 0]))
  }
  else if(metric %in% c('shannon', 'simpson')){
    diversity(taxon_table, index = metric, MARGIN = 2)
  }
  else {
    stop(paste0('Metric, "', metric, '", is not a valid metric.\n  Valid metrics are: ',
                paste(metrics, collapse = ', '), '.'))
  }
}

#' @description Plot diversity using boxplots across levels of a given factor
#' @param variable Variable in mapping file to plot by
plot_diversity = function(input, variable, metric){
  diversity = calc_diversity(input$data_loaded, metric)
  to_plot = data_frame(diversity, input$map_loaded[, variable])
  names(to_plot) = c(metric, variable)
  ggplot(to_plot, aes_string(variable, metric, fill = variable)) +
    geom_boxplot() + theme_bw() + xlab('')
}



