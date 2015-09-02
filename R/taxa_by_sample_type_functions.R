# mctoolsr

#######################
# TAXA BY SAMPLE TYPE #
#######################

# This code will: (1) Filter the taxa summary to remove taxa that do not meet
# an abundance threshold in any factor level. This is based on mean abundance.
# (2) Calculate which taxa have differences in relative abundance among factor
# levels. This is based on either Mann-Whitney tests or Kruskal-Wallis (both
# non-parametric tests), although custom models can be used. Use Mann-Whitney 
# for 2 factor levels and K-W for more than 2. (3) Output results including 
# adjusted (Bonferroni and FDR) p-values and means.

#' @keywords internal
# get metadata values for a specific variable in the same order as the samples
# in the taxa table
.get_metadata = function(t_table, map_file, variable){
  map_file[match(names(t_table), row.names(map_file)), variable]
}

#' @keywords internal
# function to filter taxa
.filter_taxa_dit = function(t_table, map_file, f_level, f_factor, smry_fun){
  # Check if the t_table only has one sample
  if(class(t_table) == "numeric"){
    "skip"
  } else {
    factorMeta = .get_metadata(t_table, map_file, f_factor)
    rowsToKeep = c()
    for(i in 1:nrow(t_table)){
      #   in the row, calculate means for each factor level and keep if one is 
      #   greater than filter
      meanAbunds <- NULL
      meanAbunds <- aggregate(as.numeric(t(t_table[i, ])), list(factorMeta), 
                              smry_fun)
      if(max(meanAbunds$x) >= f_level){
        rowsToKeep <- c(rowsToKeep, i)
      }
    }
    t_table[rowsToKeep, ]
  }
}

#' @keywords internal
# Run Wilcoxon Rank-Sum test (Mann-Whitney U test) and return p-value
.run_MW_test = function(dependent, factor){
  # check for only two factor levels
  if(length(unique(factor)) != 2) print('Mann-Whitney test requires exacly two 
                                      factor levels.')
  wilcox.test(formula = dependent ~ factor)$p.value
}

#' @keywords internal
# Run Kruskal-Wallis test
.run_KW_test = function(dependent, factor){
  kruskal.test(formula = dependent ~ factor)$p.value
}

#' @keywords internal
# Run 2-way NP test
# uses a rank transformation then lme model
.run_2WNP_test = function(dependent, factor, g_factor){
  dep.ranked = rank(dependent, ties.method = 'average')
  model = nlme::lme(fixed = dep.ranked ~ factor, random =~ 1|g_factor)
  anova(model)[["p-value"]][2]
}

#' @keywords internal
# Run custom test
# custom function should take a vector of values for an individual taxon
# and return a p-value corresponding to the test performed
.run_custom_test = function(dependent, cust_func_name){
  cust_func_name(dependent)
}

#' @keywords internal
# run statistical test (Mann-whitney, Kruskal-Wallis, or 2-way NP) on each taxon
# in a provided taxa table
.run_test = function(t_table, map_file, fctr, type, g_fctr, cust_test, smry_fun){
  fctrMeta = as.factor(as.vector(.get_metadata(t_table, map_file, fctr)))
  if(!missing(g_fctr)) gfctrMeta = as.factor(as.vector(.get_metadata(t_table, 
                                                                     map_file, 
                                                                     g_fctr)))
  pvals = c()
  for(i in 1:nrow(t_table)){
    if(type == 'MW') pvals = c(pvals, .run_MW_test(as.vector(t(t_table[i, ])), 
                                                   fctrMeta))
    else if(type == 'KW') pvals = c(pvals, .run_KW_test(as.vector(t(
      t_table[i, ])), fctrMeta))
    else if(type == '2WNP') pvals = c(pvals, .run_2WNP_test(as.vector(t(
      t_table[i, ])), fctrMeta, gfctrMeta))
    else if(type == 'custom') pvals = c(pvals, .run_custom_test(as.vector(t(
      t_table[i, ])), cust_test))
    else print('Invalid test type specified')
    if(i == 1){
      meanAbunds = aggregate(as.numeric(t(t_table[i, ])), list(fctrMeta), smry_fun)
    } else{
      means = aggregate(as.numeric(t(t_table[i, ])), list(fctrMeta), smry_fun)[, 2]
      meanAbunds = cbind(meanAbunds, means)
    }
  }
  # generate bonforroni corrected pvals
  pvalsBon = pvals * length(pvals)
  # generate FDR corrected pvals (taken from otu_category_significance.py)
  # Ranks p-values low to high and multiplies each p-value by the number of
  # comparisons divided by the rank.
  pvalsFDR = pvals * (length(pvals) / rank(pvals, ties.method="average"))
  # prep means to be added
  factorLevels = as.character(meanAbunds[, 1])
  meanAbunds[, 1] = NULL
  # make result df
  result = as.data.frame(cbind(pvals, pvalsBon, pvalsFDR, t(meanAbunds)))
  row.names(result) = row.names(t_table)
  colnames(result) = c("pvals", "pvalsBon", "pvalsFDR", factorLevels)
  result
}

#' @title Further summarize output from summarize_taxonomy by sample type
#' @details Function to show contributions of specific taxa to variation among 
#'  communities using Mann-Whitney (2 factor levels), Kruskal-Wallis (more 
#'  than 2) tests, or more complex models.
#' @param taxa_smry_df Taxa summary data frame.
#' @param metadata_map Mapping file.
#' @param out_fp (Optional) Test results output filepath. Written as a csv file.
#' @param type_header Mapping file header (in quotation marks) of factor for 
#'  which you are testing for differences.
#' @param filter_level (Optional) The minimum mean value needed in at least one.
#'  of the factor levels for a taxon to be retained in the analysis.
#' @param test_type (Optional). If omitted, no test is performed and no p-values
#'  are reported. Otherwise, either 'MW', 'KW', or 'custom' (i.e. 
#'  Wilcoxon/Mann-Whitney U for 2 factor levels or Kruskal-Wallis for more than 
#'  two factor levels). See details for custom test/model implementation.
#' @param grouping_factor (Optional) Used with 2-way tests.
#' @param custom_test_function (Optional) Name of custom test function.
#' @param smry_fun (Optional) The function to summarize values by (Default: mean).
taxa_summary_by_sample_type = function(taxa_smry_df, metadata_map, type_header, 
                                       filter_level, test_type, grouping_factor, 
                                       custom_test_function, smry_fun = mean, 
                                       out_fp){
  if(!missing(filter_level)){
    # filter taxa summary table by abundance in any/either factor level
    taxa_smry_df = .filter_taxa_dit(taxa_smry_df, metadata_map, filter_level, 
                                    type_header, smry_fun = smry_fun)
  }
  # if no tests, just summarize by factor
  if(missing(test_type)){
    results = t(apply(taxa_smry_df, 1, function(x) {
      tapply(as.numeric(x), metadata_map[, type_header], mean)
    }))
  }
  # else run test
  else{
    if(!missing(grouping_factor)){
      results = .run_test(taxa_smry_df, metadata_map, type_header, test_type, 
                               grouping_factor, smry_fun = smry_fun)
    }
    else if(!missing(custom_test_function)){
      results = .run_test(taxa_smry_df, metadata_map, type_header, test_type,
                               cust_test = custom_test_function, 
                               smry_fun = smry_fun)
    } 
    else{
      results = .run_test(taxa_smry_df, metadata_map, type_header, test_type,
                               smry_fun = smry_fun)
    }
    # Sort by pvalues 
    results = results[with(results, order(pvals)), ]
  }
  # output data
  if(!missing(out_fp)){
    write.table(results, file = out_fp, sep = ",", row.names = TRUE, 
                col.names = NA)  
  } else results
}