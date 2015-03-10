#################################################################################
### R code to find the taxa driving differences between microbial communities ###
### given a taxa summary table and mapping file from QIIME                    ###
### -- Jon Leff -- January 16, 2015 --                                        ###
#################################################################################

# This code will: (1) Filter the taxa summary to remove taxa that do not meet
# an abundance threshold in any factor level. This is based on mean abundance.
# (2) Calculate which taxa have differences in relative abundance among factor
# levels. This is based on either Mann-Whitney tests or Kruskal-Wallis (both
# non-parametric tests). Use Mann-Whitney for 2 factor levels and K-W for more
# than 2. (3) Output results including adjusted (Bonferroni and FDR) p-values
# and means.

#### Functions needed for analysis:      ####
#### Run all so that they will be usable ####
#### Example usage is at bottom of file  ####
library(nlme)

# get metadata values for a specific variable in the same order as the samples
# in the taxa table
get_metadata = function(t_table,map_file,variable) map_file[match(names(t_table),row.names(map_file)),variable]

# function to filter taxa
filter_taxa_dit = function(t_table,map_file,f_level,f_factor){
  # Check if the t_table only has one sample
  if(class(t_table)=="numeric"){
    "skip"
  } else {
    factorMeta = get_metadata(t_table,map_file,f_factor)
    rowsToKeep = c()
    for(i in 1:nrow(t_table)){
      #   in the row, calculate means for each factor level and keep if one is 
      #   greater than filter
      meanAbunds <- NULL
      meanAbunds <- aggregate(as.numeric(t(t_table[i,])),list(factorMeta),mean)
      if(max(meanAbunds$x) >= f_level){
        rowsToKeep <- c(rowsToKeep,i)
      }
    }
    t_table[rowsToKeep,]
  }
}

# run Wilcoxon Rank-Sum test (Mann-Whitney U test) and return p-value
run_MW_test = function(dependent,factor){
  # check for only two factor levels
  if(length(unique(factor))!=2) print('Mann-Whitney test requires exacly two factor levels.')
  wilcox.test(formula=dependent~factor)$p.value
}

# run Kruskal-Wallis test
run_KW_test = function(dependent,factor){
  kruskal.test(formula=dependent~factor)$p.value
}

# run 2-way NP test
#trying a rank transformation then lme model
run_2WNP_test = function(dependent,factor,g_factor){
  dep.ranked=rank(dependent,ties.method='average')
  model=lme(fixed=dep.ranked~factor,random=~1|g_factor)
  anova(model)[["p-value"]][2]
}

# run statistical test (Mann-whitney, Kruskal-Wallis, or 2-way NP) on each taxon
# in a provided taxa table
run_test = function(t_table,map_file,fctr,type,g_fctr){
  fctrMeta = as.factor(as.vector(get_metadata(t_table,map_file,fctr)))
  if(!missing(g_fctr)) gfctrMeta = as.factor(as.vector(get_metadata(t_table,map_file,g_fctr)))
  pvals = c()
  for(i in 1:nrow(t_table)){
    if(type=='MW') pvals = c(pvals,run_MW_test(as.vector(t(t_table[i,])),fctrMeta))
    else if(type=='KW') pvals = c(pvals,run_KW_test(as.vector(t(t_table[i,])),fctrMeta))
    else if(type=='2WNP') pvals = c(pvals,run_2WNP_test(as.vector(t(t_table[i,])),fctrMeta,gfctrMeta))
    else print('Invalid test type specified')
    if(i==1){
      meanAbunds = aggregate(as.numeric(t(t_table[i,])),list(fctrMeta),mean)
    } else{
      means = aggregate(as.numeric(t(t_table[i,])),list(fctrMeta),mean)[,2]
      meanAbunds = cbind(meanAbunds,means)
    }
  }
  # generate bonforroni corrected pvals
  pvalsBon = pvals*length(pvals)
  # generate FDR corrected pvals (taken from otu_category_significance.py)
  # Ranks p-values low to high and multiplies each p-value by the number of
  # comparisons divided by the rank.
  pvalsFDR = pvals*(length(pvals)/rank(pvals,ties.method="average"))
  # prep means to be added
  factorLevels = as.character(meanAbunds[,1])
  meanAbunds[,1] = NULL
  # make result df
  result = as.data.frame(cbind(pvals,pvalsBon,pvalsFDR,t(meanAbunds)))
  row.names(result) = row.names(t_table)
  colnames(result) = c("pvals","pvalsBon","pvalsFDR",factorLevels)
  result
}

# filter out blanks code (not currently used)
# if(omitBlanks){
#   taxa_table <- taxa_table[,factorMeta!=""]
#   factorMeta <- as.factor(as.character(factorMeta[factorMeta!=""]))
# }

add_taxonomy = function(results_table, taxonomy_fp){
  tax = read.table(taxonomy_fp,header=FALSE,sep="\t",row.names=1,comment.char="",check.names=FALSE)
  matchedTax = tax[match(row.names(results_table), row.names(tax)),]
  results_table.new = cbind(results_table, taxonomy = matchedTax[,1])
  results_table.new
}

# function to show contributions of specific taxa to variation among communities
# using Mann-Whitney (2 factor levels) or Kruskal-Wallis (more than 2) tests
# PARAMETERS:
# ts_fp=taxa summary filepath
# map_fp=mapping file filepath
# out_fp=test results output filepath
# factor=mapping file header (in quotation marks) of factor for which you are testing for differences
# filterLevel=number from 0 to 1--the minimum mean relative abundance needed in at least one of the 
#             factor levels for a taxon to be retained in the analysis
# testType=either 'MW' or 'KW' (i.e. Wilcoxon/Mann-Whitney U for 2 factor levels or Kruskal-Wallis 
#          for more than two factor levels)
differences_in_taxa1 = function(ts_fp,map_fp,out_fp,factor,filterLevel,testType,gfactor,tax_fp){
  # import taxa summary and mapping file
  ts = read.table(ts_fp,header=TRUE,sep="\t",row.names=1,comment.char="",check.names=FALSE)
  map = read.table(map_fp,header=TRUE,sep="\t",row.names=1,comment.char="",check.names=FALSE)
  # match up data from both
  samplesInBoth=intersect(row.names(map),names(ts))
  ts.use=ts[,match(samplesInBoth,names(ts))]
  map.use=map[match(samplesInBoth,row.names(map)),]
  # filter taxa summary table by abundance in any/either factor level
  taxa.use.filt <- filter_taxa_dit(ts.use,map.use,filterLevel,factor)
  if(!missing(gfactor)){
    testResults <- run_test(taxa.use.filt,map.use,factor,testType,gfactor)
  }
  else testResults <- run_test(taxa.use.filt,map.use,factor,testType)
  # Sort by pvalues 
  testResults <- testResults[with(testResults,order(pvals)),]
  # Optionally add taxonomy
  if(!missing(tax_fp)){
    testResults <- add_taxonomy(testResults, tax_fp)
  }
  # output data
  write.table(x=testResults,file=out_fp,sep="\t",row.names=TRUE,col.names=NA)
}

differences_in_taxa2 = function(data, out_fp, factor, filterLevel, testType, gfactor, tax_fp){
  # import taxa summary and mapping file
  ts = data$data_loaded
  map = data$map_loaded
  # filter taxa summary table by abundance in any/either factor level
  taxa.filt = filter_taxa_dit(ts, map, filterLevel, factor)
  if(!missing(gfactor)){
    testResults = run_test(taxa.filt, map, factor, testType, gfactor)
  }
  else testResults = run_test(taxa.filt, map, factor, testType)
  # Sort by pvalues 
  testResults = testResults[with(testResults, order(pvals)),]
  # Optionally add taxonomy
  if(!missing(tax_fp)){
    testResults = add_taxonomy(testResults, tax_fp)
  }
  # output data
  write.table(x=testResults, file=out_fp, sep="\t", row.names=TRUE, col.names=NA)
}

differences_in_taxa = function(data, out_fp, factor, filterLevel, testType, gfactor, tax_fp, ts_fp, map_fp){
  if(!missing(data)) differences_in_taxa2(data, out_fp, factor, filterLevel, testType, gfactor, tax_fp)
  else differences_in_taxa1(ts_fp,map_fp,out_fp,factor,filterLevel,testType,gfactor,tax_fp)
}

#######################
#### Example usage ####
#######################

# settings for all cultivars together
# map_fp='/Users/leffj/Box\ Sync/sequencing_run_2014-08-12/sym_map_2014-08-12.txt'
# factor='crop_type'
# filterLevel=1
# testType='KW'
# # gfactor=''
# # bact
# taxonomy_fp='/Users/leffj/Box\ Sync/sequencing_run_2014-08-12/16S/clustering_output/rdp_assigned_taxonomy/rep_set_numbered_tax_assignments.txt'
# ts_fp='/Users/leffj/Box\ Sync/sequencing_run_2014-08-12/16S/clustering_output/dit/otu_table_wTax_noChloroMito_1000_dit.txt'
# out_fp='/Users/leffj/Box\ Sync/sequencing_run_2014-08-12/16S/clustering_output/dit/taxa_differences.txt'
# differences_in_taxa(ts_fp,map_fp,out_fp,factor,filterLevel,testType,
#                     tax_fp=taxonomy_fp)
