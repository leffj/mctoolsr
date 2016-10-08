
generate_function_list = function(package) {
  # get list of functions and concepts if available
  funs = help.search('', fields = c('name', 'concept'),
                     package = package)$matches
  # summarize and group functions without concept into Miscellaneous
  funs_smry = reshape2::dcast(funs, Name + Title ~ Field, value.var = 'Entry')
  funs_smry$concept[is.na(funs_smry$concept)] = 'Miscellaneous'
  # filter out datasets
  funs_smry = funs_smry[!funs_smry$concept %in%
                          c('dont include', 'Datasets available by data()'),]
  # Sort by concept
  funs_smry = funs_smry[order(funs_smry$concept),]
  # return with Misc last
  rows = c(
    which(funs_smry$concept != 'Miscellaneous'),
    which(funs_smry$concept == 'Miscellaneous')
  )
  funs_smry[rows, c('concept', 'Name', 'Title')]
}

# generate_function_list('mctoolsr')


## write function list markdown file
devtools::load_all(".")

fundf = generate_function_list('mctoolsr')

outdir = system.file('..', package = 'mctoolsr')
outfp = paste0(outdir, '/function_list.md')
# if(file.exists(outfp)) file.remove(outfp)

# write first line
write('# mctoolsr Function List\n\n', file = outfp)

# write functions by category

for(cat in unique(fundf$concept)) {
  write(paste0('\n## ', cat, '\n'), file = outfp, append = TRUE)
  tmp = fundf[fundf$concept == cat,]
  apply(tmp, 1, function(x) {
    write(paste0('`', x['Name'], '()`: ', x['Title'], '\n'),
          file = outfp, append = TRUE)
  })
}

