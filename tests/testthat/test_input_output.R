library(mctoolsr)
context("Input/output")

map_fp = system.file('extdata', 'fruits_veggies_metadata.txt', 
                     package = 'mctoolsr')
tab_fp_biom = system.file('extdata', 'fruits_veggies_taxa_table_wTax.biom', 
                     package = 'mctoolsr')
tab_fp_txt = system.file('extdata', 'fruits_veggies_taxa_table_wTax.txt', 
                          package = 'mctoolsr')
dm_fp = system.file('extdata', 'fruits_veggies_dm.txt', package = 'mctoolsr')

test_that("Example taxa table (biom) loads correctly.", {
  # Doesn't throw error
  expect_that(load_taxa_table(tab_fp = tab_fp_biom, map_fp = map_fp), 
              not(throws_error()))
  # Load taxa table
  tmp = load_taxa_table(tab_fp = tab_fp_biom, map_fp = map_fp)
  # First column header in mapping file is correct
  expect_match(colnames(tmp$map_loaded)[1], 'Sample_type')
  # Last sample ID in mapping file is loaded
  expect_true('ProA34' %in% row.names(tmp$map_loaded))
  # Sample IDs in taxa table are same order as in mapping file
  expect_identical(names(tmp$data_loaded), row.names(tmp$map_loaded))
  # Taxa IDs in the taxonomy dataframe are same as in taxa table
  if(!is.null(tmp$taxonomy_loaded)){
    expect_identical(row.names(tmp$data_loaded), row.names(tmp$taxonomy_loaded))
  }
  # Loading the correct number of sequences per sample
  expect_equal(as.numeric(colSums(tmp$data_loaded)[1:5]), 
               c(1199, 2819, 3390, 3291, 2312))
})

test_that("Example taxa table (txt) loads correctly.", {
  # Doesn't throw error
  expect_that(load_taxa_table(tab_fp = tab_fp_txt, map_fp = map_fp), 
              not(throws_error()))
  # Load taxa table
  tmp = load_taxa_table(tab_fp = tab_fp_txt, map_fp = map_fp)
  # First column header in mapping file is correct
  expect_match(colnames(tmp$map_loaded)[1], 'Sample_type')
  # Last sample ID in mapping file is loaded
  expect_true('ProA34' %in% row.names(tmp$map_loaded))
  # Sample IDs in taxa table are same order as in mapping file
  expect_identical(names(tmp$data_loaded), row.names(tmp$map_loaded))
  # Taxa IDs in the taxonomy dataframe are same as in taxa table
  if(!is.null(tmp$taxonomy_loaded)){
    expect_identical(row.names(tmp$data_loaded), row.names(tmp$taxonomy_loaded))
  }
  # Loading the correct number of sequences per sample
  expect_equal(as.numeric(colSums(tmp$data_loaded)[1:5]), 
               c(1199, 2819, 3390, 3291, 2312))
})

# test_that("Dissimilarity matrix loads without error.", {
#   expect_that(load_dm(dm_fp = dm_fp, map_fp = map_fp), not(throws_error()))
# })

test_that("Taxa table can be exported and reloaded.", {
  basedir = system.file('extdata', package = 'mctoolsr')
  outfp = paste0(basedir, '/exported_taxa_table.txt')
  mapoutfp = paste0(basedir, '/exported_map_file.txt')
  export_otu_table(fruits_veggies, outfp, mapoutfp)
  expect_that(load_taxa_table(tab_fp = outfp, map_fp = mapoutfp), 
              not(throws_error()))
  unlink(outfp)
  unlink(mapoutfp)
})
