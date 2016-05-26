library(mctoolsr)
context("Input/output")

map_fp = system.file('extdata', 'fruits_veggies_metadata.txt', 
                     package = 'mctoolsr')
map_fp_2col = system.file('extdata', 'fruits_veggies_metadata_2col.txt', 
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
  # Loads a 2 column mapping file
  tmp = load_taxa_table(tab_fp = tab_fp_biom, map_fp = map_fp_2col, 
                        filter_cat = 'Sample_type', filter_vals = 'Mushrooms')
  expect_identical(class(tmp$map_loaded), 'data.frame')
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
  if (!is.null(tmp$taxonomy_loaded)) {
    expect_identical(row.names(tmp$data_loaded), row.names(tmp$taxonomy_loaded))
  }
  # Loading the correct number of sequences per sample
  expect_equal(as.numeric(colSums(tmp$data_loaded)[1:5]),
               c(1199, 2819, 3390, 3291, 2312))
})

test_that("Dissimilarity matrix loads without error.", {
  expect_that(load_dm(dm_fp = dm_fp, map_fp = map_fp), not(throws_error()))
})

test_that("Matched datasets have the same order of samples.", {
  fv_tt_mod = fruits_veggies$data_loaded[, c(8, 1:3)]
  fv_mod = .match_data_components(fv_tt_mod, fruits_veggies$map_loaded, 
                                  fruits_veggies$taxonomy_loaded)
  matched = match_datasets(fv_mod, fruits_veggies)
  expect_identical(row.names(matched$ds1$map_loaded), 
                   row.names(matched$ds2$map_loaded))
})

test_that("Matched datasets have the same order of taxa.", {
  fv_tt_mod = fruits_veggies$data_loaded[, c(8, 1:3)]
  fv_mod = .match_data_components(fv_tt_mod, fruits_veggies$map_loaded, 
                                  fruits_veggies$taxonomy_loaded)
  fv_mod = filter_taxa_from_input(fv_mod, 1)
  matched = match_datasets(fv_mod, fruits_veggies, match_taxa = TRUE)
  expect_identical(row.names(matched$ds1$data_loaded), 
                   row.names(matched$ds2$data_loaded))
})

test_that("Taxa table can be exported and reloaded.", {
  basedir = system.file('extdata', package = 'mctoolsr')
  outfp = paste0(basedir, '/exported_taxa_table.txt')
  mapoutfp = paste0(basedir, '/exported_map_file.txt')
  export_taxa_table(fruits_veggies, outfp, mapoutfp)
  expect_that(load_taxa_table(tab_fp = outfp, map_fp = mapoutfp),
              not(throws_error()))
  unlink(outfp)
  unlink(mapoutfp)
})

test_that("Taxa table can be exported and reloaded when sample IDs start with
          number.", {
            basedir = system.file('extdata', package = 'mctoolsr')
            outfp = paste0(basedir, '/exported_taxa_table.txt')
            mapoutfp = paste0(basedir, '/exported_map_file.txt')
            # add number to beginning of sample IDs
            names(fruits_veggies$data_loaded) =
              sapply(names(fruits_veggies$data_loaded), function(x)
                paste0('18', x))
            export_taxa_table(fruits_veggies, outfp, mapoutfp)
            expect_that(load_taxa_table(tab_fp = outfp, map_fp = mapoutfp),
                        not(throws_error()))
            unlink(outfp)
            unlink(mapoutfp)
          })
