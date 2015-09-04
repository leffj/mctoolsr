library(mctoolsr)
context("Taxonomy functions")

test_that("Taxa summary is produced correctly.", {
  # samples sum to 1 if relative abundances
  expect_equal(as.numeric(colSums(summarize_taxonomy(fruits_veggies, 4, 
                                                     relative = TRUE, 
                                          report_higher_tax = FALSE))), 
              rep(1, ncol(fruits_veggies$data_loaded)))
  # test case
  tmp = summarize_taxonomy(fruits_veggies, 2)
  expect_equal(tmp[grep('Verruco', row.names(tmp)), 'ProA15'], 2 / 3291)
})

test_that("Filtering taxonomy from table works as expected.", {
  tmp = summarize_taxonomy(fruits_veggies, 2, report_higher_tax = FALSE)
  expect_equal(row.names(filter_taxa_from_table(tmp, 0.1)), 
               c('p__Bacteroidetes', 'p__Firmicutes', 'p__Proteobacteria'))
  expect_equal(row.names(filter_taxa_from_table(tmp, taxa_to_keep = 
                                                  'p__Proteobacteria')), 
                         'p__Proteobacteria')
  expect_equal(row.names(filter_taxa_from_table(tmp, filter_thresh = 0.1, 
                                                taxa_to_remove = 
                                                  'p__Proteobacteria')), 
               c('p__Bacteroidetes', 'p__Firmicutes'))
  expect_equal(row.names(filter_taxa_from_table(tmp, taxa_to_keep = 
                                                  'p__Proteobacteria', 
                                                taxa_to_remove = 
                                                  'p__Proteobacteria')), 
               character(0))
})

test_that("Filtering taxonomy from input dataset works as expected.", {
  fruits_veggies_rar = single_rarefy(fruits_veggies, 100)
  expect_equal(row.names(filter_taxa_from_input(fruits_veggies_rar, 
                                                filter_thresh = 10)$data_loaded), 
               c('OTU_5339', 'OTU_6975'))
  expect_equal({
    tmp = filter_taxa_from_input(fruits_veggies_rar, 
                                 taxa_to_keep = 'Proteobacteria')
    row.names(summarize_taxonomy(tmp, 2, report_higher_tax = FALSE))},
    'p__Proteobacteria')
  # filter OTU IDs
  expect_equal({
    tmp = filter_taxa_from_input(fruits_veggies, 
                                 taxa_IDs_to_keep = c('OTU_3', 'OTU_64'))
    row.names(tmp$data_loaded)},
    c('OTU_3', 'OTU_64'))
  expect_equal({
    tmp = filter_taxa_from_input(fruits_veggies, 5,
                                 taxa_IDs_to_keep = c('OTU_3', 'OTU_64'))
    row.names(tmp$data_loaded)},
    c('OTU_64'))
  expect_equal({
    tmp = filter_taxa_from_input(fruits_veggies, 10, taxa_to_keep = 'Proteoba', 
                                 taxa_IDs_to_remove = 'OTU_64')
    nrow(tmp$data_loaded)},
    14)
})



