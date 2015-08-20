library(mctoolsr)
context("Input output")

map_fp = system.file('extdata', 'fruits_veggies_metadata.txt', 
                     package = 'mctoolsr')
tab_fp_biom = system.file('extdata', 'fruits_veggies_taxa_table_wTax.biom', 
                     package = 'mctoolsr')
tab_fp_txt = system.file('extdata', 'fruits_veggies_taxa_table_wTax.txt', 
                          package = 'mctoolsr')

test_that("Example taxa table (biom) loads correctly.", {
  expect_that(load_taxa_table(tab_fp = tab_fp_biom, map_fp = map_fp), 
              not(throws_error()))
  expect_match(colnames(load_taxa_table(tab_fp = tab_fp_biom, 
                                        map_fp = map_fp)$map_loaded)[1], 
               'Sample_type')
})

test_that("Example taxa table (txt) loads correctly.", {
  expect_that(load_taxa_table(tab_fp = tab_fp_txt, map_fp = map_fp), 
              not(throws_error()))
  expect_match(colnames(load_taxa_table(tab_fp = tab_fp_txt, 
                                        map_fp = map_fp)$map_loaded)[1], 
               'Sample_type')
})