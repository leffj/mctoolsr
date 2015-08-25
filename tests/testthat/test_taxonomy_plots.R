library(mctoolsr)
context("Taxonomy plots")

map_fp = system.file('extdata', 'fruits_veggies_metadata.txt', 
                     package = 'mctoolsr')
tab_fp_biom = system.file('extdata', 'fruits_veggies_taxa_table_wTax.biom', 
                     package = 'mctoolsr')
tab_fp_txt = system.file('extdata', 'fruits_veggies_taxa_table_wTax.txt', 
                          package = 'mctoolsr')

test_that("Heatmap figure is generated without error.", {
  expect_that(plot_ts_heatmap(sumtax_example_data, fruits_veggies$map_loaded, 
                              0.05, 'Sample_type', scale_by = 'all'), 
              not(throws_error()))
})

