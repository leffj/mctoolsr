library(mctoolsr)
context("Taxonomy plots")

test_that("Heatmap figure is generated without error.", {
  expect_that(plot_ts_heatmap(sumtax_example_data, fruits_veggies$map_loaded, 
                              0.05, 'Sample_type', scale_by = 'all'), 
              not(throws_error()))
})

test_that("Taxa bar plot is generated without error.", {
  expect_that(plot_taxa_bars(sumtax_example_data, fruits_veggies$map_loaded, 
                              'Sample_type', 8), 
              not(throws_error()))
})
