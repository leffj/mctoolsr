library(mctoolsr)
context("Filter data")

test_that("Filtering data by sample type removes specified samples.", {
  expect_false('Mushrooms' %in% 
                 filter_data(fruits_veggies, 'Sample_type', 
                             'Mushrooms')$map_loaded$Sample_type)
})

test_that("Filtering data by sample type keeps specified samples.", {
  expect_true('Mushrooms' %in% 
                 filter_data(fruits_veggies, 'Sample_type', 
                             keep_vals = 'Mushrooms')$map_loaded$Sample_type)
  expect_false('Lettuce' %in% 
                 filter_data(fruits_veggies, 'Sample_type', 
                             keep_vals = 'Mushrooms')$map_loaded$Sample_type)
})

test_that("Filtering samples by number of sequences removes correct samples.", {
  seq_counts = colSums(fruits_veggies$data_loaded)
  filtered = names(fruits_veggies$data_loaded)[seq_counts < 1211]
  notfiltered = names(fruits_veggies$data_loaded)[seq_counts >= 1211]
  filt_input = filter_samples_by_counts(fruits_veggies, 1211)
  expect_identical(row.names(filt_input$map_loaded), notfiltered)
  expect_equal(sum(filtered %in% row.names(filt_input$map_loaded)), 0)
})