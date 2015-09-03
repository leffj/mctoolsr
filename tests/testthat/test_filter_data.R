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