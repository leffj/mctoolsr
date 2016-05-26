library(mctoolsr)
context("Diversity calculations")

test_that("Diversity calculations are correct", {
  expect_identical(calc_diversity(fruits_veggies$data_loaded, "shannon"), 
                   vegan::diversity(fruits_veggies$data_loaded, "shannon", 2))
  expect_identical(calc_diversity(fruits_veggies$data_loaded, "simpson"),
                   vegan::diversity(fruits_veggies$data_loaded, "simpson", 2))
  expect_equal(calc_diversity(fruits_veggies$data_loaded, "richness"),
               colSums(fruits_veggies$data_loaded !=0))
})

test_that("Plot is generated properly", {
  expect_that(plot_diversity(fruits_veggies, "Sample_type", "shannon"),
              not(throws_error()))
})