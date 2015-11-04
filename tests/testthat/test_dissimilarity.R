library(mctoolsr)
context("Dissimilarity-based functions")

test_that("Dissimilarity calculation doesn't throw error.", {
  expect_that(calc_dm(fruits_veggies$data_loaded), not(throws_error()))
})

test_that("Ordination calculation doesn't throw error.", {
  expect_that({
    dm = calc_dm(fruits_veggies$data_loaded)
    tmp = capture.output(ord <- calc_ordination(dm, ord_type = 'nmds'))
    ord}, 
    not(throws_error()))
})

test_that("Plotting ordination doesn't throw error.", {
  expect_that({
    dm = calc_dm(fruits_veggies$data_loaded)
    tmp = capture.output(ord <- calc_ordination(dm, ord_type = 'nmds'))
    plot_ordination(fruits_veggies, ord, 'Sample_type')}, 
    not(throws_error()))
})

test_that("Calculating mean dissimilarities doesn't throw error.", {
  expect_that({
    dm = calc_dm(fruits_veggies$data_loaded)
    calc_mean_dissimilarities(dm, fruits_veggies$map_loaded, 'Sample_type')}, 
    not(throws_error()))
})

test_that("Mean dissimilarities are calculated correctly.", {
  testthat::expect_equivalent({
    known_dm = data.frame(Lettuce = c(NA, 0.9497335, 0.9015809, 0.7917143), 
                          Mushrooms = c(NA, NA, 0.9702367, 0.9411904),
                          Strawberries = c(NA, NA, NA, 0.8558376),
                          Spinach = c(NA, NA, NA, NA))
    row.names(known_dm) = c('Lettuce', 'Mushrooms', 'Strawberries', 'Spinach')
    as.dist(known_dm)
    }, {
      dm = calc_dm(fruits_veggies$data_loaded)
      round(calc_mean_dissimilarities(dm, fruits_veggies$map_loaded, 
                                'Sample_type'), 7)
    })
})