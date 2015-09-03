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