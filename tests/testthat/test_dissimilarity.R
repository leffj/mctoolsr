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
  expect_equivalent({
    known_dm = data.frame(Lettuce = c(NA, 0.9497335, 0.8986858, 0.7917143), 
                          Mushrooms = c(NA, NA, 0.9699253, 0.9411904),
                          Strawberries = c(NA, NA, NA, 0.8495972),
                          Spinach = c(NA, NA, NA, NA))
    row.names(known_dm) = c('Lettuce', 'Mushrooms', 'Strawberries', 'Spinach')
    as.dist(known_dm)
    }, {
      dm = calc_dm(fruits_veggies$data_loaded)
      round(calc_mean_dissimilarities(dm, fruits_veggies$map_loaded, 
                                'Sample_type'), 7)
    })
})

test_that("Dendrogram plots without error.", {
  expect_that(plot_dendrogram(calc_dm(fruits_veggies$data_loaded), 
                              fruits_veggies$map_loaded, "Sample_type", 
                              "Farm_type"), 
              not(throws_error()))
})

test_that("plot_nmds does not throw error when plotting.", {
  expect_that({
    capture.output(plot_nmds(
      calc_dm(fruits_veggies$data_loaded),
      fruits_veggies$map_loaded, "Sample_type", "Farm_type"
    ), file = NULL)
  },
  not(throws_error()))
  unlink('Rplots.pdf')
})

test_that("Converting to 3 column format does not throw error.", {
  expect_that(convert_dm_to_3_column(calc_dm(fruits_veggies$data_loaded)),
              not(throws_error()))
})

test_that("Metadata adds without error.", {
  expect_that(add_metadata_to_dm_clmns(convert_dm_to_3_column(calc_dm(
    fruits_veggies$data_loaded)), fruits_veggies$map_loaded, "Sample_type"),
    not(throws_error()))
})

test_that("Pairwise PERMANOVA calculations go without error.", {
  expect_that(calc_pairwise_permanovas(calc_dm(fruits_veggies$data_loaded),
                                       fruits_veggies$map_loaded, 
                                       "Sample_type"),
              not(throws_error()))
})