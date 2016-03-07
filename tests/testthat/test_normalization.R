library(mctoolsr)
context("Normalization")


test_that("Calculating taxa means produces correct values.", {
  sts = fruits_veggies$map_loaded$Sample_type
  test_otu = t(fruits_veggies$data_loaded['OTU_508', ])
  ref_means = as.vector(tapply(test_otu, sts, mean))
  means_by_st = calc_taxa_means(input = fruits_veggies, 
                                summarize_by_factor = 'Sample_type')
  test_means = as.numeric(means_by_st$data_loaded['OTU_508',])
  expect_equal(test_means, ref_means)
  # test using different input format
  means_by_st = calc_taxa_means(input = fruits_veggies$data_loaded, 
                                summarize_by_factor = 'Sample_type', 
                                metadata_map = fruits_veggies$map_loaded)
  test_means = as.numeric(means_by_st$data_loaded['OTU_508',])
  expect_equal(test_means, ref_means)
})

