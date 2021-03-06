library(mctoolsr)
context("Miscellaneous analyses")

test_that("Core taxa calculation is correct.", {
  test_ds = data.frame(samp1 = c(10, 0, 0, 9, 1291, 0), 
                       samp2 = c(0, 5, 789, 0, 0, 1), 
                       samp3 = c(5, 6, 843, 0, 87, 4),
                       samp4 = c(0, 0, 0, 0, 0, 101), 
                       samp5 = c(0, 5, 6, 0, 0, 0),
                       samp6 = c(0, 20, 1, 0, 0 , 0), 
                       row.names = paste0('OTU_', seq(6)))
  test_input = list(data_loaded = test_ds, 
                    map_loaded = data.frame(row.names = colnames(test_ds), 
                                            type = c('one', 'one', 'one', 'two', 
                                                     'two', 'two')))
  core_otus = c('OTU_2', 'OTU_3')
  core_out = core_taxa(input = test_input, type_header = 'type')
  expect_equal(core_otus, as.character(core_out$OTU_ID))
})

test_that("Prop shared calculation is correct.", {
  test_ds = data.frame(samp1 = c(10, 0, 0, 9, 1291, 0), 
                       samp2 = c(0, 5, 789, 0, 0, 1), 
                       samp3 = c(5, 6, 843, 0, 87, 4),
                       samp4 = c(0, 0, 0, 0, 0, 101), 
                       row.names = paste0('OTU_', seq(6)))
  pairs = as.data.frame(t(combn(colnames(test_ds), 2)))
  paired_data = apply(pairs, 1, function(x) test_ds[, x])
  no_taxa = sapply(paired_data, function(x) sum(x[, 1] > 0 | x[, 2] > 0))
  shared_taxa = sapply(paired_data, function(x) sum(x[, 1] * x[, 2] > 0))
  prop_shared = shared_taxa / no_taxa
  mean_prop_shared = mean(prop_shared)
  test_input = list(data_loaded = test_ds, 
                    map_loaded = data.frame(row.names = colnames(test_ds), 
                                            type = c('one', 'one', 'two', 
                                                     'two')))
  # across all samples
  expect_equal(calc_prop_shared_taxa(test_input), mean_prop_shared)
  # within between sample types
  pairs$between_types = c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE)
  mean_prop_shared_bt = mean(prop_shared[pairs$between_types])
  expect_equal(calc_prop_shared_taxa(test_input, 'type', c('one', 'two')), 
               mean_prop_shared_bt)
})


test_that("Prop taxa of another sample type calculation is correct.", {
  test_ds = data.frame(samp1 = c(10, 0, 0, 9, 1291, 0), 
                       samp2 = c(0, 5, 789, 0, 0, 1), 
                       samp3 = c(5, 6, 843, 0, 87, 4),
                       samp4 = c(0, 0, 0, 0, 0, 101), 
                       row.names = paste0('OTU_', seq(6)))
  test_input = list(data_loaded = test_ds, 
                    map_loaded = data.frame(row.names = colnames(test_ds), 
                                            type = c('one', 'one', 'two', 
                                                     'two')))
  # proportion of taxa in type two is also in type one
  # prop samp3 also in 1: 2/5
  # prop samp3 also in 2: 3/5
  # prop samp4 also in 1: 0/1
  # prop samp4 also in 2: 1/1
  expected_result = mean(c(2/5, 3/5, 0, 1/1))
  expect_equal(calc_prop_taxa_from_sample_type(test_input, 'type', 'two', 
                                               'one'), 
               expected_result)
})
