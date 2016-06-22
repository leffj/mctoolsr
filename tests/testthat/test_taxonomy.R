library(mctoolsr)
context("Taxonomy functions")

test_that("Taxa summary is produced correctly.", {
  # samples sum to 1 if relative abundances
  expect_equal(as.numeric(colSums(summarize_taxonomy(fruits_veggies, 4, 
                                                     relative = TRUE, 
                                          report_higher_tax = FALSE))), 
              rep(1, ncol(fruits_veggies$data_loaded)))
  # test case
  tmp = summarize_taxonomy(fruits_veggies, 2)
  expect_equal(tmp[grep('Verruco', row.names(tmp)), 'ProA15'], 2 / 3291)
})

test_that("Filtering taxonomy from table works as expected.", {
  tmp = summarize_taxonomy(fruits_veggies, 2, report_higher_tax = FALSE)
  expect_equal(row.names(filter_taxa_from_table(tmp, 0.1)), 
               c('p__Bacteroidetes', 'p__Firmicutes', 'p__Proteobacteria'))
  expect_equal(row.names(filter_taxa_from_table(tmp, taxa_to_keep = 
                                                  'p__Proteobacteria')), 
                         'p__Proteobacteria')
  expect_equal(row.names(filter_taxa_from_table(tmp, filter_thresh = 0.1, 
                                                taxa_to_remove = 
                                                  'p__Proteobacteria')), 
               c('p__Bacteroidetes', 'p__Firmicutes'))
  expect_equal(row.names(filter_taxa_from_table(tmp, taxa_to_keep = 
                                                  'p__Proteobacteria', 
                                                taxa_to_remove = 
                                                  'p__Proteobacteria')), 
               character(0))
})

test_that("Filtering taxonomy from input dataset works as expected.", {
  fruits_veggies_rar = single_rarefy(fruits_veggies, 100)
  expect_equal(row.names(filter_taxa_from_input(fruits_veggies_rar, 
                                                filter_thresh = 10)$data_loaded), 
               c('OTU_5339', 'OTU_6975'))
  expect_equal({
    tmp = filter_taxa_from_input(fruits_veggies_rar, 
                                 taxa_to_keep = 'Proteobacteria')
    row.names(summarize_taxonomy(tmp, 2, report_higher_tax = FALSE))},
    'p__Proteobacteria')
  # filter OTU IDs
  expect_equal({
    tmp = filter_taxa_from_input(fruits_veggies, 
                                 taxa_IDs_to_keep = c('OTU_3', 'OTU_64'))
    row.names(tmp$data_loaded)},
    c('OTU_3', 'OTU_64'))
  expect_equal({
    tmp = filter_taxa_from_input(fruits_veggies, 5,
                                 taxa_IDs_to_keep = c('OTU_3', 'OTU_64'))
    row.names(tmp$data_loaded)},
    c('OTU_64'))
  expect_equal({
    tmp = filter_taxa_from_input(fruits_veggies, 10, taxa_to_keep = 'Proteoba', 
                                 taxa_IDs_to_remove = 'OTU_64')
    nrow(tmp$data_loaded)},
    14)
})

test_that("Filtering taxonomy from taxa summary works as expected.", {
  ts = summarize_taxonomy(fruits_veggies, 2)
  only_proteos = filter_taxa_from_table(ts, filter_thresh = 0.5)
  expect_equal(row.names(only_proteos), 'k__Bacteria; p__Proteobacteria')
  firm_act = filter_taxa_from_table(ts, taxa_to_keep = 
                                      c('Firmicutes', 'Actinobacteria'))
  expect_equal(row.names(firm_act), c("k__Bacteria; p__Actinobacteria", 
                                      "k__Bacteria; p__Firmicutes"))
  no_firm_act = filter_taxa_from_table(ts, taxa_to_remove = 
                                         c('Firmicutes', 'Actinobacteria'))
  expect_true(all(c('Firmicutes', 'Actinobacteria') %in% 
                    row.names(no_firm_act) == FALSE))
  # try filtering multiple taxa with one that doesn't exist
  filtered = suppressWarnings(filter_taxa_from_table(
    tax_table = ts,
    taxa_to_keep = c('Proteobacteria',
                     'Firmites',
                     'Actinobacteria')
  ))
  expect_equal(
    row.names(filtered),
    c(
      "k__Bacteria; p__Actinobacteria",
      "k__Bacteria; p__Proteobacteria"
    )
  )
})

test_that("Taxa changes calculated correctly.", {
  ts = summarize_taxonomy(fruits_veggies, 2)
  ts_filt = filter_taxa_from_table(ts, filter_thresh = 0.01)
  st = fruits_veggies$map_loaded$Sample_type
  ft = fruits_veggies$map_loaded$Farm_type
  prot = t(ts_filt['k__Bacteria; p__Proteobacteria', ])[,1]
  testdf = data.frame(st, ft, prot)[order(st, ft), ]
  testdf_straw = testdf[testdf$st == 'Strawberries',]
  conv_mean = tapply(testdf_straw$prot, testdf_straw$ft, mean)['Conventional']
  org = testdf_straw[testdf_straw$ft == 'Organic',]
  testvals = (org$prot - conv_mean) / conv_mean * 100
  row.names(org)
  calcvals_df = calc_taxa_changes(
    ts = ts_filt, metadata_map = fruits_veggies$map_loaded,
    block_header = 'Sample_type',
    treatment_header = 'Farm_type',
    control_label = 'Conventional')
  calcvals = dplyr::filter(
    calcvals_df, block == 'Strawberries', Tx == 'Organic',
    Taxon == 'k__Bacteria; p__Proteobacteria'
  )$pct_change
  expect_identical(calcvals, testvals)
})
