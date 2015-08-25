library(mctoolsr)

map_fp = system.file('extdata', 'fruits_veggies_metadata.txt', 
                     package = 'mctoolsr')
tab_fp_biom = system.file('extdata', 'fruits_veggies_taxa_table_wTax.biom', 
                          package = 'mctoolsr')

input = load_taxa_table(tab_fp = tab_fp_biom, map_fp = map_fp)

## fruits_veggies
fruits_veggies = input
devtools::use_data(fruits_veggies)

## sumtax example data
sumtax_example_data = summarize_taxonomy(input, 2)
devtools::use_data(sumtax_example_data)