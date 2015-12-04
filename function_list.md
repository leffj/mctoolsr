# mctoolsr Function List

## Load external data

`load_2_dms()`: Load two dissimilarity matrices for use with mctoolsr

`load_dm()`: Load a dissimilarity matrix for use with mctoolsr

`load_taxa_table()`: Load a taxa table (e.g. an OTU table) with corresponding metadata for use with mctoolsr


## Taxa table manipulation

`export_otu_table()`: Export an OTU table as a text file

`filter_data()`: Filter samples from dataset

`filter_taxa_from_input()`: Filter taxa from a loaded dataset

`match_datasets()`: Match up samples from two datasets

`single_rarefy()`: Rarefy samples in a taxa table

`convert_to_relative_abundances()`: Convert taxon table to relative abundances


## Diversity (Alpha diversity)

`calc_diversity()`: Calculate diversity values


## Dissimilarity calculation & manipulation

`calc_dm()`: Calculate a dissimilarity matrix from a taxa table

`filter_dm()`: Filter samples from dissimilarity matrix

`calc_mean_dissimilarities()`: Calculate mean dissimilarities using a metadata factor

`convert_dm_to_3_column()`: Convert dissimilarity matrix to 3 column format

`add_metadata_to_dm_clmns()`: Add metadata to an additional column in column formatted dissimilarities dataframe


## Disimilarity-based Analyses

`calc_ordination()`: Calculate point coordinates in an ordination

`calc_pairwise_permanovas()`: Calculate pairwise PERMANOVA results

`calc_taxa_changes()`: Calculate the changes in taxonomic relative abundances compared to controls


## Taxonomy-based Analyses

`summarize_taxonomy()`: Calculate values for coarser taxonomic groupings

`taxa_summary_by_sample_type()`: Further summarize output from summarize_taxonomy by sample type


## Misc. Analyses

`return_top_taxa()`: Return the most abundant taxa in a dataset


## Plots

`plot_diversity()`: Plot diversity values

`plot_nmds()`: Generate an NMDS plot quickly

`plot_ordination()`: Generate an ordination plot

`plot_taxa_bars()`: Plot stacked bar plots to represent taxa compompositions

`plot_ts_heatmap()`: Plot taxa summary heatmap

`plot_venn_diagram()`: Plot venn diagrams

`plot_dendrogram()`: Generate a dendrogram based on a dissimilarity matrix


## Example Data

`fruits_veggies()`: Fruits and Vegetables Bacterial Community Data

`sumtax_example_data()`: Taxa Summary from Produce Dataset
