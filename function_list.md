# mctoolsr Function List



## Dissimilarity calculation and manipulation

`add_metadata_to_dm_clmns()`: Add metadata to an additional column in column formatted dissimilarities dataframe

`calc_dm()`: Calculate a dissimilarity matrix from a taxa table

`calc_mean_dissimilarities()`: Calculate mean dissimilarities using a metadata factor

`calc_ordination()`: Calculate Point Coordinates in an Ordination

`calc_pairwise_permanovas()`: Calculate pairwise PERMANOVA results

`convert_dm_to_3_column()`: Convert dissimilarity matrix to 3 column format

`filter_dm()`: Filter samples from dissimilarity matrix


## Diversity

`calc_diversity()`: Calculate diversity values


## Load external data

`load_2_dms()`: Load two dissimilarity matrices for use with mctoolsr

`load_dm()`: Load a dissimilarity matrix for use with mctoolsr

`load_taxa_table()`: Load a taxa table for use with mctoolsr

`load_tree()`: Load in a phylogenetic tree


## Misc analyses

`calc_prop_shared_taxa()`: Calculate mean proportions of shared taxa between paired samples

`calc_prop_taxa_from_sample_type()`: Calculate the proportion of taxa in a set of samples that are also observed in another sample type

`calc_taxa_changes()`: Calculate the changes in taxonomic relative abundances compared to controls

`core_taxa()`: Determine taxa that are common (core) across sample types

`return_top_taxa()`: Return the most abundant taxa in a dataset


## Phylogeny

`filter_tree()`: Filter tips in phylogenetic tree


## Plots

`plot_dendrogram()`: Generate a dendrogram based on a dissimilarity matrix

`plot_diversity()`: Plot diversity values

`plot_nmds()`: Generate an NMDS Plot Quickly

`plot_ordination()`: Generate an Ordination Plot

`plot_taxa_bars()`: Plot stacked bar plots to represent taxa compompositions

`plot_ts_heatmap()`: Plot taxa summary heatmap

`plot_venn_diagram()`: Plot Venn Diagrams


## Taxa table manipulation

`calc_taxa_means()`: Calculate mean taxa values across a specified factor

`export_taxa_table()`: Export a taxa table as a text file

`filter_data()`: Filter samples from dataset

`filter_samples_by_counts()`: Filter samples from a dataset based on number of sequences

`match_datasets()`: Match up samples from two datasets

`rename_samples()`: Rename samples in an mctoolsr dataset


## Taxa table normalization

`convert_to_relative_abundances()`: Convert taxa table to relative abundances

`single_rarefy()`: Rarefy samples in a taxa table


## Taxonomy-based analyses

`collapse_taxonomy()`: Collapse taxonomy dataframe to character vector

`filter_taxa_from_input()`: Filter taxa from a loaded dataset

`filter_taxa_from_table()`: Filter taxa from an individual taxa summary table

`summarize_taxonomy()`: Calculate values for coarser taxonomic groupings

`taxa_summary_by_sample_type()`: Further summarize output from summarize_taxonomy by sample type

