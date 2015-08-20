# Getting Started with mctoolsr
Jonathan Leff  
2015-08-20  

**mctoolsr** is an **R** package developed to facilitate microbial community analyses. The current functions are meant to handle an input taxa (OTU) table in either biom or tab-delimited format and help streamline common (and more specialized) downstream analyses. It is under active development, so please submit bug reports and feature requests as indicated below.

This document serves as a brief introduction to using **mctoolsr**. This document will go through getting the package working and a few examples using the most popular functions.

### Getting and using **mctoolsr**

**mctoolsr** is available on Github at: https://github.com/leffj/mctoolsr

To use:


```r
install.packages("devtools")
devtools::install_github("leffj/mctoolsr", ref = "package_prep")
```


```r
library(mctoolsr)
```


### Updating **mctoolsr**

**mctoolsr** is under active devlopment. Update regularly to use the latest features and for bug fixes.


```r
devtools::install_github("leffj/mctoolsr", ref = "package_prep")
```

### To report bugs or request features, go to:
https://github.com/leffj/mctoolsr/issues


### Examples

Note that the following examples use an example dataset taken from a study examining the bacterial communities associated with fruits and vegetables ([Leff et al. 2013][ref1]). You can find this dataset in the mctoolsr/examples directory.    


#### Loading OTU tables and metadata

You can load a taxa (i.e. OTU) table in biom format using the following approach. Note that filepaths are specific to your system, so they will likely have to be altered depending on where you cloned **mctoolsr** into.

One of the nice things about loading your data this way is that all the sample IDs will be matched between your taxon table and metadata so that they will be in the same order and any sample IDs not present in one or the other will be dropped.

You can optionally filter out samples of a specific type during this step, but this can also be done separately as shown here.


```r
tax_table_fp = system.file('extdata', 'fruits_veggies_taxa_table_wTax.biom', 
                     package = 'mctoolsr')
map_fp = system.file('extdata', 'fruits_veggies_metadata.txt', 
                     package = 'mctoolsr')
input = load_taxa_table(tax_table_fp, map_fp)
```

The loaded data will consist of three parts:

1. The taxon table itself: "data_loaded"
2. The metadata: "map_loaded"
3. The taxonomic classifiers (if provided in the biom file): "taxonomy_loaded"

Any of these components can be quickly accessed using the '$' sign notation as shown in the next example.  


#### Returning numbers of sequences per sample

This can be achieved simply by calculating column sums on the taxon table:


```r
sort(colSums(input$data_loaded))
```

```
## ProA37 ProB70 ProC66 ProB39 ProC40 ProB57 ProA12 ProC38 ProA58 ProB34 
##   1009   1011   1068   1152   1179   1192   1199   1211   1216   1265 
## ProC36 ProA65  ProB9 ProB60 ProC65 ProB35 ProB10 ProB67 ProA66 ProB36 
##   1313   1367   1371   1395   1409   1492   1599   1611   1614   1642 
## ProB71 ProB40 ProA36 ProB58 ProB33 ProA16 ProA35 ProA33 ProB12 ProA13 
##   1745   1771   1802   1860   2257   2312   2530   2585   2642   2819 
## ProA34 ProA15 ProA14 
##   2982   3291   3390
```


#### Rarefying

As you can see from the previous example, we can rarefy (i.e. normalize for variable sequence depths) to 1000 sequences per sample without losing any samples. This can be done using the following command:


```r
input_rar = single_rarefy(input, 1000)
colSums(input_rar$data_loaded)
```

```
## ProA12 ProA13 ProA14 ProA15 ProA16 ProA33 ProA34 ProA35 ProA36 ProA37 
##   1000   1000   1000   1000   1000   1000   1000   1000   1000   1000 
## ProA58 ProA65 ProA66 ProB10 ProB12 ProB33 ProB34 ProB35 ProB36 ProB39 
##   1000   1000   1000   1000   1000   1000   1000   1000   1000   1000 
## ProB40 ProB57 ProB58 ProB60 ProB67 ProB70 ProB71  ProB9 ProC36 ProC38 
##   1000   1000   1000   1000   1000   1000   1000   1000   1000   1000 
## ProC40 ProC65 ProC66 
##   1000   1000   1000
```


#### Summarize taxonomic relative abundances at a higher taxonomic level

It is useful to get a feel for the taxonomic composition of your samples early on in the exploratory data analysis process. This can quickly be done by calculating taxonomic summaries at higher taxonomic levels - in this case at the phylum level. The values represent the sum of all the relative abundances for OTUs classified as belonging to the indicated phylum. In this example just the first few phyla and samples are shown.


```r
tax_sum_phyla = summarize_taxonomy(input_rar, level = 2, report_higher_tax = FALSE)
tax_sum_phyla[1:5, 1:8]
```

```
##                    ProA12 ProA13 ProA14 ProA15 ProA16 ProA33 ProA34 ProA35
## p__                 0.001  0.000  0.000  0.000  0.000  0.000  0.000  0.000
## p__[Thermi]         0.001  0.000  0.000  0.000  0.000  0.000  0.004  0.000
## p__Acidobacteria    0.001  0.000  0.000  0.000  0.000  0.000  0.000  0.000
## p__Actinobacteria   0.386  0.009  0.002  0.020  0.002  0.055  0.075  0.087
## p__Armatimonadetes  0.000  0.000  0.000  0.004  0.000  0.000  0.000  0.000
```


#### Calculating a dissimilarity matrix

For dissimilarity-based analyses such as ordinations and PERMANOVA, it is necessary to calculate a dissimilarity matrix. There is currently support for Bray-Curtis dissimilarities based on square-root transformed data. This is a widely used dissimilarity metric for these analyses, but others will be added as requested.


```r
dm = calc_dm(input_rar$data_loaded)
```


#### Plotting an ordination

There are two ways to plot ordinations in **mctoolsr**. The multistep way is shown here, but there is also a shortcut using the `plot_nmds()` function.


```r
ord = calc_ordination(dm, 'nmds')
```

```
## Run 0 stress 0.1182119 
## Run 1 stress 0.1182121 
## ... procrustes: rmse 9.73424e-05  max resid 0.0004377437 
## *** Solution reached
```

```r
plot_ordination(input_rar, ord, 'Sample_type', 'Farm_type', hulls = TRUE)
```

![](README_files/figure-html/unnamed-chunk-10-1.png) 


#### Filtering samples

It is easy to filter samples from your dataset in **mctoolsr**. You can specify to remove samples meeting a specified condition in the metadata or keep those samples. In the example below, lettuce samples are removed, and the ordination is plotted again.


```r
input_rar_filt = filter_data(input_rar, 'Sample_type', filter_vals = 'Lettuce')
dm = calc_dm(input_rar_filt$data_loaded)
ord = calc_ordination(dm, 'nmds')
```

```
## Run 0 stress 0.1074726 
## Run 1 stress 0.1263015 
## Run 2 stress 0.1074711 
## ... New best solution
## ... procrustes: rmse 0.0008032965  max resid 0.003320957 
## *** Solution reached
```

```r
plot_ordination(input_rar_filt, ord, 'Sample_type', 'Farm_type', hulls = TRUE)
```

![](README_files/figure-html/unnamed-chunk-11-1.png) 


#### Filtering taxa

There are multiple taxa filtering options in **mctoolsr**. This example shows how to explore the proteobacteria sequences across the samples. Taxa can also be filtered based on their relative abundance.


```r
input_proteobact = filter_taxa_from_data(input, taxa_to_keep = 'p__Proteobacteria')
sort(colSums(input_proteobact$data_loaded))
```

```
## ProC38 ProC40 ProA14 ProA12 ProA58 ProA16 ProB70 ProB10  ProB9 ProA37 
##    219    231    435    459    614    780    843    889    940    980 
## ProB12 ProC66 ProC36 ProB57 ProA13 ProA36 ProB39 ProB34 ProB60 ProA65 
##   1019   1032   1051   1063   1080   1082   1082   1196   1196   1261 
## ProA35 ProC65 ProB71 ProA34 ProB67 ProA66 ProB35 ProB58 ProB36 ProA15 
##   1336   1369   1376   1406   1422   1436   1472   1544   1594   1601 
## ProB40 ProA33 ProB33 
##   1748   1768   2126
```

```r
input_proteobact_rar = single_rarefy(input_proteobact, 219)
plot_nmds(calc_dm(input_proteobact_rar$data_loaded), metadata_map = input_proteobact_rar$map_loaded, 
          color_cat = 'Sample_type')
```

```
## Run 0 stress 0.1336329 
## Run 1 stress 0.1176976 
## ... New best solution
## ... procrustes: rmse 0.1150264  max resid 0.281871 
## Run 2 stress 0.1176979 
## ... procrustes: rmse 4.351562e-05  max resid 0.0001168773 
## *** Solution reached
```

![](README_files/figure-html/unnamed-chunk-12-1.png) 


#### Taxa based exploration

It is often useful to determine the taxa driving differences between the community compositions of different sample types. This example shows one way to do this to determine taxa driving differences between sample types.


```r
tax_sum_families = summarize_taxonomy(input_rar_filt, level = 5, report_higher_tax = FALSE)
taxa_summary_by_sample_type(tax_sum_families, input_rar_filt$map_loaded, 
                            type_header = 'Sample_type', filter_level = 0.05, test_type = 'KW')
```

```
##                               pvals     pvalsBon     pvalsFDR Mushrooms
## f__Pseudomonadaceae    9.926019e-06 6.948214e-05 6.948214e-05  0.145875
## f__Sphingobacteriaceae 5.248041e-05 3.673629e-04 1.836814e-04  0.287250
## f__[Weeksellaceae]     1.843942e-04 1.290760e-03 4.302532e-04  0.108875
## f__Enterobacteriaceae  5.606372e-04 3.924461e-03 9.811152e-04  0.031750
## unclassified           3.143748e-03 2.200624e-02 4.401247e-03  0.074250
## f__Bacillaceae         8.929942e-03 6.250959e-02 1.041827e-02  0.002250
## f__Sphingomonadaceae   2.542440e-02 1.779708e-01 2.542440e-02  0.007875
##                             Spinach Strawberries
## f__Pseudomonadaceae    0.0614285714 0.0014285714
## f__Sphingobacteriaceae 0.0028571429 0.0006428571
## f__[Weeksellaceae]     0.0020000000 0.0016428571
## f__Enterobacteriaceae  0.7322857143 0.5732142857
## unclassified           0.0272857143 0.0215714286
## f__Bacillaceae         0.0048571429 0.2032857143
## f__Sphingomonadaceae   0.0008571429 0.0605000000
```

This analysis demonstrates that Pseudomonadaceae and Sphingobacteriaceae tend to have higher relative abundances on mushrooms than spinach and strawberries. The p values are based on Kruskal-Wallis tests and two different corrections are reported to deal with the multiple comparisons (Bonferroni and FDR). Rare taxa are filtered out using the `filter_level` peramter. The values indicated under the sample types are mean relative abundances.    


#### Calculating mean dissimilarities

Sometimes it is necessary to calculate mean dissimilarities. This is important in cases where sample types are pseudoreplecated. This is not the case here, but this example demonstrates this functionality.


```r
dm = calc_dm(input_rar_filt$data_loaded)
dm_aggregated = calc_mean_dissimilarities(dm, input_rar_filt$map_loaded, 
                                          'Sample_Farming', return_map = TRUE)
ord = calc_ordination(dm_aggregated$dm, ord_type = 'nmds')
```

```
## Run 0 stress 0 
## Run 1 stress 7.255037e-05 
## ... procrustes: rmse 0.01888761  max resid 0.03154863 
## Run 2 stress 0 
## ... procrustes: rmse 0.2016847  max resid 0.3368155 
## Run 3 stress 8.453686e-05 
## ... procrustes: rmse 0.2339414  max resid 0.3256835 
## Run 4 stress 0 
## ... procrustes: rmse 0.219606  max resid 0.3304342 
## Run 5 stress 8.647223e-05 
## ... procrustes: rmse 0.2376754  max resid 0.3381609 
## Run 6 stress 0.1620205 
## Run 7 stress 0 
## ... procrustes: rmse 0.184446  max resid 0.3232259 
## Run 8 stress 0 
## ... procrustes: rmse 0.06935467  max resid 0.1218684 
## Run 9 stress 0 
## ... procrustes: rmse 0.07756823  max resid 0.1553512 
## Run 10 stress 7.703313e-05 
## ... procrustes: rmse 0.2032422  max resid 0.3169305 
## Run 11 stress 5.460904e-05 
## ... procrustes: rmse 0.1798411  max resid 0.311924 
## Run 12 stress 0 
## ... procrustes: rmse 0.1685774  max resid 0.3028243 
## Run 13 stress 8.1398e-05 
## ... procrustes: rmse 0.224312  max resid 0.3489612 
## Run 14 stress 9.595804e-05 
## ... procrustes: rmse 0.1502636  max resid 0.2768574 
## Run 15 stress 0 
## ... procrustes: rmse 0.2459439  max resid 0.3438293 
## Run 16 stress 0.1842126 
## Run 17 stress 9.992642e-05 
## ... procrustes: rmse 0.02346964  max resid 0.03950049 
## Run 18 stress 9.178795e-07 
## ... procrustes: rmse 0.2442173  max resid 0.389403 
## Run 19 stress 0 
## ... procrustes: rmse 0.07802475  max resid 0.1597618 
## Run 20 stress 0 
## ... procrustes: rmse 0.2091285  max resid 0.3589257
```

```
## Warning in vegan::metaMDS(dm, k = 2): Stress is (nearly) zero - you may
## have insufficient data
```

```r
plot_ordination(dm_aggregated, ord, color_cat = 'Sample_type')
```

![](README_files/figure-html/unnamed-chunk-14-1.png) 


#### Exporting an OTU table

OTU tables that have been edited in **mctoolsr* can be exported in text (tab delimited) format for later use. Use this function:


```r
export_otu_table(input_rar_filt, "export/path/otu_table_export.txt")
```



[ref1]: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0059310
