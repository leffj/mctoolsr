# Getting Started with mctoolsr
Jonathan W. Leff  
2015-12-03  

**mctoolsr** (prononounced, M-C-tools-R) is an **R** package developed to facilitate microbial community analyses. The current functions are meant to handle an input taxa (OTU) table in either biom or tab-delimited format and help streamline common (and more specialized) downstream analyses. It is under active development, so please submit bug reports and feature requests as indicated below.

This document serves as a brief introduction to using **mctoolsr**. This document will go through getting the package working and a few examples using the most popular functions.

**See here for a list of current functions**: https://github.com/leffj/mctoolsr/blob/master/function_list.md

### Getting and using **mctoolsr**

**mctoolsr** is available on Github at: https://github.com/leffj/mctoolsr

To install:


```r
install.packages("devtools")
devtools::install_github("leffj/mctoolsr")
```

To use in a script, load and attach the package using:

```r
library(mctoolsr)
```


### Updating **mctoolsr**

**mctoolsr** is under active devlopment. Update regularly to use the latest features and for bug fixes. You will have to load and attach as above after updating.


```r
devtools::install_github("leffj/mctoolsr")
```

### To report bugs or request features, go to:
https://github.com/leffj/mctoolsr/issues


### Citing
The package is in initial development, so for now, please simply cite the Github URL: `https://github.com/leffj/mctoolsr/`


### Examples

Note that the following examples use an example dataset taken from a study examining the bacterial communities associated with fruits and vegetables ([Leff et al. 2013][ref1]). You can find this dataset in the mctoolsr/examples directory.    


#### Loading OTU tables and metadata

You can load a taxa (i.e. OTU) table in biom format using the following approach. Note that you don't need to use `system.file()` when using your own filepaths -- just use the filepath directly. This function just helps your system find the example files.

One of the nice things about loading your data this way is that all the sample IDs will be matched between your taxon table and metadata so that they will be in the same order and any sample IDs not present in one or the other will be dropped.

You can optionally filter out samples of a specific type during this step, but this can also be done separately as shown here.


```r
tax_table_fp = system.file('extdata', 'fruits_veggies_taxa_table_wTax.biom', 
                     package = 'mctoolsr')
map_fp = system.file('extdata', 'fruits_veggies_metadata.txt', 
                     package = 'mctoolsr')
input = load_taxa_table(tax_table_fp, map_fp)
```

```
## 32 samples loaded
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
## ProA37 ProB70 ProC66 ProB39 ProC40 ProB57 ProA12 ProA58 ProB34 ProC36 
##   1009   1011   1068   1152   1179   1192   1199   1216   1265   1313 
## ProA65  ProB9 ProB60 ProC65 ProB35 ProB10 ProB67 ProA66 ProB36 ProB71 
##   1367   1371   1395   1409   1492   1599   1611   1614   1642   1745 
## ProB40 ProA36 ProB58 ProB33 ProA16 ProA35 ProA33 ProB12 ProA13 ProA34 
##   1771   1802   1860   2257   2312   2530   2585   2642   2819   2982 
## ProA15 ProA14 
##   3291   3390
```


#### Rarefying

As you can see from the previous example, we can rarefy (i.e. normalize for variable sequence depths) to 1000 sequences per sample without losing any samples. This can be done using the following command:


```r
input_rar = single_rarefy(input, 1000)
```

```
## 32 samples remaining
```

```r
colSums(input_rar$data_loaded)
```

```
## ProA12 ProA13 ProA14 ProA15 ProA16 ProA33 ProA34 ProA35 ProA36 ProA37 
##   1000   1000   1000   1000   1000   1000   1000   1000   1000   1000 
## ProA58 ProA65 ProA66 ProB10 ProB12 ProB33 ProB34 ProB35 ProB36 ProB39 
##   1000   1000   1000   1000   1000   1000   1000   1000   1000   1000 
## ProB40 ProB57 ProB58 ProB60 ProB67 ProB70 ProB71  ProB9 ProC36 ProC40 
##   1000   1000   1000   1000   1000   1000   1000   1000   1000   1000 
## ProC65 ProC66 
##   1000   1000
```


#### Summarize taxonomic relative abundances at a higher taxonomic level

It is useful to get a feel for the taxonomic composition of your samples early on in the exploratory data analysis process. This can quickly be done by calculating taxonomic summaries at higher taxonomic levels - in this case at the phylum level. The values represent the sum of all the relative abundances for OTUs classified as belonging to the indicated phylum. In this example just the first few phyla and samples are shown.


```r
tax_sum_phyla = summarize_taxonomy(input_rar, level = 2, report_higher_tax = FALSE)
tax_sum_phyla[1:5, 1:8]
```

```
##                    ProA12 ProA13 ProA14 ProA15 ProA16 ProA33 ProA34 ProA35
## p__                 0.001  0.000  0.000  0.000  0.000   0.00  0.000  0.000
## p__[Thermi]         0.001  0.000  0.000  0.000  0.000   0.00  0.007  0.000
## p__Acidobacteria    0.001  0.000  0.000  0.000  0.000   0.00  0.000  0.000
## p__Actinobacteria   0.392  0.006  0.002  0.018  0.002   0.05  0.068  0.087
## p__Armatimonadetes  0.000  0.001  0.000  0.002  0.000   0.00  0.000  0.000
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
## Run 0 stress 0.112383 
## Run 1 stress 0.1123835 
## ... procrustes: rmse 0.0008635581  max resid 0.00344912 
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
```

```
## 28 samples remaining
```

```r
dm = calc_dm(input_rar_filt$data_loaded)
ord = calc_ordination(dm, 'nmds')
```

```
## Run 0 stress 0.101338 
## Run 1 stress 0.1013378 
## ... New best solution
## ... procrustes: rmse 0.0003283851  max resid 0.001170659 
## *** Solution reached
```

```r
plot_ordination(input_rar_filt, ord, 'Sample_type', 'Farm_type', hulls = TRUE)
```

![](README_files/figure-html/unnamed-chunk-11-1.png) 


#### Filtering taxa

There are multiple taxa filtering options in **mctoolsr**. This example shows how to explore the proteobacteria sequences across the samples. Taxa can also be filtered based on their relative abundance.


```r
input_proteobact = filter_taxa_from_input(input, taxa_to_keep = 'p__Proteobacteria')
```

```
## 2265 taxa removed
```

```r
sort(colSums(input_proteobact$data_loaded))
```

```
## ProC40 ProA14 ProA12 ProA58 ProA16 ProB70 ProB10  ProB9 ProA37 ProB12 
##    231    435    459    614    780    843    889    940    980   1019 
## ProC66 ProC36 ProB57 ProA13 ProA36 ProB39 ProB34 ProB60 ProA65 ProA35 
##   1032   1051   1063   1080   1082   1082   1196   1196   1261   1336 
## ProC65 ProB71 ProA34 ProB67 ProA66 ProB35 ProB58 ProB36 ProA15 ProB40 
##   1369   1376   1406   1422   1436   1472   1544   1594   1601   1748 
## ProA33 ProB33 
##   1768   2126
```

```r
input_proteobact_rar = single_rarefy(input_proteobact, 219)
```

```
## 32 samples remaining
```

```r
plot_nmds(calc_dm(input_proteobact_rar$data_loaded), metadata_map = input_proteobact_rar$map_loaded, 
          color_cat = 'Sample_type')
```

```
## Run 0 stress 0.1400091 
## Run 1 stress 0.1363107 
## ... New best solution
## ... procrustes: rmse 0.02851887  max resid 0.1416946 
## Run 2 stress 0.1400091 
## Run 3 stress 0.1363107 
## ... New best solution
## ... procrustes: rmse 2.690236e-05  max resid 8.827398e-05 
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
## f__Pseudomonadaceae    1.573687e-05 0.0001101581 0.0001101581  0.146000
## f__Sphingobacteriaceae 7.765991e-05 0.0005436194 0.0002718097  0.288000
## f__Enterobacteriaceae  3.830648e-04 0.0026814537 0.0008938179  0.033125
## f__[Weeksellaceae]     3.962317e-04 0.0027736221 0.0006934055  0.106375
## unclassified           4.702537e-03 0.0329177617 0.0065835523  0.072750
## f__Bacillaceae         1.312857e-02 0.0918999590 0.0153166598  0.002750
## f__Sphingomonadaceae   4.606166e-02 0.3224316222 0.0460616603  0.007125
##                             Spinach Strawberries
## f__Pseudomonadaceae    0.0618571429 0.0012307692
## f__Sphingobacteriaceae 0.0028571429 0.0007692308
## f__Enterobacteriaceae  0.7380000000 0.6155384615
## f__[Weeksellaceae]     0.0025714286 0.0016923077
## unclassified           0.0282857143 0.0231538462
## f__Bacillaceae         0.0037142857 0.1624615385
## f__Sphingomonadaceae   0.0007142857 0.0595384615
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
## Run 1 stress 0 
## ... procrustes: rmse 0.2034204  max resid 0.355068 
## Run 2 stress 0.1842128 
## Run 3 stress 0.2465821 
## Run 4 stress 0 
## ... procrustes: rmse 0.04681936  max resid 0.09369443 
## Run 5 stress 0.1842127 
## Run 6 stress 0 
## ... procrustes: rmse 0.2405042  max resid 0.3483202 
## Run 7 stress 0 
## ... procrustes: rmse 0.1811332  max resid 0.3081502 
## Run 8 stress 0.1842126 
## Run 9 stress 9.82482e-05 
## ... procrustes: rmse 0.2235992  max resid 0.303423 
## Run 10 stress 9.812367e-05 
## ... procrustes: rmse 0.2092321  max resid 0.3300858 
## Run 11 stress 9.939148e-05 
## ... procrustes: rmse 0.1922165  max resid 0.3217053 
## Run 12 stress 0 
## ... procrustes: rmse 0.04003455  max resid 0.07309192 
## Run 13 stress 0.2465469 
## Run 14 stress 0 
## ... procrustes: rmse 0.02953637  max resid 0.05566956 
## Run 15 stress 0 
## ... procrustes: rmse 0.2624305  max resid 0.3667243 
## Run 16 stress 0 
## ... procrustes: rmse 0.2237332  max resid 0.3510673 
## Run 17 stress 7.875637e-05 
## ... procrustes: rmse 0.1522022  max resid 0.2695034 
## Run 18 stress 0.1842126 
## Run 19 stress 9.203126e-05 
## ... procrustes: rmse 0.2406716  max resid 0.3134842 
## Run 20 stress 0 
## ... procrustes: rmse 0.2012379  max resid 0.2997708
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
