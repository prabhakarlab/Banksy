
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Banksy

<!-- badges: start -->

<!-- badges: end -->

Banksy is an R package that incorporates spatial information to cluster
cells in a feature space (e.g. gene expression). Spatial information is
incorporated by averaging the features of the k nearest neighbours to
generate new ‘neighbour’ features for a given cell. This is concatenated
to the cell’s own features to generate a combined feature matrix which
is used for constructing a nearest neighbour network. Leiden clustering
is used to obtain spatially-informed clusters.

# Installation

``` r
# Install from GitHub
remotes::install_github("jleechung/Banksy")

# Install from Bioconductor (sooon)
```

# Documentation

NOTE: These links will not work at the moment since the repo is still
private. Refer to the
[vignettes](https://github.com/jleechung/Banksy/tree/main/vignettes) for
the time being.

*Banksy* comes installed with
[documentation](https://jleechung.github.io/Banksy/reference/index.html)
of main functions and their usage, along with several vignettes which
detail different use cases:

  - [Working with Banksy
    objects](https://jleechung.github.io/Banksy/articles/banksy-object.html):
    Introduction to the *BanksyObject* class which serves as a container
    for *Banksy*.

  - [Finding optimal clustering
    parameters](https://jleechung.github.io/Banksy/articles/grid-search.html):
    Illustrates a grid search of parameters which best cluster cells.

  - [Analysing multiple
    datasets](https://jleechung.github.io/Banksy/articles/multiple-datasets.html):
    Illustrates analysis of multiple spatial transcriptomic datasets.

  - [Integrating *Banksy* with other
    packages](https://jleechung.github.io/Banksy/articles/integrating-pkgs.html):
    Illustrates how users can integrate *Banksy* with packages like
    *Seurat* and *Giotto*.

# Basic Usage

## Input data

Inputs consist of an expression matrix and cell locations. Sample data
is provided with the package:

``` r
library(Banksy)

expr <- readRDS(system.file('/extdata/expression.rds', package = 'Banksy'))
locs <- readRDS(system.file('/extdata/locations.rds', package = 'Banksy'))
```

The gene expression matrix for cells should be a `matrix`:

``` r
class(expr)
#> [1] "matrix" "array"
head(expr[,1:5])
#>          cell_4 cell_5 cell_6 cell_7 cell_8
#> Slc1a2        4      2     13     26      6
#> Scn4b         3      6     17      0      1
#> Itpr1         3      0     15      2      3
#> Slc25a23      3      0      5      3      1
#> Slc1a3        1      1      0     31      2
#> Nfib          1      1      1      4      0
```

while cell locations should be supplied as a `data.frame`:

``` r
class(locs)
#> [1] "data.frame"
head(locs)
#>            sdimx     sdimy
#> cell_4  68.49701 13951.186
#> cell_5  73.80242 18085.385
#> cell_6 119.55364  3143.897
#> cell_7 105.25295  2191.132
#> cell_8  96.64224  4806.681
#> cell_9 110.59598 14124.008
```

We store the total number of detect genes and total counts as metadata:

``` r
metadata <- data.frame(num_genes = colSums(expr > 0),
                       total_count = colSums(expr))
```

## Spatial clustering

First, create a *BanksyObject* with the expression matrix and cell
locations.

``` r
bank <- BanksyObject(own.expr = expr,
                     cell.locs = locs,
                     meta.data = metadata)
```

Next, normalize the expression matrix, compute the neighbour matrix, and
scale the resulting gene-cell and neighbour feature-cell matrix.

``` r
bank <- NormalizeBanksy(bank, normFactor = 100)
bank <- ComputeBanksy(bank)
#> Computing neighbors...
#> Computing neighbor matrix...
#> Done
bank <- ScaleBanksy(bank)
```

At this point, the joint expression matrix (gene-cell matrix and
neighbour feature-cell matrix) can be extracted with *getBanksyMatrix*,
which returns the joint matrix and cell locations:

``` r
joint <- getBanksyMatrix(bank, lambda = 0.25)

joint$expr[1:5,1:5]
#>              cell_4      cell_5     cell_6     cell_7     cell_8
#> Slc1a2   0.41040297 -0.36165667 -0.4364615  0.8607711  0.3392922
#> Scn4b    0.85950683  2.02739671  0.2884143 -0.7755390 -0.4313188
#> Itpr1    0.77398900 -1.00281439  0.0173598 -0.7598327  0.1193772
#> Slc25a23 2.51032763 -1.15087647 -0.4501676 -0.3998602 -0.3800967
#> Slc1a3   0.07620877 -0.01765429 -0.5808326  3.5972767  0.2491144
joint$locs[1:5,]
#>            sdimx     sdimy
#> cell_4  68.49701 13951.186
#> cell_5  73.80242 18085.385
#> cell_6 119.55364  3143.897
#> cell_7 105.25295  2191.132
#> cell_8  96.64224  4806.681
```

Next, we obtain cluster assignments for the following parameters:

  - `lambda`. A mixing parameter from 0 to 1 which determines how much
    spatial information is incorporated.  
  - `resolution`. Leiden clustering resolution.  
  - `kneighbours`. Number of k neighbours to use for constructing sNN.

<!-- end list -->

``` r
set.seed(1234)
bank <- ClusterBanksy(bank, lambda = 0.5,
                            resolution = 0.8,
                            kneighbours = 40)
#> Iteration 1 out of 1
#> Finished clustering for Lambda=0.5, Resolution=0.8, K Neighbours=40
```

## Visualization

UMAP visualization:

``` r
plotUMAP(bank, by = 'res0.8_lam0.5_k40',  reduction = 'umap_0.5', pt.size = 0.02, main.size = 10)
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

Spatial plot:

``` r
plotSpatialDims(bank, by = 'res0.8_lam0.5_k40', pt.size = 0.8, main.size = 10)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

Heatmap:

``` r
plotHeatmap(bank, assay = 'banksy',
            lambda = 0.5,
            annotate = TRUE,
            annotate.by = 'res0.8_lam0.5_k40',
            barplot.by = c('num_genes', 'total_count'))
#> `use_raster` is automatically set to TRUE for a matrix with more than
#> 2000 columns You can control `use_raster` argument by explicitly
#> setting TRUE/FALSE to it.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
```

<img src="man/figures/README-unnamed-chunk-13-1.png" width="100%" />

## Subsetting

The *BanksyObject* can be subset by cells, features or metadata columns.

``` r
bankSubset1 <- SubsetBanksy(bank,
                            dims = (sdimx > 1000 & sdimx < 4000) & 
                                   (sdimy > 2000 & sdimy < 3000),
                            cells = c('cell_364', 'cell_367', 'cell_384','cell_389'),
                            features = sample(rownames(own.expr(bank)), 10))
head(bankSubset1)
#> own expression:
#>           cell_364   cell_367   cell_384   cell_389
#> Cacna1d -0.4665858 -0.4665858 -0.4665858 -0.4665858
#> Smad3   -0.6349799 -0.6349799 -0.6349799  1.2142839
#> Tsc22d1 -0.6742251 -0.6742251  0.6137452 -0.6742251
#> Traf1   -0.3479753 -0.3479753 -0.3479753 -0.3479753
#> Adamts2 -0.4211489 -0.4211489 -0.4211489 -0.4211489
#> 
#> neighbour expression:
#>               cell_364   cell_367   cell_384   cell_389
#> Cacna1d.nbr -0.4170890 -0.3539951 -0.6708647 -0.9097485
#> Smad3.nbr    0.8297271  0.5464924  0.1457932  0.4353906
#> Tsc22d1.nbr -0.3273127 -0.3267730 -0.2396807  0.4324479
#> Traf1.nbr   -0.3782101 -0.4923421  0.4701541  1.1875814
#> Adamts2.nbr -0.4919150 -0.6060863 -0.2951297 -0.2370119
#> 
#> cell locations:
#>             sdimx    sdimy
#> cell_364 1019.050 2660.030
#> cell_367 1040.727 2585.758
#> cell_384 1082.720 2383.199
#> cell_389 1084.240 2254.025
#> 
#> metadata:
#>           cell_ID num_genes total_count res0.8_lam0.5_k40
#> cell_364 cell_364        15          25                 3
#> cell_367 cell_367        19          35                 2
#> cell_384 cell_384        35          99                 6
#> cell_389 cell_389        35          60                 6
```

``` r
bankSubset2 <- SubsetBanksy(bank,
                            metadata = res0.8_lam0.5_k40 %in% c(3,6))
plotSpatialDims(bankSubset2, by = 'res0.8_lam0.5_k40', 
                pt.size = 0.5, main.size = 10)
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

# Session information

``` r
sessionInfo()
#> R version 4.0.3 (2020-10-10)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 19043)
#> 
#> Matrix products: default
#> 
#> locale:
#> [1] LC_COLLATE=English_Singapore.1252  LC_CTYPE=English_Singapore.1252   
#> [3] LC_MONETARY=English_Singapore.1252 LC_NUMERIC=C                      
#> [5] LC_TIME=English_Singapore.1252    
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] Banksy_0.99.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] maps_3.3.0               assertthat_0.2.1         highr_0.9               
#>  [4] stats4_4.0.3             yaml_2.2.1               pillar_1.6.1            
#>  [7] lattice_0.20-41          glue_1.4.2               RcppEigen_0.3.3.9.1     
#> [10] digest_0.6.27            RColorBrewer_1.1-2       colorspace_2.0-0        
#> [13] htmltools_0.5.1.1        Matrix_1.3-4             plyr_1.8.6              
#> [16] pkgconfig_2.0.3          GetoptLong_1.0.5         magick_2.7.2            
#> [19] purrr_0.3.4              scales_1.1.1             sccore_0.1.3            
#> [22] leidenAlg_0.1.1          ggalluvial_0.12.3        RSpectra_0.16-0         
#> [25] collapse_1.5.3           tibble_3.1.2             farver_2.0.3            
#> [28] generics_0.1.0           dbscan_1.1-8             IRanges_2.24.1          
#> [31] grr_0.9.5                ggplot2_3.3.5            ellipsis_0.3.2          
#> [34] BiocGenerics_0.36.1      magrittr_2.0.1           crayon_1.4.1            
#> [37] evaluate_0.14            fansi_0.5.0              RcppArmadillo_0.10.5.0.0
#> [40] Cairo_1.5-12.2           tools_4.0.3              data.table_1.13.6       
#> [43] GlobalOptions_0.1.2      lifecycle_1.0.0          matrixStats_0.59.0      
#> [46] ComplexHeatmap_2.6.2     stringr_1.4.0            S4Vectors_0.28.1        
#> [49] munsell_0.5.0            cluster_2.1.0            irlba_2.3.3             
#> [52] compiler_4.0.3           rlang_0.4.11             grid_4.0.3              
#> [55] Matrix.utils_0.9.8       dichromat_2.0-0          rjson_0.2.20            
#> [58] RcppAnnoy_0.0.18         circlize_0.4.13          igraph_1.2.6            
#> [61] labeling_0.4.2           rmarkdown_2.8            gtable_0.3.0            
#> [64] codetools_0.2-16         DBI_1.1.1                R6_2.5.0                
#> [67] knitr_1.33               dplyr_1.0.6              uwot_0.1.10             
#> [70] utf8_1.2.1               zeallot_0.1.0            clue_0.3-59             
#> [73] pals_1.7                 shape_1.4.6              stringi_1.5.3           
#> [76] parallel_4.0.3           Rcpp_1.0.5               vctrs_0.3.8             
#> [79] mapproj_1.2.7            png_0.1-7                tidyselect_1.1.1        
#> [82] xfun_0.23
```

## To-dos

  - Pass devtools::check: S3/4 documentation
  - Normalization methods for multiple datasets
  - Clustering methods
  - Optimization for handling large matrices (Rcpp)
  - Subsetting and Connecting clusters for multiple datasets
  - Unit tests
