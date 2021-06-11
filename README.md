
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

## Installation

``` r
remotes::install_github("jleechung/Banksy")
```

Load the package:

``` r
library(Banksy)
```

## Quick Start

### Input data

Inputs consist of an expression matrix and cell locations. Sample data
provided with the
package:

``` r
expr <- readRDS(system.file('/extdata/expression.rds', package = 'Banksy'))
locs <- readRDS(system.file('/extdata/locations.rds', package = 'Banksy'))
```

The gene expression matrix for cells should either be a `data.frame` or
sparse matrix of class `dgCMatrix`:

``` r
head(expr[,1:5])
#>          cell_4 cell_5 cell_6 cell_7 cell_8
#> Slc1a2        4      2     13     26      6
#> Scn4b         3      6     17      0      1
#> Itpr1         3      0     15      2      3
#> Slc25a23      3      0      5      3      1
#> Slc1a3        1      1      0     31      2
#> Nfib          1      1      1      4      0
```

The locations of each cell as a `data.frame`:

``` r
head(locs)
#>            sdimx     sdimy
#> cell_4  68.49701 13951.186
#> cell_5  73.80242 18085.385
#> cell_6 119.55364  3143.897
#> cell_7 105.25295  2191.132
#> cell_8  96.64224  4806.681
#> cell_9 110.59598 14124.008
```

### Spatial clustering

First, create a *BanksyObject* with the expression matrix and cell
locations.

``` r
bank <- BanksyObject(own.expr = expr,
                     cell.locs = locs)
```

Compute the neighbour matrix, and scale and normalize matrix. The
options below are defaults and listed for clarity.

``` r
bank <- ComputeBanksy(bank,
                      normalizeColumns = TRUE,
                      normalizeColumnsTo = 100,
                      zScaleRows = TRUE,
                      zScaleBeforeAveraging = FALSE,
                      zScaleOwnAfterAveraging = TRUE,
                      zScaleNbrAfterAveraging = TRUE)
#> Performing normalization...
#> Computing Banksy matrices...
#> Spatial mode is kNN_r, k_geom = 10
#> Banksy matrix: 62.397 sec elapsed
```

At this point, the joint expression matrix (gene-cell matrix and
neighbour feature-cell matrix) can be extracted with *getBanksyMatrix*:

``` r
joint <- getBanksyMatrix(bank, lambda = 0.25)
```

Obtain clusters for selected parameters.

  - `lambda`. A mixing parameter from 0 to 1 which determines how much
    spatial information is incorporated.  
  - `resolution`. Leiden clustering resolution.  
  - `kneighbours`. Number of k neighbours to use for constructing sNN.

<!-- end list -->

``` r
bank <- ClusterBanksy(bank, lambda = 0.25,
                            resolution = 1.2,
                            kneighbours = 30)
#> 
#>  external python path provided and will be used
#> Iteration 1 out of 1
#> Consider to install these (optional) packages to run all possible Giotto commands:  RTriangle FactoMiner
#>  Giotto does not automatically install all these packages as they are not absolutely required and this reduces the number of dependencieshvg  was not found in the gene metadata information, all genes will be used
#> Finished clustering for Lambda=0.25, Resolution=1.2, K Neighbours=30
#> 29.666 sec elapsed
```

To explore the parameter space, one can provide sequences for the
parameters:

``` r
## Not run
bank <- ClusterBanksy(bank, 
                      lambda = seq(0, 0.75, by = 0.25),
                      resolution = seq(0.8, 1.2, by = 0.2),
                      kneighbours = seq(30, 40, by = 10))
```

### Visualization

UMAP visualization:

``` r
plotUMAP(bank, params = 'res1.2_lam0.25_k30', pt.size = 0.02)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

Spatial plot:

``` r
plotSpatialDims(bank, params = 'res1.2_lam0.25_k30', pt.size = 0.5)
```

<img src="man/figures/README-unnamed-chunk-13-1.png" width="100%" />

### Subsetting

The *BanksyObject* can be subset by cells, features or metadata columns.

``` r
bankSubset1 <- SubsetBanksy(bank,
                              dimx = sdimx > 1000 & sdimx < 4000,
                              dimy = sdimy > 2000 & sdimy < 3000,
                              cells = c(cell_364, cell_367, cell_384,cell_389),
                              features = sample(rownames(own.expr(bank)), 10))
#> Before filtering: 6725 cells and 140 genes.
#> After filtering: 4 cells and 10 genes.
#> Filtered 6721 cells and 130 genes.
head(bankSubset1)
#> own expression:
#>         cell_364 cell_367 cell_384 cell_389
#> Slc1a2         6        5       33        6
#> Grm4           0        0        0        0
#> Flt1           0        0        1        0
#> Cacna1e        0        0        0        1
#> Adamts2        0        0        0        0
#> 
#> neighbour expression:
#>               cell_364   cell_367   cell_384   cell_389
#> Slc1a2.nbr  11.3054557 14.9663042 12.9835741 11.3417701
#> Grm4.nbr     0.9118782  0.9519174  0.7268138  0.4679934
#> Flt1.nbr     0.3687012  1.0479533  1.3290719  0.7870465
#> Cacna1e.nbr  0.9909432  0.9872081  0.9295227  0.3704704
#> Adamts2.nbr  0.1261979  0.1027164  0.1666704  0.1786234
#> 
#> own normalized scaled expression:
#>           cell_364   cell_367   cell_384   cell_389
#> Slc1a2   1.1603457  0.2510181  2.0340135 -0.1501558
#> Grm4    -0.6799992 -0.6799992 -0.6799992 -0.6799992
#> Flt1    -0.3653916 -0.3653916  0.0225537 -0.3653916
#> Cacna1e -0.6094535 -0.6094535 -0.6094535  1.0911288
#> Adamts2 -0.4211489 -0.4211489 -0.4211489 -0.4211489
#> 
#> neighbour normalized scaled expression:
#>                cell_364    cell_367    cell_384    cell_389
#> Slc1a2.nbr  -0.04780222  0.55959368  0.23062565 -0.04177705
#> Grm4.nbr     0.30465882  0.37637906 -0.02683804 -0.49045043
#> Flt1.nbr    -0.56636798  0.08228811  0.35074399 -0.16686647
#> Cacna1e.nbr  0.75835950  0.75109842  0.63895757 -0.44784377
#> Adamts2.nbr -0.49191500 -0.60608629 -0.29512972 -0.23701190
#> 
#> cell locations:
#>             sdimx    sdimy
#> cell_364 1019.050 2660.030
#> cell_367 1040.727 2585.758
#> cell_384 1082.720 2383.199
#> cell_389 1084.240 2254.025
#> 
#> metadata:
#>      cell_ID res1.2_lam0.25_k30
#> 352 cell_364                  7
#> 355 cell_367                  3
#> 372 cell_384                  4
#> 377 cell_389                  4
```

``` r
bankSubset2 <- SubsetBanksy(bank,
                            metadata = res1.2_lam0.25_k30 %in% 1:4)
#> Before filtering: 6725 cells and 140 genes.
#> After filtering: 3877 cells and 140 genes.
#> Filtered 2848 cells and 0 genes.
plotSpatialDims(bankSubset2, params = 'res1.2_lam0.25_k30', pt.size = 0.5)
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

## Session information

``` r
sessionInfo()
#> R version 3.6.0 (2019-04-26)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 18.04.2 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
#> 
#> locale:
#> [1] C
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] Banksy_0.99.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_1.0.5         RSpectra_0.16-0    pracma_2.3.3       RColorBrewer_1.1-2
#>  [5] plyr_1.8.5         pillar_1.4.1       compiler_3.6.0     tools_3.6.0       
#>  [9] uwot_0.1.10        digest_0.6.19      jsonlite_1.6       evaluate_0.14     
#> [13] lifecycle_0.1.0    tibble_2.1.3       gtable_0.3.0       lattice_0.20-38   
#> [17] pkgconfig_2.0.2    rlang_0.4.11       igraph_1.2.4.2     Matrix_1.2-17     
#> [21] parallel_3.6.0     ggrepel_0.8.1      yaml_2.2.0         xfun_0.12         
#> [25] Giotto_1.0.2       stringr_1.4.0      dplyr_0.8.4        knitr_1.27        
#> [29] rappdirs_0.3.1     tictoc_1.0.1       grid_3.6.0         tidyselect_0.2.5  
#> [33] cowplot_1.0.0      reticulate_1.14    glue_1.3.1         data.table_1.12.8 
#> [37] R6_2.4.0           RcppAnnoy_0.0.18   rmarkdown_2.1      irlba_2.3.3       
#> [41] farver_2.0.3       ggplot2_3.3.2      purrr_0.3.3        magrittr_1.5      
#> [45] codetools_0.2-16   scales_1.1.0       htmltools_0.5.1.1  assertthat_0.2.1  
#> [49] colorspace_1.4-1   labeling_0.3       stringi_1.4.5      munsell_0.5.0     
#> [53] crayon_1.3.4       dbscan_1.1-5
```
