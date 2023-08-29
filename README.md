
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/jleechung/Banksy/actions/workflows/check-standard.yml/badge.svg)](https://github.com/jleechung/Banksy/actions/workflows/check-standard.yml)
[![codecov](https://codecov.io/gh/jleechung/Banksy/branch/main/graph/badge.svg?token=OZZK4EDVH9)](https://codecov.io/gh/jleechung/Banksy)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/jleechung/Banksy/graphs/commit-activity)
<!-- badges: end -->

## Overview

BANKSY is a method for clustering spatial transcriptomic data by
augmenting the transcriptomic profile of each cell with an average of
the transcriptomes of its spatial neighbors. By incorporating
neighborhood information for clustering, BANKSY is able to

-   improve cell-type assignment in noisy data
-   distinguish subtly different cell-types stratified by
    microenvironment
-   identify spatial zones sharing the same microenvironment

BANKSY is applicable to a wide array of spatial technologies (e.g. 10x
Visium, Slide-seq, MERFISH) and scales well to large datasets. For more
details, check out:

-   the
    [preprint](https://www.biorxiv.org/content/10.1101/2022.04.14.488259v1),
-   a
    [tweetorial](https://twitter.com/vipul1891/status/1515323372535644166?s=20&t=Bc6rz8VeWWptF67FejGYfQ)
    on BANKSY,
-   and a [Python version](https://github.com/prabhakarlab/Banksy_py) of
    this package.

## Installation

The *Banksy* package can be installed via `remotes`:

``` r
remotes::install_github("prabhakarlab/Banksy", dependencies = TRUE)
```

Installation should take less than three minutes.

**Known installation issues**

1.  Installation of `leidenAlg` has non-zero exit status

-   Refer to the [`leidenAlg` package
    website](https://github.com/kharchenkolab/leidenAlg#installation)
    for installation details.

## Documentation

Detailed description of *Banksy* functionality and example analyses are
available at the [package
webpage](https://prabhakarlab.github.io/Banksy/).

*Banksy* comes installed with
[documentation](https://prabhakarlab.github.io/Banksy/reference/index.html)
of main functions and their usage, along with several vignettes which
detail different use cases:

-   [Working with Banksy
    objects](https://prabhakarlab.github.io/Banksy/articles/banksy-object.html):
    Introduction to the *BanksyObject* class which serves as a container
    for *Banksy*.

-   [Mouse hippocampus VeraFISH
    dataset](https://prabhakarlab.github.io/Banksy/articles/hippocampus-analysis.html):
    Illustrates a grid search of parameters which best cluster cells.

-   [Human dorsolateral prefrontal cortex 10x Visium
    dataset](https://prabhakarlab.github.io/Banksy/articles/dlpfc-analysis.html):
    Illustrates analysis of multiple spatial transcriptomic datasets.

-   [Mouse hypothalamus MERFISH
    dataset](https://prabhakarlab.github.io/Banksy/articles/hypothalamus-analysis.html):
    Illustrates visualization functionality with a dataset with 3
    spatial dimensions.

-   [Interoperability with
    SingleCellExperiment](https://prabhakarlab.github.io/Banksy/articles/single-cell-exp.html):
    Illustrates BANKSY interoperability with Bioconductor
    [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
    framework for interfacing with packages like
    [scran](https://bioconductor.org/packages/release/bioc/html/scran.html)
    or
    [scater](https://bioconductor.org/packages/release/bioc/html/scater.html).

*Banksy* is also interoperable with
[Seurat](https://satijalab.org/seurat/) via *SeuratWrappers*.
Documentation on how to run BANKSY on Seurat objects can be found
[here](https://github.com/jleechung/seurat-wrappers/blob/feat-aft/docs/banksy.md).
For installation of *SeuratWrappers* with BANKSY version `>= 0.1.4`, run

``` r
remotes::install_github('jleechung/seurat-wrappers@feat-aft')
```

## Quick start

*Banksy* takes as input an expression matrix and cell centroids. Example
datasets are provided with the package:

``` r
library(Banksy)

data(hippocampus)
expr <- hippocampus$expression
locs <- hippocampus$locations
```

The gene expression matrix for cells should be a `matrix`:

``` r
class(expr)
#> [1] "matrix" "array"
head(expr[,1:5])
#>         cell_1276 cell_8890 cell_691 cell_396 cell_9818
#> Sparcl1        45         0       11       22         0
#> Slc1a2         17         0        6        5         0
#> Map            10         0       12       16         0
#> Sqstm1         26         0        0        2         0
#> Atp1a2          0         0        4        3         0
#> Tnc             0         0        0        0         0
```

while cell locations should be supplied as a `data.frame`:

``` r
class(locs)
#> [1] "data.frame"
head(locs)
#>                 sdimx    sdimy
#> cell_1276  -13372.899 15776.37
#> cell_8890    8941.101 15866.37
#> cell_691   -14882.899 15896.37
#> cell_396   -15492.899 15835.37
#> cell_9818   11308.101 15846.37
#> cell_11310  14894.101 15810.37
```

Next, create a *BanksyObject* with the expression matrix and cell
locations.

``` r
bank <- BanksyObject(own.expr = expr,
                     cell.locs = locs)

head(bank)
#> own expression:
#>         cell_1276 cell_8890 cell_691 cell_396 cell_9818
#> Sparcl1        45         0       11       22         0
#> Slc1a2         17         0        6        5         0
#> Map            10         0       12       16         0
#> Sqstm1         26         0        0        2         0
#> Atp1a2          0         0        4        3         0
#> 
#> neighbour expression:
#> NULL
#> 
#> cell locations:
#>                sdimx    sdimy
#> cell_1276 -13372.899 15776.37
#> cell_8890   8941.101 15866.37
#> cell_691  -14882.899 15896.37
#> cell_396  -15492.899 15835.37
#> cell_9818  11308.101 15846.37
#> 
#> metadata:
#>             cell_ID nCount NODG
#> cell_1276 cell_1276    266   51
#> cell_8890 cell_8890     13    3
#> cell_691   cell_691    132   36
#> cell_396   cell_396     95   27
#> cell_9818 cell_9818     10    5
```

Apply basic QC by keeping only cells with total counts within the 5th
and 98th percentile:

``` r
bank <- SubsetBanksy(bank, metadata = nCount > quantile(nCount, 0.05) &
                                      nCount < quantile(nCount, 0.98))
```

We first normalize the expression matrix, compute the neighbor matrix,
and scale the resulting matrices.

``` r
bank <- NormalizeBanksy(bank)
bank <- ComputeBanksy(bank, k_geom = c(15, 30))
#> Computing neighbors...
#> Spatial mode is kNN median
#> Parameters: k_geom = 15
#> Done
#> Computing neighbors...
#> Spatial mode is kNN median
#> Parameters: k_geom = 30
#> Done
#> Computing harmonic m = 0
#> Using 15 neighbors
#> Done
#> Computing harmonic m = 1
#> Using 30 neighbors
#> Centering
#> Done
bank <- ScaleBanksy(bank)
```

Run PCA on the BANKSY matrix for `lambda=0` (no spatial information) and
`lambda=0.2`.

``` r
bank <- RunBanksyPCA(bank, lambda = c(0, 0.2))
#> Running PCA for M=1 lambda=0
#> BANKSY matrix with own.expr, F0, F1
#> Squared lambdas: 1, 0, 0
#> Running PCA for M=1 lambda=0.2
#> BANKSY matrix with own.expr, F0, F1
#> Squared lambdas: 0.8, 0.1333, 0.0667
```

Next, we obtain cluster assignments using graph-based clustering with
the Leiden algorithm on the first 20 PCs.

``` r
set.seed(42)
bank <- ClusterBanksy(bank, lambda = c(0, 0.2), pca = TRUE, npcs = 20,
                      method = 'leiden', resolution = 1.2)
```

Different clustering runs can be harmonised with `ConnectClusters`:

``` r
bank <- ConnectClusters(bank, map.to = clust.names(bank)[1])
```

Visualise the clustering output for non-spatial clustering (`lambda=0`)
and BANKSY clustering (`lambda=0.2`).

``` r
features <- clust.names(bank)
feature.types <- rep('discrete', 2)
main <- c('Non-spatial', 'BANKSY')
plotSpatialFeatures(bank, by = features, type = feature.types, main = main, 
                    pt.size = 1.5, main.size = 15, nrow = 1, ncol = 2)
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

For clarity, we can visualise each of the clusters separately with
`wrap = TRUE`:

``` r
plotSpatialFeatures(bank, by = features, type = feature.types, main = main, 
                    pt.size = 1, main.size = 15, nrow = 1, ncol = 2, 
                    wrap = TRUE)
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

<details>
<summary>
Runtime for analysis
</summary>

<<<<<<< HEAD
    #> Time difference of 1.915979 mins
=======
    #> Time difference of 38.77595 secs
>>>>>>> main

</details>
<details>
<summary>
Session information
</summary>

``` r
sessionInfo()
#> R version 4.2.1 (2022-06-23 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 19043)
#> 
#> Matrix products: default
#> 
#> locale:
#> [1] LC_COLLATE=English_Singapore.utf8  LC_CTYPE=English_Singapore.utf8   
#> [3] LC_MONETARY=English_Singapore.utf8 LC_NUMERIC=C                      
#> [5] LC_TIME=English_Singapore.utf8    
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] Banksy_0.1.5
#> 
#> loaded via a namespace (and not attached):
#>  [1] bitops_1.0-7                matrixStats_0.62.0         
#>  [3] doParallel_1.0.17           RColorBrewer_1.1-3         
#>  [5] progress_1.2.2              GenomeInfoDb_1.34.4        
#>  [7] tools_4.2.1                 utf8_1.2.2                 
#>  [9] R6_2.5.1                    irlba_2.3.5.1              
#> [11] uwot_0.1.14                 DBI_1.1.3                  
#> [13] BiocGenerics_0.44.0         colorspace_2.0-3           
#> [15] GetoptLong_1.0.5            withr_2.5.0                
#> [17] tidyselect_1.2.0            gridExtra_2.3              
#> [19] prettyunits_1.1.1           compiler_4.2.1             
#> [21] cli_3.4.1                   Biobase_2.58.0             
#> [23] DelayedArray_0.24.0         labeling_0.4.2             
#> [25] scales_1.2.1                stringr_1.5.0              
#> [27] digest_0.6.30               dbscan_1.1-11              
#> [29] rmarkdown_2.19              XVector_0.38.0             
#> [31] dichromat_2.0-0.1           pkgconfig_2.0.3            
#> [33] htmltools_0.5.4             MatrixGenerics_1.10.0      
#> [35] highr_0.10                  fastmap_1.1.0              
#> [37] maps_3.4.1                  rlang_1.1.1                
#> [39] GlobalOptions_0.1.2         pals_1.7                   
#> [41] rstudioapi_0.14             farver_2.1.1               
#> [43] shape_1.4.6                 generics_0.1.3             
#> [45] mclust_6.0.0                dplyr_1.0.10               
#> [47] RCurl_1.98-1.9              magrittr_2.0.3             
#> [49] GenomeInfoDbData_1.2.9      Matrix_1.5-3               
#> [51] Rcpp_1.0.9                  munsell_0.5.0              
#> [53] S4Vectors_0.36.0            fansi_1.0.3                
#> [55] lifecycle_1.0.3             stringi_1.7.8              
#> [57] leidenAlg_1.1.0             yaml_2.3.6                 
#> [59] ggalluvial_0.12.3           SummarizedExperiment_1.28.0
#> [61] zlibbioc_1.44.0             plyr_1.8.8                 
#> [63] grid_4.2.1                  parallel_4.2.1             
#> [65] crayon_1.5.2                lattice_0.20-45            
#> [67] sccore_1.0.2                mapproj_1.2.9              
#> [69] circlize_0.4.15             hms_1.1.2                  
#> [71] knitr_1.41                  ComplexHeatmap_2.14.0      
#> [73] pillar_1.8.1                igraph_1.3.5               
#> [75] GenomicRanges_1.50.1        rjson_0.2.21               
#> [77] codetools_0.2-18            stats4_4.2.1               
#> [79] glue_1.6.2                  evaluate_0.19              
#> [81] data.table_1.14.6           png_0.1-7                  
#> [83] vctrs_0.5.1                 foreach_1.5.2              
#> [85] gtable_0.3.1                clue_0.3-62                
#> [87] assertthat_0.2.1            ggplot2_3.4.0              
#> [89] xfun_0.36                   tibble_3.1.8               
#> [91] RcppHungarian_0.2           iterators_1.0.14           
#> [93] IRanges_2.32.0              cluster_2.1.4              
#> [95] ellipsis_0.3.2
```

</details>

## Contributing

Bug reports, questions, request for enhancements or other contributions
can be raised at the [issue
page](https://github.com/prabhakarlab/Banksy/issues).
