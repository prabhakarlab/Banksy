
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/jleechung/Banksy/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jleechung/Banksy/actions/workflows/R-CMD-check.yaml)
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

-   Refer to the [package
    website](https://github.com/kharchenkolab/leidenAlg#installation)
    for *leidenAlg* installation details. Otherwise, users may also
    install a separate branch of *Banksy* with

``` r
remotes::install_github("prabhakarlab/Banksy@feat-igraph-leiden")
```

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

-   [Figure 4 data
    analysis](https://prabhakarlab.github.io/Banksy/articles/Fig4_vignette.html):
    Shows how the results shown in Fig. 4 of the paper were generated.

*Banksy* is also interoperable with
[Seurat](https://satijalab.org/seurat/) via *SeuratWrappers*.
Documentation on how to run BANKSY on Seurat objects can be found
[here](https://github.com/satijalab/seurat-wrappers/blob/master/docs/banksy.md).

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

-   `resolution`. Leiden clustering resolution.

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

    #> Time difference of 46.51171 secs

</details>
<details>
<summary>
Session information
</summary>

``` r
sessionInfo()
#> R version 4.1.2 (2021-11-01)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 19045)
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
#> [1] Banksy_0.1.4
#> 
#> loaded via a namespace (and not attached):
#>  [1] bitops_1.0-7                matrixStats_0.61.0         
#>  [3] progress_1.2.2              doParallel_1.0.17          
#>  [5] RColorBrewer_1.1-3          GenomeInfoDb_1.30.1        
#>  [7] tools_4.1.2                 utf8_1.2.2                 
#>  [9] R6_2.5.1                    irlba_2.3.5                
#> [11] uwot_0.1.11                 DBI_1.1.2                  
#> [13] BiocGenerics_0.40.0         colorspace_2.0-2           
#> [15] GetoptLong_1.0.5            prettyunits_1.1.1          
#> [17] tidyselect_1.1.2            gridExtra_2.3              
#> [19] compiler_4.1.2              cli_3.1.0                  
#> [21] Biobase_2.54.0              DelayedArray_0.20.0        
#> [23] labeling_0.4.2              scales_1.2.0               
#> [25] stringr_1.4.0               digest_0.6.29              
#> [27] dbscan_1.1-10               rmarkdown_2.13             
#> [29] XVector_0.34.0              dichromat_2.0-0.1          
#> [31] pkgconfig_2.0.3             htmltools_0.5.2            
#> [33] MatrixGenerics_1.6.0        highr_0.9                  
#> [35] fastmap_1.1.0               maps_3.4.0                 
#> [37] rlang_1.0.2                 GlobalOptions_0.1.2        
#> [39] pals_1.7                    rstudioapi_0.13            
#> [41] farver_2.1.0                shape_1.4.6                
#> [43] generics_0.1.2              mclust_5.4.9               
#> [45] dplyr_1.0.7                 RCurl_1.98-1.6             
#> [47] magrittr_2.0.1              GenomeInfoDbData_1.2.7     
#> [49] Matrix_1.3-4                Rcpp_1.0.9                 
#> [51] munsell_0.5.0               S4Vectors_0.32.3           
#> [53] fansi_0.5.0                 lifecycle_1.0.1            
#> [55] stringi_1.7.6               leidenAlg_1.0.3            
#> [57] yaml_2.2.1                  ggalluvial_0.12.3          
#> [59] SummarizedExperiment_1.24.0 zlibbioc_1.40.0            
#> [61] plyr_1.8.6                  grid_4.1.2                 
#> [63] parallel_4.1.2              crayon_1.5.1               
#> [65] lattice_0.20-45             sccore_1.0.1               
#> [67] hms_1.1.1                   mapproj_1.2.8              
#> [69] circlize_0.4.15             knitr_1.37                 
#> [71] ComplexHeatmap_2.10.0       pillar_1.7.0               
#> [73] igraph_1.3.4                GenomicRanges_1.46.1       
#> [75] rjson_0.2.21                codetools_0.2-18           
#> [77] stats4_4.1.2                glue_1.6.0                 
#> [79] evaluate_0.15               data.table_1.14.2          
#> [81] png_0.1-7                   vctrs_0.3.8                
#> [83] foreach_1.5.2               gtable_0.3.0               
#> [85] grr_0.9.5                   purrr_0.3.4                
#> [87] clue_0.3-60                 assertthat_0.2.1           
#> [89] ggplot2_3.3.6               xfun_0.29                  
#> [91] tibble_3.1.6                RcppHungarian_0.2          
#> [93] iterators_1.0.14            Matrix.utils_0.9.8         
#> [95] IRanges_2.28.0              cluster_2.1.2              
#> [97] ellipsis_0.3.2
```

</details>

## Contributing

Bug reports, questions, request for enhancements or other contributions
can be raised at the [issue
page](https://github.com/prabhakarlab/Banksy/issues).
