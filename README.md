
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

  - improve cell-type assignment in noisy data
  - distinguish subtly different cell-types stratified by
    microenvironment
  - identify spatial zones sharing the same microenvironment

BANKSY is applicable to a wide array of spatial technologies (e.g. 10x
Visium, Slide-seq, MERFISH) and scales well to large datasets. For more
details, check out:

  - the
    [preprint](https://www.biorxiv.org/content/10.1101/2022.04.14.488259v1),
  - a
    [tweetorial](https://twitter.com/vipul1891/status/1515323372535644166?s=20&t=Bc6rz8VeWWptF67FejGYfQ)
    on BANKSY,
  - and a [Python version](https://github.com/prabhakarlab/Banksy_py) of
    this package.

## Installation

The *Banksy* package can be installed via `remotes`:

``` r
remotes::install_github("prabhakarlab/Banksy", dependencies = TRUE)
```

Installation should take less than three minutes.

**Known installation issues**

1.  Installation of `leidenAlg` has non-zero exit status

<!-- end list -->

  - Refer to the [package
    website](https://github.com/kharchenkolab/leidenAlg#installation)
    for *leidenAlg* installation details. Otherwise, users may also
    install a separate branch of *Banksy* with

<!-- end list -->

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

  - [Working with Banksy
    objects](https://prabhakarlab.github.io/Banksy/articles/banksy-object.html):
    Introduction to the *BanksyObject* class which serves as a container
    for *Banksy*.

  - [Mouse hippocampus VeraFISH
    dataset](https://prabhakarlab.github.io/Banksy/articles/hippocampus-analysis.html):
    Illustrates a grid search of parameters which best cluster cells.

  - [Human dorsolateral prefrontal cortex 10x Visium
    dataset](https://prabhakarlab.github.io/Banksy/articles/dlpfc-analysis.html):
    Illustrates analysis of multiple spatial transcriptomic datasets.

  - [Mouse hypothalamus MERFISH
    dataset](https://prabhakarlab.github.io/Banksy/articles/hypothalamus-analysis.html):
    Illustrates visualization functionality with a dataset with 3
    spatial dimensions.

  - [Interoperability with
    SingleCellExperiment](https://prabhakarlab.github.io/Banksy/articles/single-cell-exp.html):
    Illustrates BANKSY interoperability with Bioconductor
    [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
    framework for interfacing with packages like
    [scran](https://bioconductor.org/packages/release/bioc/html/scran.html)
    or
    [scater](https://bioconductor.org/packages/release/bioc/html/scater.html).

  - [Figure 4 data
    analysis](https://prabhakarlab.github.io/Banksy/articles/Fig4-vignette.html):
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

We store the total counts for each cell and the number of expressed
genes as metadata `data.frame`, which can optionally be supplied:

``` r
total_count <- colSums(expr)
num_genes <- colSums(expr > 0)
meta <- data.frame(total_count = total_count, num_genes = num_genes)
```

Next, create a *BanksyObject* with the expression matrix and cell
locations (metadata is optional).

``` r
bank <- BanksyObject(own.expr = expr,
                     cell.locs = locs,
                     meta.data = meta)
```

Apply basic QC by keeping only cells with total counts within the 5th
and 98th percentile:

``` r
bank <- SubsetBanksy(bank, metadata = total_count > quantile(total_count, 0.05) &
                                      total_count < quantile(total_count, 0.98))
```

We first normalize the expression matrix, compute the neighbor matrix,
and scale the resulting matrices.

``` r
bank <- NormalizeBanksy(bank, normFactor = 100)
bank <- ComputeBanksy(bank, k_geom = 10, spatialMode = 'kNN_r')
#> Computing neighbors...
#> Computing neighbor matrix...
#> Done
bank <- ScaleBanksy(bank)
```

Run PCA on the BANKSY matrix for `lambda = 0` (no spatial information)
and `lambda = 0.3`.

``` r
bank <- RunPCA(bank, lambda = c(0, 0.3), npcs = 30)
#> Running PCA for lambda=0
#> Running PCA for lambda=0.3
```

Next, we obtain cluster assignments using graph-based clustering with
the Leiden algorithm on the first 20 PCs. Specify the following
parameters:

  - `resolution`. Leiden clustering resolution.  
  - `k.neighbours`. Number of k neighbours to use for constructing sNN.

<!-- end list -->

``` r
set.seed(42)
bank <- ClusterBanksy(bank, lambda = c(0, 0.3), pca = TRUE, npcs = 20,
                      method = 'leiden', k.neighbors = 50, resolution = 1.2)
#> Iteration 1 out of 2
#> Iteration 2 out of 2
```

Different clustering runs can be harmonised with `ConnectClusters`:

``` r
bank <- ConnectClusters(bank, map.to = clust.names(bank)[1])
```

Visualise the clustering output for non-spatial clustering (`lambda=0`)
and BANKSY clustering (`lambda = 0.3`).

``` r
features <- clust.names(bank)
feature.types <- rep('discrete', 2)
main <- c('Non-spatial', 'BANKSY')
plotSpatialFeatures(bank, by = features, type = feature.types, main = main, 
                    pt.size = 1.5, main.size = 15, nrow = 1, ncol = 2)
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

For clarity, we can visualise each of the clusters separately with `wrap
= TRUE`:

``` r
plotSpatialFeatures(bank, by = features, type = feature.types, main = main, 
                    pt.size = 0.5, main.size = 15, nrow = 1, ncol = 2, 
                    wrap = TRUE)
```

<img src="man/figures/README-unnamed-chunk-16-1.png" width="100%" />

<details>

<summary>Runtime for analysis</summary>

    #> Time difference of 34.56732 secs

</details>

<details>

<summary>Session information</summary>

``` r
sessionInfo()
#> R version 4.0.2 (2020-06-22)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 18.04.6 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
#> LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] Banksy_0.1.3
#> 
#> loaded via a namespace (and not attached):
#>  [1] bitops_1.0-6                matrixStats_0.59.0         
#>  [3] RColorBrewer_1.1-2          GenomeInfoDb_1.26.7        
#>  [5] tools_4.0.2                 utf8_1.2.1                 
#>  [7] R6_2.5.0                    irlba_2.3.3                
#>  [9] uwot_0.1.10                 DBI_1.1.0                  
#> [11] BiocGenerics_0.36.1         colorspace_2.0-2           
#> [13] GetoptLong_1.0.5            tidyselect_1.1.1           
#> [15] gridExtra_2.3               compiler_4.0.2             
#> [17] Biobase_2.50.0              Cairo_1.5-12.2             
#> [19] DelayedArray_0.16.3         labeling_0.4.2             
#> [21] scales_1.1.1                stringr_1.4.0              
#> [23] digest_0.6.27               dbscan_1.1-8               
#> [25] rmarkdown_2.3               XVector_0.30.0             
#> [27] dichromat_2.0-0             pkgconfig_2.0.3            
#> [29] htmltools_0.5.0             MatrixGenerics_1.2.1       
#> [31] highr_0.8                   maps_3.3.0                 
#> [33] rlang_0.4.11                GlobalOptions_0.1.2        
#> [35] pals_1.7                    shape_1.4.6                
#> [37] generics_0.1.0              farver_2.1.0               
#> [39] mclust_5.4.7                dplyr_1.0.7                
#> [41] RCurl_1.98-1.2              magrittr_2.0.1             
#> [43] GenomeInfoDbData_1.2.4      Matrix_1.3-4               
#> [45] Rcpp_1.0.7                  munsell_0.5.0              
#> [47] S4Vectors_0.28.1            fansi_0.5.0                
#> [49] lifecycle_1.0.0             stringi_1.4.6              
#> [51] leidenAlg_0.1.1             yaml_2.2.1                 
#> [53] ggalluvial_0.12.3           SummarizedExperiment_1.20.0
#> [55] zlibbioc_1.36.0             plyr_1.8.6                 
#> [57] grid_4.0.2                  blob_1.2.1                 
#> [59] parallel_4.0.2              crayon_1.4.1               
#> [61] lattice_0.20-41             sccore_0.1.3               
#> [63] mapproj_1.2.7               circlize_0.4.13            
#> [65] knitr_1.36                  ComplexHeatmap_2.6.2       
#> [67] pillar_1.6.1                igraph_1.2.6               
#> [69] GenomicRanges_1.42.0        rjson_0.2.20               
#> [71] stats4_4.0.2                glue_1.4.2                 
#> [73] evaluate_0.14               data.table_1.14.0          
#> [75] png_0.1-7                   vctrs_0.3.8                
#> [77] gtable_0.3.0                grr_0.9.5                  
#> [79] purrr_0.3.4                 clue_0.3-59                
#> [81] assertthat_0.2.1            ggplot2_3.3.5              
#> [83] xfun_0.28                   tibble_3.1.3               
#> [85] RcppHungarian_0.1           Matrix.utils_0.9.8         
#> [87] IRanges_2.24.1              cluster_2.1.0              
#> [89] ellipsis_0.3.2
```

</details>

## Contributing

Bug reports, questions, request for enhancements or other contributions
can be raised at the [issue
page](https://github.com/prabhakarlab/Banksy/issues).
