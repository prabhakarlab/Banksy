
<!-- README.md is generated from README.Rmd. Please edit that file -->

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
remotes::install_github('prabhakarlab/Banksy')
```

## Quick start

Load *BANKSY*. We’ll also load *SpatialExperiment* and
*SummarizedExperiment* for containing and manipulating the data,
*scuttle* for normalization and quality control, and *scater*, *ggplot2*
and *cowplot* for visualisation.

``` r
library(Banksy)

library(SummarizedExperiment)
library(SpatialExperiment)
library(scuttle)

library(scater)
library(cowplot)
library(ggplot2)
```

Here, we’ll run *BANKSY* on mouse hippocampus data.

``` r
data(hippocampus)
gcm <- hippocampus$expression
locs <- as.matrix(hippocampus$locations)
```

Initialize a SpatialExperiment object and perform basic quality control
and normalization.

``` r
se <- SpatialExperiment(assay = list(counts=gcm), spatialCoords = locs)

# QC based on total counts
qcstats <- perCellQCMetrics(se)
thres <- quantile(qcstats$total, c(0.05,0.98))
keep <- (qcstats$total > thres[1]) & (qcstats$total < thres[2])
se <- se[, keep]

# Normalization to mean library size
se <- computeLibraryFactors(se)
aname <- 'normcounts'
assay(se, aname) <- normalizeCounts(se, log=FALSE)
```

Compute the neighborhood matrices for *BANKSY*. Setting
`compute_agf=TRUE` computes both the weighted neighborhood mean
($\mathcal{M}$) and the azimuthal Gabor filter ($\mathcal{G}$). The
number of spatial neighbors used to compute $\mathcal{M}$ and
$\mathcal{G}$ are `k_geom[1]=15` and `k_geom[2]=30` respectively. We run
*BANKSY* at `lambda=0` corresponding to non-spatial clustering, and
`lambda=0.2` corresponding to *BANKSY* for cell-typing.

``` r
lambda <- c(0, 0.2)
k_geom <- c(15, 30)

se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)
#> Computing neighbors...
#> Spatial mode is kNN_median
#> Parameters: k_geom=15
#> Done
#> Computing neighbors...
#> Spatial mode is kNN_median
#> Parameters: k_geom=30
#> Done
#> Computing harmonic m = 0
#> Using 15 neighbors
#> Done
#> Computing harmonic m = 1
#> Using 30 neighbors
#> Centering
#> Done
```

Next, run PCA on the BANKSY matrix and perform clustering. Setting
`use_agf=TRUE` uses both $\mathcal{M}$ and $\mathcal{G}$ to construct
the BANKSY matrix.

``` r
set.seed(1000)
se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = 1.2)
```

Different clustering runs can be relabeled to minimise their differences
with `connectClusters`:

``` r
se <- Banksy::connectClusters(se)
#> clust_M1_lam0.2_k50_res1.2 --> clust_M1_lam0_k50_res1.2
```

Visualise the clustering output for non-spatial clustering (`lambda=0`)
and BANKSY clustering (`lambda=0.2`).

``` r
cnames <- colnames(colData(se))
cnames <- cnames[grep('^clust', cnames)]
colData(se) <- cbind(colData(se), spatialCoords(se))

plot_nsp <- plotColData(
    se,
    x = 'sdimx',
    y = 'sdimy',
    point_size = 0.6,
    colour_by = cnames[1]
) 
plot_bank <- plotColData(
    se,
    x = 'sdimx',
    y = 'sdimy',
    point_size = 0.6,
    colour_by = cnames[2]
)


plot_grid(plot_nsp + coord_equal(), plot_bank + coord_equal(), ncol = 2)
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

For clarity, we can visualise each of the clusters separately:

``` r
plot_grid(
    plot_nsp + facet_wrap(~colour_by), 
    plot_bank + facet_wrap(~colour_by), 
    ncol = 2
)
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

<details>
<summary>
Runtime for analysis
</summary>

    #> Time difference of 56.69139 secs

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
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] cowplot_1.1.1               scater_1.26.1              
#>  [3] ggplot2_3.4.0               scuttle_1.8.3              
#>  [5] SpatialExperiment_1.8.0     SingleCellExperiment_1.20.0
#>  [7] SummarizedExperiment_1.28.0 Biobase_2.58.0             
#>  [9] GenomicRanges_1.50.1        GenomeInfoDb_1.34.4        
#> [11] IRanges_2.32.0              S4Vectors_0.36.0           
#> [13] BiocGenerics_0.44.0         MatrixGenerics_1.10.0      
#> [15] matrixStats_0.62.0          Banksy_0.2.5               
#> 
#> loaded via a namespace (and not attached):
#>  [1] bitops_1.0-7              tools_4.2.1              
#>  [3] utf8_1.2.2                R6_2.5.1                 
#>  [5] irlba_2.3.5.1             vipor_0.4.5              
#>  [7] HDF5Array_1.26.0          uwot_0.1.14              
#>  [9] DBI_1.1.3                 colorspace_2.0-3         
#> [11] rhdf5filters_1.10.0       withr_2.5.0              
#> [13] gridExtra_2.3             tidyselect_1.2.0         
#> [15] compiler_4.2.1            cli_3.4.1                
#> [17] BiocNeighbors_1.16.0      DelayedArray_0.24.0      
#> [19] labeling_0.4.2            scales_1.2.1             
#> [21] stringr_1.5.0             digest_0.6.30            
#> [23] dbscan_1.1-11             rmarkdown_2.19           
#> [25] R.utils_2.12.2            aricode_1.0.2            
#> [27] XVector_0.38.0            pkgconfig_2.0.3          
#> [29] htmltools_0.5.4           sparseMatrixStats_1.10.0 
#> [31] highr_0.10                fastmap_1.1.0            
#> [33] limma_3.54.0              rlang_1.1.1              
#> [35] rstudioapi_0.14           DelayedMatrixStats_1.20.0
#> [37] farver_2.1.1              generics_0.1.3           
#> [39] BiocParallel_1.32.5       dplyr_1.0.10             
#> [41] R.oo_1.25.0               RCurl_1.98-1.9           
#> [43] magrittr_2.0.3            BiocSingular_1.14.0      
#> [45] GenomeInfoDbData_1.2.9    Matrix_1.5-3             
#> [47] ggbeeswarm_0.7.1          Rcpp_1.0.9               
#> [49] munsell_0.5.0             Rhdf5lib_1.20.0          
#> [51] fansi_1.0.3               viridis_0.6.2            
#> [53] lifecycle_1.0.3           R.methodsS3_1.8.2        
#> [55] stringi_1.7.8             leidenAlg_1.1.0          
#> [57] yaml_2.3.6                edgeR_3.40.1             
#> [59] zlibbioc_1.44.0           rhdf5_2.42.0             
#> [61] grid_4.2.1                ggrepel_0.9.2            
#> [63] parallel_4.2.1            dqrng_0.3.0              
#> [65] lattice_0.20-45           sccore_1.0.2             
#> [67] beachmat_2.14.0           locfit_1.5-9.6           
#> [69] magick_2.7.3              knitr_1.41               
#> [71] pillar_1.8.1              igraph_1.3.5             
#> [73] rjson_0.2.21              codetools_0.2-18         
#> [75] ScaledMatrix_1.6.0        glue_1.6.2               
#> [77] evaluate_0.19             data.table_1.14.6        
#> [79] vctrs_0.5.1               gtable_0.3.1             
#> [81] assertthat_0.2.1          xfun_0.36                
#> [83] rsvd_1.0.5                DropletUtils_1.18.1      
#> [85] viridisLite_0.4.1         tibble_3.1.8             
#> [87] RcppHungarian_0.2         beeswarm_0.4.0
```

</details>
