
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->
<!--![BANKSY-panel](man/figures/banksy-panel-subtitle.png)-->

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
details on use-cases and methods, see the [preprint]().

## Installation

The *Banksy* package can be installed via `remotes`:

``` r
remotes::install_github("jleechung/Banksy", dependencies = TRUE)
```

## Documentation

Detailed description of *Banksy* functionality and example analyses are
available at the [package webpage]().

*Banksy* comes installed with [documentation]() of main functions and
their usage, along with several vignettes which detail different use
cases:

-   [Working with Banksy objects](): Introduction to the *BanksyObject*
    class which serves as a container for *Banksy*.

-   [Mouse hippocampus MERFISH dataset](): Illustrates a grid search of
    parameters which best cluster cells.

-   [Human dorsolateral prefrontal cortex 10x Visium dataset]():
    Illustrates analysis of multiple spatial transcriptomic datasets.

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
#>         cellID_110 cellID_124 cellID_128 cellID_150 cellID_160
#> Dcx              0          0          0          0          0
#> Sqstm1           3          0          0          1          0
#> Rgs5             0          0          0          0          0
#> Slc1a3           0          0          0          0          0
#> Sparcl1          7          0          1          9          1
#> Notch1           1          0          0          0          0
```

while cell locations should be supplied as a `data.frame`:

``` r
class(locs)
#> [1] "data.frame"
head(locs)
#>                  sdimx    sdimy
#> cellID_110 -13039.7282 16162.19
#> cellID_124  -8470.7282 16153.19
#> cellID_128  -8009.7282 16185.19
#> cellID_150  -2664.7282 16164.19
#> cellID_160   -797.7282 16223.19
#> cellID_166   2200.2718 16208.19
```

Next, create a *BanksyObject* with the expression matrix and cell
locations.

``` r
bank <- BanksyObject(own.expr = expr,
                     cell.locs = locs)
```

We first normalize the expression matrix, compute the neighbor matrix,
and scale the resulting matrices.

``` r
bank <- NormalizeBanksy(bank, normFactor = 100)
bank <- ComputeBanksy(bank)
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
the Leiden algorithm. Specify the following parameters:

-   `resolution`. Leiden clustering resolution.  
-   `k.neighbours`. Number of k neighbours to use for constructing sNN.

``` r
bank <- ClusterBanksy(bank, lambda = c(0, 0.3), method = 'leiden',
                      k.neighbors = 50, resolution = 1.2, seed = 1234)
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

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

For clarity, we can visualise each of the clusters separately with
`wrap = TRUE`:

``` r
plotSpatialFeatures(bank, by = features, type = feature.types, main = main, 
                    pt.size = 0.5, main.size = 15, nrow = 1, ncol = 2, 
                    wrap = TRUE)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

## Contributing

Bug reports, questions, request for enhancements or other contributions
can be raised at the [issue page]().
