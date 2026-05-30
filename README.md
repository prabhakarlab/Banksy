
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Overview

BANKSY is a method for clustering spatial omics data by augmenting the
features of each cell with both an average of the features of its
spatial neighbors along with neighborhood feature gradients. By
incorporating neighborhood information for clustering, BANKSY is able to

- improve cell-type assignment in noisy data
- distinguish subtly different cell-types stratified by microenvironment
- identify spatial domains sharing the same microenvironment

BANKSY is applicable to a wide array of spatial technologies (e.g. 10x
Visium, Slide-seq, MERFISH, CosMX, CODEX) and scales well to large
datasets. For more details, check out:

- the [paper](https://www.nature.com/articles/s41588-024-01664-3),
- the [peer review
  file](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-024-01664-3/MediaObjects/41588_2024_1664_MOESM3_ESM.pdf),
- a
  [tweetorial](https://x.com/shyam_lab/status/1762648072360792479?s=20)
  on BANKSY,
- a set of [vignettes](https://prabhakarlab.github.io/Banksy) showing
  basic usage,
- usage compatibility with Seurat
  ([here](https://github.com/satijalab/seurat-wrappers/blob/master/docs/banksy.md)
  and
  [here](https://satijalab.org/seurat/articles/visiumhd_analysis_vignette#identifying-spatially-defined-tissue-domains)),
- a [Python version](https://github.com/prabhakarlab/Banksy_py) of this
  package,
- a [Zenodo archive](https://zenodo.org/records/10258795) containing
  scripts to reproduce the analyses in the paper, and the corresponding
  [GitHub Pages](https://github.com/jleechung/banksy-zenodo) (and
  [here](https://github.com/prabhakarlab/Banksy_py/tree/Banksy_manuscript)
  for analyses done in Python).

**BANKSY now includes a lazy PCA mode (`lazy=TRUE` in `runBanksyPCA`)
that computes PCA directly via an implicit linear operator without
materializing the full BANKSY matrix. This is the default and
recommended mode, scaling to millions of cells with low memory usage.
For analysis on large datasets (\> 1 million samples), we recommend
using BANKSY via SeuratWrappers: see [this
vignette](https://github.com/jleechung/seurat-wrappers/blob/feat-sparse-matmul/docs/banksy.md#scaling-to-large-datasets).**

## Installation

The *Banksy* package can be installed via Bioconductor. This currently
requires R `>= 4.4.0`.

``` r
BiocManager::install('Banksy')
```

To install directly from GitHub instead, use

``` r
remotes::install_github("prabhakarlab/Banksy")
```

To use the legacy version of *Banksy* utilising the `BanksyObject`
class, use

``` r
remotes::install_github("prabhakarlab/Banksy@legacy")
```

*Banksy* is also interoperable with
[*Seurat*](https://satijalab.org/seurat/) via
[*SeuratWrappers*](https://github.com/satijalab/seurat-wrappers).
Documentation on how to run BANKSY on Seurat objects can be found
[here](https://github.com/satijalab/seurat-wrappers/blob/master/docs/banksy.md).
For installation of *SeuratWrappers* with BANKSY version `>= 0.1.6`, run

``` r
remotes::install_github('satijalab/seurat-wrappers')
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
se <- SpatialExperiment(assay = list(counts = gcm), spatialCoords = locs)

# QC based on total counts
qcstats <- perCellQCMetrics(se)
thres <- quantile(qcstats$total, c(0.05, 0.98))
keep <- (qcstats$total > thres[1]) & (qcstats$total < thres[2])
se <- se[, keep]

# Normalization to mean library size
se <- computeLibraryFactors(se)
aname <- "normcounts"
assay(se, aname) <- normalizeCounts(se, log = FALSE)
```

Run `runBanksyPCA` to compute BANKSY PCA embeddings. By default, this
uses a lazy linear operator that computes PCA without materializing the
full BANKSY matrix — computing the kNN graph and applying the BANKSY
transform on the fly during the iterative PCA solver. This is
memory-efficient and scales to millions of cells.

We run *BANKSY* at `lambda=0` corresponding to non-spatial clustering,
and `lambda=0.2` corresponding to *BANKSY* for cell-typing. The number
of spatial neighbors is set by `k_geom`.

> **An important note about choosing the `lambda` parameter for the
> older [Visium v1 / v2 55um
> datasets](https://doi.org/10.1038/s41593-020-00787-0) or the original
> [ST 100um technology](https://doi.org/10.1038/s41596-018-0045-2):**
>
> **Modern high resolution technologies** (Xenium, Visium HD, StereoSeq,
> MERFISH, STARmap PLUS, SeqFISH+, SlideSeq v2, and CosMx, and others):
> we recommend the usual defaults for `lambda`. For cell typing, use
> `lambda = 0.2` (as shown below, or in [this
> vignette](https://prabhakarlab.github.io/Banksy/articles/parameter-selection.html))
> and for [domain
> segmentation](https://prabhakarlab.github.io/Banksy/articles/domain-segment.html),
> use `lambda = 0.8`. These technologies are either imaging based,
> having true single-cell resolution (e.g., MERFISH), or are sequencing
> based, having barcoded spots on the scale of single-cells (e.g.,
> [Visium
> HD](https://www.10xgenomics.com/products/visium-hd-spatial-gene-expression)).
> We find that the usual defaults work well at this measurement
> resolution.
>
> **Older Visium v1/v2 or ST technologies** (with their much lower
> resolution spots, 55um and 100um diameter, respectively): we find that
> `lambda = 0.2` seems to work best for domain segmentation. This could
> be because each spot already measures the average transcriptome of
> several cells in a neighbourhood. It seems that `lambda = 0.2` shares
> enough information between these neighbourhoods to lead to good domain
> segmentation performance. For example, in the [human DLPFC
> vignette](https://prabhakarlab.github.io/Banksy/articles/multi-sample.html),
> we use `lambda = 0.2` on a Visium v1/v2 dataset. Also note that in
> these lower resolution technologies, each spot can have multiple cells
> of different types, and as such *cell-typing* is not defined for them.

``` r
lambda <- c(0, 0.2)

se <- Banksy::runBanksyPCA(se, assay_name = aname, lambda = lambda, k_geom = 15)
#> Computing neighbors...
#> Spatial mode is kNN_median
#> Parameters: k_geom=15
#> Done
#> --- lambda = 0 ---
#> Building sparse weight matrix
#> Computing scaling parameters for own expression
#> Computing clipping excess for own expression
#> Computing scaling params and clipping for H0
#> H0 genes requiring clipping: 0 / 120
#> Clipping corrections: own=0 H0=0 entries
#> Computing BANKSY PCA (20 PCs) via C++ irlba (work=27)
#>   iter=1  mprod=54  sv[20]=7.9454e+01  t=0s
#>   iter=2  mprod=68  sv[20]=9.5131e+01  t=0s
#>   iter=5  mprod=110  sv[20]=1.0435e+02  t=1s
#>   iter=10  mprod=180  sv[20]=1.0568e+02  t=1s
#>   iter=11  mprod=194  sv[20]=1.0568e+02  t=1s
#>   Converged: iter=11, mprod=194
#> --- lambda = 0.2 ---
#> Building sparse weight matrix
#> Computing scaling parameters for own expression
#> Computing clipping excess for own expression
#> Computing scaling params and clipping for H0
#> H0 genes requiring clipping: 0 / 120
#> Clipping corrections: own=0 H0=0 entries
#> Computing BANKSY PCA (20 PCs) via C++ irlba (work=27)
#>   iter=1  mprod=54  sv[20]=6.2079e+01  t=0s
#>   iter=2  mprod=68  sv[20]=8.2518e+01  t=0s
#>   iter=5  mprod=110  sv[20]=9.5141e+01  t=1s
#>   iter=9  mprod=166  sv[20]=9.6004e+01  t=1s
#>   Converged: iter=9, mprod=166
#> Done.
```

Next, compute UMAP and perform clustering:

``` r
set.seed(1000)
se <- Banksy::runBanksyUMAP(se, use_agf = FALSE, lambda = lambda)
se <- Banksy::clusterBanksy(se, use_agf = FALSE, lambda = lambda, resolution = 1.2)
```

Different clustering runs can be relabeled to minimise their differences
with `connectClusters`:

``` r
se <- Banksy::connectClusters(se)
#> clust_M0_lam0.2_k50_res1.2 --> clust_M0_lam0_k50_res1.2
```

Visualise the clustering output for non-spatial clustering (`lambda=0`)
and BANKSY clustering (`lambda=0.2`).

``` r
cnames <- colnames(colData(se))
cnames <- cnames[grep("^clust", cnames)]
colData(se) <- cbind(colData(se), spatialCoords(se))

plot_nsp <- plotColData(se,
    x = "sdimx", y = "sdimy",
    point_size = 0.6, colour_by = cnames[1]
)
plot_bank <- plotColData(se,
    x = "sdimx", y = "sdimy",
    point_size = 0.6, colour_by = cnames[2]
)


plot_grid(plot_nsp + coord_equal(), plot_bank + coord_equal(), ncol = 2)
```

<img src="man/figures/README-unnamed-chunk-13-1.png" width="100%" />

For clarity, we can visualise each of the clusters separately:

``` r
plot_grid(
    plot_nsp + facet_wrap(~colour_by),
    plot_bank + facet_wrap(~colour_by),
    ncol = 2
)
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

Visualize UMAPs of the non-spatial and BANKSY embedding:

``` r
rdnames <- reducedDimNames(se)

umap_nsp <- plotReducedDim(se,
    dimred = grep("UMAP.*lam0$", rdnames, value = TRUE),
    colour_by = cnames[1]
)
umap_bank <- plotReducedDim(se,
    dimred = grep("UMAP.*lam0.2$", rdnames, value = TRUE),
    colour_by = cnames[2]
)
plot_grid(
    umap_nsp,
    umap_bank,
    ncol = 2
)
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

<details>
<summary>
Runtime for analysis
</summary>

    #> Time difference of 35.0907 secs

</details>
<details>
<summary>
Session information
</summary>

``` r
sessionInfo()
#> R version 4.5.1 (2025-06-13)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Sonoma 14.6.1
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: America/Los_Angeles
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] cowplot_1.2.0               scater_1.37.0              
#>  [3] ggplot2_3.5.2               scuttle_1.19.0             
#>  [5] SpatialExperiment_1.18.1    SingleCellExperiment_1.30.1
#>  [7] SummarizedExperiment_1.39.1 Biobase_2.69.0             
#>  [9] GenomicRanges_1.61.1        Seqinfo_0.99.1             
#> [11] IRanges_2.43.0              S4Vectors_0.47.0           
#> [13] BiocGenerics_0.55.0         generics_0.1.4             
#> [15] MatrixGenerics_1.21.0       matrixStats_1.5.0          
#> [17] Banksy_1.9.2               
#> 
#> loaded via a namespace (and not attached):
#>  [1] beeswarm_0.4.0      gtable_0.3.6        rjson_0.2.23       
#>  [4] xfun_0.52           ggrepel_0.9.6       lattice_0.22-7     
#>  [7] vctrs_0.6.5         tools_4.5.1         parallel_4.5.1     
#> [10] tibble_3.3.0        sccore_1.0.6        pkgconfig_2.0.3    
#> [13] BiocNeighbors_2.3.1 Matrix_1.7-3        data.table_1.17.6  
#> [16] RColorBrewer_1.1-3  lifecycle_1.0.4     compiler_4.5.1     
#> [19] farver_2.1.2        aricode_1.0.3       codetools_0.2-20   
#> [22] vipor_0.4.7         GenomeInfoDb_1.45.7 htmltools_0.5.8.1  
#> [25] yaml_2.3.10         pillar_1.11.0       crayon_1.5.3       
#> [28] BiocParallel_1.43.4 uwot_0.2.3          DelayedArray_0.35.2
#> [31] dbscan_1.2.2        viridis_0.6.5       magick_2.8.7       
#> [34] abind_1.4-8         mclust_6.1.1        rsvd_1.0.5         
#> [37] tidyselect_1.2.1    digest_0.6.37       BiocSingular_1.24.0
#> [40] dplyr_1.1.4         labeling_0.4.3      fastmap_1.2.0      
#> [43] grid_4.5.1          cli_3.6.5           SparseArray_1.9.0  
#> [46] magrittr_2.0.3      leidenAlg_1.1.5     S4Arrays_1.9.1     
#> [49] dichromat_2.0-0.1   withr_3.0.2         scales_1.4.0       
#> [52] UCSC.utils_1.5.0    ggbeeswarm_0.7.2    rmarkdown_2.29     
#> [55] XVector_0.49.0      httr_1.4.7          igraph_2.1.4       
#> [58] gridExtra_2.3       ScaledMatrix_1.16.0 beachmat_2.25.1    
#> [61] evaluate_1.0.4      RcppHungarian_0.3   knitr_1.50         
#> [64] RcppAnnoy_0.0.22    viridisLite_0.4.2   irlba_2.3.5.1      
#> [67] rlang_1.1.6         Rcpp_1.1.0          glue_1.8.0         
#> [70] jsonlite_2.0.0      R6_2.6.1
```

</details>
