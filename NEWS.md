
# Version 1.9.2

+ Lazy PCA mode in runBanksyPCA (lazy=TRUE, now the default) computes PCA via an implicit linear operator without materializing the full BANKSY matrix, scaling to millions of cells with low memory usage
+ C++ backend (pca_backend="cpp", default) implements irlba with OpenMP-parallelized sparse matrix-vector products and column-wise H0 scaling that computes neighborhood statistics without materializing the gene-by-cell product matrix
+ The lazy path supports multi-sample analysis with per-group kNN, within-group scaling, optional parallel kNN via mclapply, and multiple lambda values in a single call with kNN computed once and reused. 
+ At 96% sparsity and 5k features with average group size ~400k, expected runtime / peak mem for runBanksyPCA is 1.9 min / 27 GB at 1M cells, 2.8 min / 39 GB at 2M, 8.2 min / 84 GB at 5M, and 13.0 min / 163 GB at 10M cells, scaling approximately linearly in both runtime and memory.

# Version 0.99.8

+ Add feature scaling options for PCA and UMAP

# Version 0.99.7

+ SCE/SPE-compatible

# Version 0.1.6

+ Fix compatibility with SeuratWrappers

# Version 0.1.5

+ Implemented SmoothLabels for k-nearest neighbors cluster label smoothing
+ Parallel clustering for Leiden graph-based clustering
+ Version depedency on leidenAlg (>= 1.1.0) for compatibility with igraph (>= 1.5.0)
+ Seed setting for clustering
+ Neighborhood sampling for computing neighborhood feature matrices. See arguments 
`sample_size`, `sample_renorm` and `seed` in function ComputeBanksy

# Version 0.1.4

+ Implemented Azimuthal Gabor filters in ComputeBanksy, with the number of 
harmonics determined by the `M` argument. To obtain similar results as version 
0.1.3, use `M=0`
+ RunPCA and RunUMAP renamed to RunBanksyPCA and RunBanksyUMAP respectively to 
avoid namespace collisions with other packages
+ Changed some function argument names (e.g. `spatialMode` to `spatial_mode`)

# Version 0.1.3

+ Interoperability with SingleCellExperiment with asBanksyObject
+ Passing Bioc and R CMD checks

# Version 0.1.2

+ BanksyObject: dimensionality reduction slot name changed to reduction
+ RunPCA and RunUMAP
+ Plotting is generalised (plotReduction instead of plotUMAP / plotPCA)
+ ConnectClusters improved and returns BanksyObject
+ ARI computation and plotting 

# Version 0.0.9 

+ Legacy version
