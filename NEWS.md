
# Version 1.9.4

+ The lazy path in `runBanksyPCA` is now deterministic. Both `pca_backend` options start the iterative solver from a fixed vector instead of a random one, so repeated runs give identical embeddings and the caller's RNG stream is left untouched. Previously the C++ backend called `set.seed(42)` internally, overriding the user's `seed`, and the R backend was unseeded

# Version 1.9.2

+ Lazy PCA mode in runBanksyPCA (lazy=TRUE, now the default) computes PCA via an implicit linear operator without materializing the full BANKSY matrix, scaling to millions of cells with low memory usage
+ C++ backend (pca_backend="cpp", default) implements irlba with OpenMP-parallelized sparse matrix-vector products and column-wise H0 scaling that computes neighborhood statistics without materializing the gene-by-cell product matrix
+ The lazy path supports multi-sample analysis with per-group kNN, within-group scaling, optional parallel kNN via mclapply, and multiple lambda values in a single call with kNN computed once and reused. 
+ The lazy path reads on-disk matrices (e.g. BPCells) without coercing them to memory, allowing datasets beyond R's 2^31 sparse non-zero limit
+ OpenMP regions in the C++ backend are gated on a minimum work size, so small datasets are not slowed by thread creation when the thread count exceeds the available cores. Override with `BANKSY_OMP_MIN_WORK`
+ `runBanksyPCA` with the lazy path has been benchmarked to 25 million cells. Tiled Xenium Prime 5K, `lambda=0.2`, `k_geom=30`, `npcs=50`, per-sample kNN and scaling, on 8 CPUs:

| Cells | Features | Sparsity | Samples | Runtime | Peak memory |
|---:|---:|---:|---:|---:|---:|
| 1,221,372 | 5,101 | 96.18% | 3 | 1.6 min | 22.1 GB |
| 2,035,620 | 5,101 | 96.18% | 5 | 2.5 min | 38.9 GB |
| 5,292,612 | 5,101 | 96.18% | 13 | 6.6 min | 83.9 GB |
| 10,178,100 | 5,101 | 96.18% | 25 | 12.5 min | 163.0 GB |
| 25,241,688 | 5,101 | 96.18% | 62 | 34.6 min | 162.4 GB |

Runtime scales approximately linearly in cell count (`1.3 * n^0.97` minutes, n in millions). Memory scales linearly up to 10M; the 25M run uses on-disk BPCells storage, which caps resident memory at the cost of runtime.

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
