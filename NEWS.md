
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
