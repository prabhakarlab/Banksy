
# BanksyObject single dataset --------------------------------------------------

library(testthat)
library(matrixStats)
library(plyr)

data <- Banksy::simulateDataset()
gcm <- data$gcm
locs <- data$locs
mdata <- data$meta
mdata$total_count <- colSums(gcm)

## Construction

test_that('Construct Banksy without expression', {
    expect_error(BanksyObject(cell.locs = locs, meta.data = mdata))
})


test_that('Construct Banksy without locations', {
    expect_error(BanksyObject(own.expr = gcm, meta.data = mdata))
})


test_that('Construct Banksy with errant metadata', {
    expect_error(BanksyObject(own.expr = gcm,
                              cell.locs = locs,
                              meta.data = mdata[seq_len(10),,drop=FALSE]))
})

test_that('Construct Banksy', {
    expect_s4_class(BanksyObject(own.expr = gcm,
                                 cell.locs = locs), 'BanksyObject')
    expect_s4_class(BanksyObject(own.expr = gcm,
                              cell.locs = locs,
                              meta.data = mdata), 'BanksyObject')
})

bank <- BanksyObject(own.expr = gcm,
                     cell.locs = locs,
                     meta.data = mdata)

test_that('Show', {
    expect_output(show(bank))
})

test_that('Head', {
    expect_output(head(bank))
})

# Normalization own

bank <- NormalizeBanksy(bank, normFactor = 100)

test_that('Normalization for own', {
    
    expect_error(NormalizeBanksy(bank, assay = 'hello'))
    
    expect_equal(as.numeric(round(colSums(bank@own.expr))), rep(100,501))
    expect_equal(class(NormalizeBanksy(bank, logNorm = TRUE)@own.expr), 
                 c('matrix', 'array'))
    
    bankOwn <- NormalizeBanksy(bank, assay = 'own')
    expect_equal(as.numeric(round(colSums(bankOwn@own.expr))), rep(100,501))
})

# Compute neighbors

bank <- ComputeBanksy(bank, verbose = TRUE)

test_that('Compute neighbors', {
    
    expect_equal(dim(bank@nbr.expr), c(100, 501))
    bank <- ComputeBanksy(bank, spatialMode = 'kNN_rank', verbose = TRUE)
    expect_equal(dim(bank@nbr.expr), c(100, 501))
    bank <- ComputeBanksy(bank, spatialMode = 'kNN_unif', verbose = TRUE)
    expect_equal(dim(bank@nbr.expr), c(100, 501))
    bank <- ComputeBanksy(bank, spatialMode = 'rNN_gauss', verbose = TRUE)
    expect_equal(dim(bank@nbr.expr), c(100, 501))
    bank <- ComputeBanksy(bank, spatialMode = 'kNN_rn', verbose = TRUE)
    expect_equal(dim(bank@nbr.expr), c(100, 501))
})

test_that('Compute neighbors with invalid dimensions', {
    expect_error(ComputeBanksy(bank, dimensions = c('')))
})

test_that('Compute neighbors with less than 2 dimensions', {
    bankDummy <- bank
    bankDummy@cell.locs <- bank@cell.locs[,1,drop=FALSE]
    expect_error(ComputeBanksy(bankDummy))
})

test_that('Compute neighbors with too many neighbors', {
    expect_error(ComputeBanksy(bank, k_geom = 1000))
})

test_that('Compute neighbors with invalid method', {
    expect_error(ComputeBanksy(bank, spatialMode = 'dummy'))
})

# Normalization nbr

test_that('Normalization for neighbors', {

    bankNbr <- bank
    bankNbr <- NormalizeBanksy(bankNbr, assay = 'nbr')
    expect_equal(as.numeric(round(colSums(bankNbr@nbr.expr))), rep(100,501))
    
})

# Scale

bank <- ScaleBanksy(bank)

test_that('Scaling', {
    expect_error(ScaleBanksy(bank, assay = 'hello'))
    bankOwn <- ScaleBanksy(bank, assay = 'own')
    expect_equal(as.numeric(rowMeans2(bankOwn@own.expr)), rep(0, 100))
    expect_equal(as.numeric(rowSds(bankOwn@own.expr)), rep(1, 100))
    bankNbr <- ScaleBanksy(bankOwn, assay = 'nbr')
    expect_equal(as.numeric(rowMeans2(bankNbr@nbr.expr)), rep(0, 100))
    expect_equal(as.numeric(rowSds(bankNbr@nbr.expr)), rep(1, 100))
})

# Reduction

bank <- RunPCA(bank, lambda = 0.25, npcs = 15)

test_that('PCA', {
    expect_equal(names(bank@reduction),'pca_0.25')
    expect_equal(dim(bank@reduction$pca_0.25$x), c(501, 15))
})

test_that('Scree', {
    expect_error(plotScree(bank, lambda = 0.2))
    expect_equal(class(plotScree(bank, lambda = 0.25)), c('gg', 'ggplot'))
})

bank <- RunUMAP(bank, lambda = 0.25, pca = TRUE, npcs = 10)

test_that('UMAP', {
    expect_equal(names(bank@reduction), c('pca_0.25', 'umap_0.25'))
    expect_equal(dim(bank@reduction$umap_0.25), c(501,2))
})

test_that('UMAP on BANKSY matrix', {
    expect_equal(names(RunUMAP(bank, lambda = 0.25, pca = FALSE)@reduction),
                 c('pca_0.25', 'umap_0.25'))
})

test_that('UMAP errors', {
    expect_error(RunUMAP(bank, lambda = 0.2))
    expect_error(RunUMAP(bank, lambda = 0.25, npcs = 1))
    expect_warning(RunUMAP(bank, lambda = 0.25, npcs = 30))
})

# Cluster

test_that('Cluster method check', {
    expect_equal(length(Banksy:::checkMethod(c('leiden', 'louvain'))), 1)
})

test_that('Leiden clustering unspec. parameters', {
    expect_error(ClusterBanksy(bank, pca = TRUE, npcs = 15, method = 'leiden'))
    expect_error(ClusterBanksy(bank, pca = TRUE, npcs = 15, 
                               method = 'leiden', k.neighbors = 10))
    expect_error(ClusterBanksy(bank, pca = TRUE, npcs = 15, 
                               method = 'leiden', resolution = 10))
})

test_that('Louvain clustering unspec. parameters', {
    expect_error(ClusterBanksy(bank, pca = TRUE, npcs = 15, method = 'louvain'))
})

test_that('k means clustering unspec. parameters', {
    expect_error(ClusterBanksy(bank, pca = TRUE, npcs = 15, method = 'kmeans'))
})

test_that('Mclust clustering unspec. parameters', {
    expect_error(ClusterBanksy(bank, pca = TRUE, npcs = 15, method = 'mclust'))
})


bank <- ClusterBanksy(bank, lambda = 0.25, pca = TRUE, npcs = 5, 
                      method = 'leiden', k.neighbors = 30, 
                      resolution = c(0.5, 1, 1.5))


test_that('Cluster with different methods', {
    expect_true('clust_lam0.25_k30_res0.5' %in% clust.names(bank))
    bank <- ClusterBanksy(bank, lambda = 0.25, pca = TRUE, npcs = 5, 
                          method = 'louvain', k.neighbors = c(20,30))
    expect_true('clust_lam0.25_k30_louvain' %in% clust.names(bank))
    bank <- ClusterBanksy(bank, lambda = 0.25, pca = TRUE, npcs = 5,
                          method = 'kmeans', kmeans.centers = c(3,4))
    expect_true('clust_lam0.25_kmeans3' %in% clust.names(bank))
    expect_equal(length(unique(bank@meta.data$clust_lam0.25_kmeans3)),3)
    bank <- ClusterBanksy(bank, lambda = 0.25, pca = TRUE, npcs = 5,
                          method = 'mclust', mclust.G = c(3,4))
    expect_true('clust_lam0.25_mclust3' %in% clust.names(bank))
    expect_equal(length(unique(bank@meta.data$clust_lam0.25_kmeans3)),3)
})

test_that('Insuffient PCs', {
    expect_error(ClusterBanksy(bank, lambda = 0.2))
    expect_error(ClusterBanksy(bank, lambda = 0.25, pca = TRUE, npcs = 30))
})

test_that('Cluster with matrix', {
    expect_equal(clust.names(ClusterBanksy(bank, lambda = 0.25, pca = FALSE, 
                                           method = 'leiden', k.neighbors = 30, 
                                           resolution = c(0.5,1,1.5))),
                 clust.names(bank))
})

# Connect clusters

test_that('Connect clusters', {
    expect_error(ConnectClusters(bank, map.to = 'dummy'))
    expect_equal(clust.names(ConnectClusters(bank, map.to = NULL)), 
                 clust.names(bank))
    expect_equal(clust.names(ConnectClusters(bank, 
                                             map.to = clust.names(bank)[1])), 
                 clust.names(bank))
    expect_s4_class(ConnectClusters(bank), 'BanksyObject')
})

test_that('ARI computation', {
    bankTmp <- bank
    bankTmp@meta.data <- bank@meta.data[,seq_len(4)]
    expect_error(getARI(bankTmp))
    expect_equal(class(getARI(bank)), c('matrix', 'array'))
})


# Plotting

test_that('Plot reduction', {
    expect_error(plotReduction(bank, reduction = 'pca_0'))
    expect_error(plotReduction(bank, reduction = 'pca_0.25', 
                               by = clust.names(bank)[1]))
    expect_error(plotReduction(bank, reduction = 'pca_0.25', 
                               by = clust.names(bank)[1], 
                               type = c('discrete', 'continuous')))

    expect_warning(plotReduction(bank, reduction = 'pca_0.25', 
                                 by = clust.names(bank), type = 'discrete'))
    
    expect_equal(class(plotReduction(bank, reduction = 'pca_0.25')), 
                 c('gg', 'ggplot'))
    expect_equal(class(plotReduction(bank, reduction = 'pca_0.25', 
                                     by = clust.names(bank)[1], 
                                     type = 'discrete')), 
                 c('gg', 'ggplot'))
    expect_equal(class(plotReduction(bank, reduction = 'pca_0.25',
                                     by = 'total_count', type = 'continuous')), 
                 c('gg', 'ggplot'))
    
    expect_equal(class(plotReduction(bank, reduction = 'umap_0.25')), 
                 c('gg', 'ggplot'))
    
})


test_that('Plot spatial', {
    
    cnm <- clust.names(bank)[1]
    
    expect_equal(class(plotSpatial(bank)), c('gg', 'ggplot'))
    expect_equal(class(plotSpatial(bank, by = cnm, 
                                   type = 'discrete')), 
                 c('gg', 'ggplot'))
    
    
    palette <- seq_len(3)
    names(palette) <- seq_len(3)
    expect_equal(class(plotSpatial(bank, by = cnm, 
                                   type = 'discrete', col.discrete = palette)), 
                 c('gg', 'ggplot'))
    palette <- seq_len(2)
    expect_error(plotSpatial(bank, by = cnm, type = 'discrete', 
                             col.discrete = palette))
    palette <- seq_len(3)
    expect_error(plotSpatial(bank, by = cnm, type = 'discrete', 
                             col.discrete = palette))
    palette <- seq_len(3)
    names(palette) <- seq(5,7)
    expect_error(plotSpatial(bank, by = cnm, type = 'discrete', 
                             col.discrete = palette))
    
    bankChar <- bank
    bankChar@meta.data$LabelChar <- mapvalues(
        bank@meta.data$Label, from = seq_len(3), to = letters[seq_len(3)])
    expect_equal(class(plotSpatial(bankChar, by = 'LabelChar', 
                                   type = 'discrete')), 
                 c('gg', 'ggplot'))
    
    palette <- seq_len(3)
    names(palette) <- letters[seq_len(3)]
    expect_equal(class(plotSpatial(bankChar, by = 'LabelChar', 
                                   type = 'discrete', col.discrete = palette)), 
                 c('gg', 'ggplot'))
    
    expect_equal(class(plotSpatial(bank, by = 'total_count', 
                                   type = 'continuous')), 
                 c('gg', 'ggplot'))
    expect_equal(class(plotSpatial(bank, by = cnm, 
                                   type = 'discrete', wrap = TRUE)),
                 c('gg', 'ggplot'))
    expect_equal(class(plotSpatial(bank, by = 'gene_2', type = 'continuous')),
                 c('gg', 'ggplot'))
    expect_equal(class(plotSpatial(bank, by = 'gene_2', type = 'continuous',
                                   col.midpoint = 0, col.highpoint = 1, 
                                   col.lowpoint = -1, pt.size = 2)),
                 c('gg', 'ggplot'))
    
    bankZdim <- bank
    bankZdim@cell.locs$sdimz <- 0.26
    expect_equal(class(plotSpatial(bankZdim)), c('gg', 'ggplot'))
    expect_equal(class(plotSpatial(bankZdim, by = 'Label', type = 'discrete',
                                   wrap = TRUE)), c('gg', 'ggplot'))
    
    
})

test_that('Plot spatial features', {
    
    bankSub <- SubsetBanksy(bank, cells = sample(colnames(bank@own.expr), 100))
    expect_error(plotSpatialFeatures(bankSub, by = clust.names(bank), 
                                     type = 'discrete'))
    expect_null(plotSpatialFeatures(
        bankSub, by = clust.names(bank), 
        type = rep('discrete', length(clust.names(bank)))))
    expect_equal(class(plotSpatialFeatures(
        bankSub, by = clust.names(bank), return.plot = TRUE,
        type = rep('discrete', length(clust.names(bank))))), 'list')
})

test_that('Plot heatmap', {
    expect_error(plotHeatmap(bank, assay = 'dummy'))
    expect_error(plotHeatmap(bank, assay = 'banksy'))
    expect_error(plotHeatmap(bank,
                             col = c('blue', 'white', 'red'),
                             col.breaks = c(-1,1)))
    expect_error(plotHeatmap(bank,
                             annotate = TRUE,
                             annotate.by = 'dummy'))
    expect_error(plotHeatmap(bank,
                             annotate = TRUE,
                             annotate.by = 'Label',
                             order.by = 'dummy'))
    expect_error(plotHeatmap(bank,
                             annotate = TRUE,
                             annotate.by = 'Label',
                             barplot.by = 'dummy'))

    expect_s4_class(plotHeatmap(bank, assay = 'own.expr'), 'Heatmap')
    expect_s4_class(plotHeatmap(bank, assay = 'own.expr', max.cols = 100), 
                    'Heatmap')
    expect_s4_class(plotHeatmap(bank, assay = 'nbr.expr'), 'Heatmap')
    expect_s4_class(plotHeatmap(bank, assay = 'banksy', lambda = 0.25,
                                annotate = TRUE, annotate.by = 'Label'),
                    'Heatmap')
    expect_s4_class(plotHeatmap(bank,
                                annotate = TRUE,
                                annotate.by = c('Label', clust.names(bank)),
                                barplot.by = 'total_count'), 'HeatmapList')
})

test_that('Plot ARI', {
    expect_equal(class(plotARI(bank)), c('gg', 'ggplot'))
})

test_that('Plot alluvia', {
    bankSub <- bank
    bankSub@meta.data <- bank@meta.data[,seq_len(4)]
    expect_error(plotAlluvia(bankSub))
    expect_error(plotAlluvia(bank, flow.colors = 'red'))
    
    
    expect_equal(class(plotAlluvia(bank, max.cells = NULL)), c('gg', 'ggplot'))
    expect_equal(class(plotAlluvia(bank, max.cells = 100, flow.colors = 1:8)), 
                 c('gg', 'ggplot'))
    expect_equal(class(plotAlluvia(bank)), c('gg', 'ggplot'))
})


# Handling

test_that('Split Banksy', {
    expect_error(SplitBanksy(bank, by = 'hello'))
    expect_error(SplitBanksy(bank, by = 'Label', names = letters[seq_len(4)]))
    expect_s4_class(SplitBanksy(bank, by = 'Label'), 'BanksyObject')
    bankSplit <- SplitBanksy(bank, by = 'Label')
    expect_equal(length(bankSplit@own.expr), 3)
    expect_equal(length(bankSplit@nbr.expr), 3)
    expect_equal(length(bankSplit@cell.locs), 3)
})

# Subset
test_that('Subsetting cells', {
    bankSubsetCells <- SubsetBanksy(
        bank, cells = sample(bank@meta.data$cell_ID, 50))
    expect_equal(dim(bankSubsetCells@own.expr), c(100, 50))
    expect_equal(dim(bankSubsetCells@nbr.expr), c(100, 50))
    expect_equal(nrow(bankSubsetCells@reduction$pca_0.25$x), 50)
    expect_equal(nrow(bankSubsetCells@reduction$umap_0.25), 50)
    expect_equal(dim(bankSubsetCells@cell.locs), c(50, 2))
    expect_equal(dim(bankSubsetCells@meta.data[,clust.names(bank)]), c(50, 3))
})

test_that('Subsetting genes', {
    bankSubsetGenes <- SubsetBanksy(
        bank, features = sample(rownames(bank@own.expr), 50))
    expect_equal(dim(bankSubsetGenes@own.expr), c(50, 501))
    expect_equal(dim(bankSubsetGenes@nbr.expr), c(50, 501))
    expect_equal(dim(bankSubsetGenes@cell.locs), c(501, 2))
    expect_equal(dim(bankSubsetGenes@meta.data[,clust.names(bank)]), c(501, 3))
})

test_that('Subsetting by metadata', {
    bankSubsetMeta <- SubsetBanksy(bank, metadata = Label == 1)
    expect_equal(dim(bankSubsetMeta@own.expr), c(100, 184))
    expect_equal(dim(bankSubsetMeta@nbr.expr), c(100, 184))
    expect_equal(nrow(bankSubsetMeta@reduction$pca_0.25$x), 184)
    expect_equal(nrow(bankSubsetMeta@reduction$umap_0.25), 184)
    expect_equal(dim(bankSubsetMeta@cell.locs), c(184, 2))
    expect_equal(dim(bankSubsetMeta@meta.data[,clust.names(bank)]), c(184, 3))
})

test_that('Subsetting by dimensions', {
    bankSubsetDims <- SubsetBanksy(
        bank, dims = (sdimx > -30 & sdimx < 30) & (sdimy > -30 & sdimy < 30))
    expect_equal(dim(bankSubsetDims@own.expr), c(100, 304))
    expect_equal(nrow(bankSubsetDims@reduction$pca_0.25$x), 304)
    expect_equal(nrow(bankSubsetDims@reduction$umap_0.25), 304)
    expect_equal(dim(bankSubsetDims@cell.locs), c(304, 2))
    expect_equal(dim(bankSubsetDims@meta.data[,clust.names(bank)]), c(304, 3))
})

# Slots

test_that('Slots', {
    expect_equal(dim(bank@own.expr), c(100, 501))
    expect_equal(dim(bank@nbr.expr), c(100, 501))
    expect_equal(dim(bank@cell.locs), c(501, 2))
    expect_equal(dim(bank@meta.data[,clust.names(bank)]), c(501, 3))
    
    expect_equal(class(own.expr(bank)), c('matrix', 'array'))
    expect_equal(class(nbr.expr(bank)), c('matrix', 'array'))
    expect_equal(class(cell.locs(bank)), 'data.frame')
    expect_equal(class(meta.data(bank)), 'data.frame')
    expect_equal(class(reduction(bank)), 'list')
    
    own.expr(bank) <- bank@own.expr
    nbr.expr(bank) <- bank@nbr.expr
    cell.locs(bank) <- bank@cell.locs
    meta.data(bank) <- bank@meta.data
    reduction(bank) <- bank@reduction
})
