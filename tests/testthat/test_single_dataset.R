
library(testthat)
library(Banksy)
suppressWarnings(library(matrixStats))
library(plyr)

data <- Banksy:::simulateDataset(seed = 42)
gcm<- data[[1]]
locs <- data[[2]]
clusters <- data[[3]]
celltype <- mapvalues(x = clusters,
                      from = unique(clusters),
                      to = letters[seq_len(length(unique(clusters)))])

mdata <- data.frame(clusters = clusters,
                     count = colSums(gcm),
                     celltype = celltype)

##### BanksyObject single dataset ----------------------------------------------

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
                              cell.locs = locs,
                              meta.data = mdata), 'BanksyObject')
})

bank <- BanksyObject(own.expr = gcm,
                     cell.locs = locs,
                     meta.data = mdata)

test_that('Slots', {
    expect_equal(dim(bank@own.expr), c(100, 100))
    expect_equal(dim(bank@nbr.expr), NULL)
    expect_equal(dim(bank@cell.locs), c(100, 2))
    expect_equal(dim(bank@meta.data), c(100, 4))
})

test_that('Show', {
    expect_output(show(bank))
})

test_that('Head', {
    expect_output(head(bank))
})

set.seed(1234)
bankSubset <- SubsetBanksy(bank,
                           cells = sample(bank@meta.data$cell_ID, 50),
                           features = sample(rownames(bank@own.expr), 40),
                           metadata = clusters == 1,
                           dims = sdimx < 800 & sdimy < 1300)

test_that('Subsetting', {
    expect_equal(dim(bankSubset@own.expr), c(40, 6))
    expect_equal(dim(bankSubset@nbr.expr), NULL)
    expect_equal(dim(bankSubset@cell.locs), c(6, 2))
    expect_equal(dim(bankSubset@meta.data), c(6, 4))
})

bank <- NormalizeBanksy(bank)

test_that('Normalization', {
    expect_equal(as.numeric(colSums(bank@own.expr)), rep(100,100))
})

bank <- ComputeBanksy(bank, spatialMode = 'kNN_r', verbose=TRUE)

test_that('Compute neighbors', {
    expect_equal(dim(bank@nbr.expr), c(100, 100))
    bank <- ComputeBanksy(bank, spatialMode = 'kNN_rank', verbose=TRUE)
    expect_equal(dim(bank@nbr.expr), c(100, 100))
    bank <- ComputeBanksy(bank, spatialMode = 'kNN_unif', verbose=TRUE)
    expect_equal(dim(bank@nbr.expr), c(100, 100))
    bank <- ComputeBanksy(bank, spatialMode = 'rNN_gauss', verbose=TRUE)
    expect_equal(dim(bank@nbr.expr), c(100, 100))
})

test_that('Compute neighbors with invalid dimensions', {
    expect_error(ComputeBanksy(bank, dimensions = c('dummy')))
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

bank <- ScaleBanksy(bank)

test_that('Scaling', {
    expect_equal(as.numeric(rowMeans2(bank@own.expr)), rep(0, 100))
    expect_equal(as.numeric(rowMeans2(bank@nbr.expr)), rep(0, 100))
    expect_equal(as.numeric(rowSds(bank@own.expr)), rep(1, 100))
    expect_equal(as.numeric(rowSds(bank@nbr.expr)), rep(1, 100))
})

bank <- ClusterBanksy(bank, npcs = 5, resolution = c(0.5, 1, 1.5))

test_that('Cluster', {
    expect_true('res1_lam0.25_k40' %in% names(bank@meta.data))
    expect_true(all(c('pca_0.25', 'umap_0.25') %in% names(bank@dim.reduction)))
    expect_equal(dim(bank@dim.reduction$pca_0.25), c(100, 5))
    expect_equal(dim(bank@dim.reduction$umap_0.25),c(100, 2))
})

test_that('Plot UMAP', {
    expect_warning(plotUMAP(bank, by = c('clusters', 'dummy'),
                            reduction = 'umap_0.25'))
    expect_warning(plotUMAP(bank, by = 'clusters',
                            reduction = c('umap_0.25', 'umap_0')))
    expect_error(plotUMAP(bank, by = 'dummy', reduction = ''))
    expect_error(plotUMAP(bank, by = 'clusters', reduction = 'umap_0'))
    expect_error(plotUMAP(bank, by = 'clusters', reduction = 'umap_0.25',
                          col = 'black'))
    expect_equal(class(plotUMAP(bank, by = 'clusters',
                                reduction = 'umap_0.25')), c('gg', 'ggplot'))
})

test_that('Plot Spatial Dimensions', {
    expect_warning(plotSpatialDims(bank, by = c('clusters', 'dummy')))
    expect_error(plotSpatialDims(bank, by = 'dummy'))
    expect_error(plotSpatialDims(bank, by = 'clusters', col = 'black'))
    expect_equal(class(plotSpatialDims(bank, by = 'clusters', col = 1:10)),
                 c('gg', 'ggplot'))
    expect_equal(class(plotSpatialDims(bank, by = 'clusters')),
                 c('gg', 'ggplot'))
    expect_equal(class(plotSpatialDims(bank, by = 'celltype')),
                 c('gg', 'ggplot'))
})

test_that('Plot own expression heatmap', {
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
                             annotate.by = 'clusters',
                             order.by = 'dummy'))
    expect_error(plotHeatmap(bank,
                             annotate = TRUE,
                             annotate.by = 'clusters',
                             barplot.by = 'dummy'))

    expect_s4_class(plotHeatmap(bank, assay = 'own.expr'), 'Heatmap')
    expect_s4_class(plotHeatmap(bank, assay = 'nbr.expr'), 'Heatmap')
    expect_s4_class(plotHeatmap(bank, assay = 'banksy', lambda = 0.25,
                                annotate = TRUE, annotate.by = 'clusters'),
                    'Heatmap')
    expect_s4_class(plotHeatmap(bank,
                                annotate = TRUE,
                                annotate.by = c('clusters', 'celltype'),
                                barplot.by = 'count'), 'HeatmapList')
})


bank <- SplitBanksy(bank, by = 'clusters')

test_that('Split Banksy', {
    expect_s4_class(bank, 'BanksyObject')
    expect_equal(length(bank@own.expr), 5)
    expect_equal(length(bank@nbr.expr), 5)
    expect_equal(length(bank@cell.locs), 5)
})

test_that('Connect clusters', {
    expect_s4_class(ConnectClusters(bank)$BanksyObject, 'BanksyObject')
})



