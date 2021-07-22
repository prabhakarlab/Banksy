
library(testthat)
library(Banksy)
suppressWarnings(library(matrixStats))
library(plyr)

data <- Banksy:::simulateDataset(n_clusters = 3,
                                 n_cells = 75,
                                 n_genes = 75,
                                 seed = 42)
gcm1 <- data[[1]]
locs1 <- data[[2]]
clusters <- data[[3]]
celltype <- mapvalues(x = clusters,
                      from = unique(clusters),
                      to = letters[seq_len(length(unique(clusters)))])
mdata1 <- data.frame(clusters = clusters,
                    count = colSums(gcm1),
                    celltype = celltype)

data <- Banksy:::simulateDataset(n_clusters = 3,
                                 n_cells = 50,
                                 n_genes = 80,
                                 seed = 1234)
gcm2 <- data[[1]]
locs2 <- data[[2]]
clusters <- data[[3]]
celltype <- mapvalues(x = clusters,
                      from = unique(clusters),
                      to = letters[seq_len(length(unique(clusters)))])
mdata2 <- data.frame(clusters = clusters,
                     count = colSums(gcm2),
                     celltype = celltype)

mdata <- rbind(mdata1, mdata2)

##### BanksyObject multiple dataset --------------------------------------------

test_that('Construct Banksy without equal number of assays', {
    expect_error(BanksyObject(own.expr = list(gcm1, gcm2),
                              cell.locs = list(locs1)))
})

test_that('Construct Banksy without corresponding assays', {
    locs3 <- locs2
    rownames(locs3) <- paste0('dummy_', seq_len(nrow(locs3)))
    expect_error(BanksyObject(own.expr = list(gcm1, gcm2),
                              cell.locs = list(locs1, locs3)))
})

test_that('Construct Banksy with errant metadata', {
    expect_error(BanksyObject(own.expr = list(gcm1, gcm2),
                              cell.locs = list(locs1, locs2),
                              meta.data = mdata[seq_len(10),,drop=FALSE]))
})


test_that('Construct Banksy', {
    expect_s4_class(BanksyObject(own.expr = list(gcm1, gcm2),
                              cell.locs = list(locs1, locs2),
                              meta.data = mdata), 'BanksyObject')
})

bank <- BanksyObject(own.expr = list(gcm1, gcm2),
                     cell.locs = list(locs1, locs2),
                     meta.data = mdata)

test_that('Show', {
    expect_output(show(bank))
})

test_that('Head', {
    expect_output(head(bank))
})

test_that('Gene filtering', {
    bankTmp <- BanksyObject(own.expr = list(gcm1, gcm2),
                            cell.locs = list(locs1, locs2),
                            genes.filter = 'union')
    expect_equal(nrow(bankTmp@own.expr[[1]]), 80)
    bankTmp <- BanksyObject(own.expr = list(gcm1, gcm2),
                            cell.locs = list(locs1, locs2),
                            genes.filter = 'intersect')
    expect_equal(nrow(bankTmp@own.expr[[2]]), 75)
})

test_that('Slots', {
    expect_equal(length(bank@own.expr), 2)
    expect_equal(length(bank@cell.locs), 2)
})

bank <- NormalizeBanksy(bank)

test_that('Normalization', {
    expect_equal(as.numeric(unlist(sapply(bank@own.expr, colSums2))),
                 rep(100, 50 + 75))
})

bank <- ComputeBanksy(bank)

test_that('Compute Banksy on multiple datasets', {
    expect_s4_class(bank, 'BanksyObject')
})

bank <- NormalizeBanksy(bank)

bank <- ScaleBanksy(bank)

test_that('Scale together', {
    expect_s4_class(ScaleBanksy(bank, separate = FALSE), 'BanksyObject')
})

set.seed(1234)
bankSubset <- SubsetBanksy(bank,
                           cells = sample(bank@meta.data$cell_ID, 50),
                           features = sample(rownames(bank@own.expr[[1]]), 50),
                           dims = sdimx < 800 & sdimy < 1300)

test_that('Subsetting', {
    expect_error(SubsetBanksy(bank, dims = sdimx < 800, dataset = 'dummy'))
    expect_s4_class(SubsetBanksy(bank, dims = sdimx < 800, dataset = 'd1'),
                    'BanksyObject')
})


test_that('Plot Spatial Dimensions', {
    expect_warning(class(plotSpatialDims(bank, by = 'clusters')))
    expect_equal(class(plotSpatialDims(bank, by = 'clusters', dataset = 'd2')),
                 c('gg', 'ggplot'))
})

test_that('Plot own expression heatmap', {
    expect_error(plotHeatmap(bank, assay = 'own.expr', dataset = 'dummy'))
    expect_error(plotHeatmap(bank, assay = 'banksy', dataset = 'dummy'))

    expect_s4_class(plotHeatmap(bank, assay = 'own.expr', dataset = 'd1'),
                    'Heatmap')
})




