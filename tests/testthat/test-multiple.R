# 
# BanksyObject single dataset --------------------------------------------------

library(testthat)
library(matrixStats)
library(plyr)

set.seed(1000)
d1 <- Banksy::simulateDataset(n_cells = c(400,100,100),
                              n_deg = c(20,20,20), de_dropout = c(1,1,1))
d2 <- Banksy::simulateDataset(n_cells = c(200,200,200),
                              n_deg = c(20,20,20), de_dropout = c(1,1,1))

gcm <- list(d1$gcm[sample(seq_len(100), 60),],
            d2$gcm[sample(seq_len(100), 70),])
locs <- list(d1$locs, d2$locs)
meta <- rbind(d1$meta, d2$meta)

## Construction

test_that('Construct Banksy with different number of inputs', {
    expect_error(BanksyObject(own.expr = list(d1$gcm, d2$gcm),
                              meta.data = list(d1$locs)))
    expect_error(BanksyObject(own.expr = gcm, cell.locs = locs,
                              meta.data = meta[seq_len(100),]))
})

test_that('Construct Banksy', {

    bank <- BanksyObject(own.expr = gcm, cell.locs = locs,
                         genes.filter = 'intersect')
    expect_equal(as.numeric(sapply(bank@own.expr, nrow)), rep(41,2))

    bank <- BanksyObject(own.expr = gcm, cell.locs = locs,
                         genes.filter = 'intersect', min.cells.expressed = 5)
    expect_equal(as.numeric(sapply(bank@own.expr, nrow)), rep(41,2))

    expect_s4_class(BanksyObject(own.expr = gcm, cell.locs = locs,
                         genes.filter = 'union'), 'BanksyObject')

})

bank <- BanksyObject(own.expr = gcm,
                     cell.locs = locs,
                     meta.data = meta)

test_that('Show', {
    expect_output(show(bank))
})

test_that('Head', {
    expect_output(head(bank))
})

# Normalization own

bank <- NormalizeBanksy(bank)

test_that('Normalization for own', {

    expect_s4_class(NormalizeBanksy(bank, assay = 'own'), 'BanksyObject')

})

# Compute neighbors

bank <- ComputeBanksy(bank)

test_that('Compute neighbors', {
    expect_s4_class(ComputeBanksy(bank), 'BanksyObject')
    expect_s4_class(ComputeBanksy(bank, M=1, k_geom=c(15,30)), 'BanksyObject')
})

# Normalization nbr

test_that('Normalization for neighbors', {
    expect_s4_class(NormalizeBanksy(bank, assay = 'nbr'), 'BanksyObject')
})

# Scale

bank <- ScaleBanksy(bank)

test_that('Scaling', {
    expect_s4_class(ScaleBanksy(bank), 'BanksyObject')
    expect_s4_class(ScaleBanksy(bank, separate = FALSE), 'BanksyObject')
})

test_that('Plot spatial', {
    expect_error(plotSpatial(bank))
    expect_equal(class(plotSpatial(bank, dataset = 'd1')), c('gg', 'ggplot'))
    expect_equal(class(plotSpatial(bank, dataset = 'd1', by = 'Label',
                                   type = 'discrete')), c('gg', 'ggplot'))
    expect_equal(class(plotSpatial(bank, dataset = 'd1', by = 'gene_64',
                                   type = 'continuous')), c('gg', 'ggplot'))
})

test_that('Plot heatmap', {

    expect_error(plotHeatmap(bank, dataset = 'hello'))

    expect_s4_class(plotHeatmap(bank, max.cols = 100), 'Heatmap')
    expect_s4_class(plotHeatmap(bank, dataset = 'd1', max.cols = 100),
                    'Heatmap')
    
    expect_s4_class(plotHeatmap(bank, M=1, max.cols = 100, assay = 'banksy',
                                lambda = 0.25), 'Heatmap')
    expect_s4_class(plotHeatmap(bank, M=1, assay = 'banksy', dataset = 'd1',
                                lambda = 0.25, max.cols = 100), 'Heatmap')
})

bank <- RunBanksyPCA(bank, lambda = 0.25, npcs = 20)
bank <- RunBanksyUMAP(bank, lambda = 0.25, npcs = 20)

# Subset
test_that('Subsetting', {

    expect_error(SubsetBanksy(bank, dataset = 'hello'))

    bankSubset <- SubsetBanksy(bank,
        cells = sample(bank@meta.data$cell_ID, 500),
        features = sample(rownames(bank@own.expr[[1]]), 20),
        dims = sdimx > 0)

    all_cellnames <- as.character(unlist(sapply(bankSubset@own.expr, colnames)))

    expect_equal(colnames(bankSubset@own.expr[[1]]),
                 rownames(bankSubset@cell.locs[[1]]))
    expect_equal(colnames(bankSubset@own.expr[[2]]),
                 rownames(bankSubset@cell.locs[[2]]))
    expect_equal(colnames(bankSubset@nbr.expr[[1]]),
                 rownames(bankSubset@cell.locs[[1]]))
    expect_equal(colnames(bankSubset@nbr.expr[[2]]),
                 rownames(bankSubset@cell.locs[[2]]))
    expect_equal(all_cellnames, bankSubset@meta.data$cell_ID)
    expect_equal(all_cellnames, rownames(bankSubset@reduction$pca_M1_lam0.25$x))
    expect_equal(all_cellnames, rownames(bankSubset@reduction$umap_M1_lam0.25))

})

# Split
test_that('Split', {
    expect_warning(SplitBanksy(bank, by = 'Label'))
})


