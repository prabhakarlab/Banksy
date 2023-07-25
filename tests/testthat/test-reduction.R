
# Test funcs. in reduction.R

library(SummarizedExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)

data(rings)
spe <- rings
spe <- computeBanksy(spe, assay_name='counts')

test_that('runBanksyPCA gives message when seeded', {
    expect_message(runBanksyPCA(spe, seed = 1000))
})

test_that('runBanksyPCA gives expected output', {
    spe <- runBanksyPCA(spe)
    expect_equal(dim(reducedDim(spe)), c(ncol(spe), 20))
    expect_in('percentVar', names(attributes(reducedDim(spe))))
})


test_that('runBanksyUMAP without PCs computed', {
    expect_error(runBanksyUMAP(spe, use_pcs = TRUE))
})

test_that('runBanksyUMAP requesting too many PCs', {
    spe <- runBanksyPCA(spe)
    expect_error(runBanksyUMAP(spe, use_pcs = TRUE, npcs = 50))
})


test_that('runBanksyUMAP gives expected output', {
    spe <- runBanksyPCA(spe)
    spe <- runBanksyUMAP(spe, use_pcs = TRUE)
    expect_equal(dim(reducedDim(spe, 'UMAP_M1_lam0.2')), c(ncol(spe), 2))

    spe <- runBanksyUMAP(spe, use_pcs = FALSE)
    expect_equal(dim(reducedDim(spe, 'UMAP_M1_lam0.2')), c(ncol(spe), 2))
})

test_that('runBanksyUMAP gives message when seed', {
    spe <- runBanksyPCA(spe)
    expect_message(runBanksyUMAP(spe, use_pcs = TRUE, seed = 1000))
})

