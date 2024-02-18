# Test funcs. in reduction.R

library(SummarizedExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)

data(rings)
spe <- rings
spe <- computeBanksy(spe, assay_name = "counts", compute_agf = TRUE)

test_that("runBanksyPCA gives message when seeded", {
    expect_message(runBanksyPCA(spe, use_agf = TRUE, seed = 1000))
})

test_that("runBanksyPCA gives expected output", {
    spe <- runBanksyPCA(spe, use_agf = TRUE)
    expect_equal(dim(reducedDim(spe)), c(ncol(spe), 20))
    expect_in("percentVar", names(attributes(reducedDim(spe))))
})


test_that("runBanksyUMAP without PCs computed", {
    expect_error(runBanksyUMAP(spe, use_agf = TRUE, use_pcs = TRUE))
})

test_that("runBanksyUMAP requesting too many PCs", {
    spe <- runBanksyPCA(spe, use_agf = TRUE)
    expect_error(runBanksyUMAP(spe, use_agf = TRUE, use_pcs = TRUE, npcs = 50))
})


test_that("runBanksyUMAP gives expected output", {
    spe <- runBanksyPCA(spe, use_agf = TRUE)
    spe <- runBanksyUMAP(spe, use_agf = TRUE, use_pcs = TRUE)
    expect_equal(dim(reducedDim(spe, "UMAP_M1_lam0.2")), c(ncol(spe), 2))

    spe <- runBanksyUMAP(spe, use_agf = TRUE, use_pcs = FALSE)
    expect_equal(dim(reducedDim(spe, "UMAP_M1_lam0.2")), c(ncol(spe), 2))
})

test_that("runBanksyUMAP gives message when seed", {
    spe <- runBanksyPCA(spe, use_agf = TRUE)
    expect_message(runBanksyUMAP(spe, use_agf = TRUE, 
                                 use_pcs = TRUE, seed = 1000))
})
