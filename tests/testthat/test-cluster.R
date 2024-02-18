library(SummarizedExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)

set.seed(1000)
spe <- simulateDataset(rate = 2)
spe <- computeBanksy(spe, assay_name = "counts", compute_agf = TRUE)
spe <- runBanksyPCA(spe, use_agf = TRUE, seed = 1000)

test_that("clusterBanksy with invalid algo", {
    expect_error(clusterBanksy(spe, algo = "?"))
})

test_that("clusterBanksy with leiden specifies correct args", {
    expect_error(clusterBanksy(spe, k_neighbors = NULL))
    expect_error(clusterBanksy(spe, resolution = NULL))
})

test_that("clusterBanksy gives message when seeded", {
    expect_message(clusterBanksy(spe, use_agf = TRUE, seed = 1000))
})

test_that("clusterBanksy gives expected output", {
    spe1 <- clusterBanksy(spe, use_agf = TRUE, use_pcs = TRUE, seed = 1000)
    expect_equal(length(clusterNames(spe1)), 1)
    spe2 <- clusterBanksy(spe, use_agf = TRUE, use_pcs = FALSE, seed = 1000)
    expect_equal(length(clusterNames(spe2)), 1)
})

spe <- clusterBanksy(spe, use_agf = TRUE, resolution = 1, seed = 1000)

test_that("connectClusters with invalid map_to", {
    expect_error(connectClusters(spe, map_to = "?"))
})

test_that("connectClusters gives expected output", {
    cnm <- tail(colnames(colData(spe)), n = 1)

    before <- sum(spe$cluster == colData(spe)[, cnm])

    spe_connect <- connectClusters(spe)
    after <- sum(spe_connect$cluster == colData(spe_connect)[, cnm])
    expect_gte(after, before)

    spe_connect <- connectClusters(spe, map_to = cnm)
    after <- sum(spe_connect$cluster == colData(spe_connect)[, cnm])
    expect_gte(after, before)
})

test_that("smoothLabels with invalid cluster_names", {
    expect_error(smoothLabels(spe, cluster_names = "?"))
})

test_that("smoothLabels with invalid cluster_names", {
    data(rings)
    rings$cluster <- NULL
    expect_error(smoothLabels(rings))
})

test_that("smoothLabels with gives expected output", {
    cnm <- tail(colnames(colData(spe)), n = 1)
    spe <- smoothLabels(spe, cnm)
    expect_true(any(grepl("smooth", colnames(colData(spe)))))
})

test_that("compareClusters with no clusters", {
    spe$cluster <- NULL
    expect_error(compareClusters(spe))
})

test_that("compareClusters gives expected output", {
    cnm <- tail(colnames(colData(spe)), n = 1)
    expect_true(is.numeric(
        compareClusters(spe, func = c("ARI", "NMI"))
    ))
})
