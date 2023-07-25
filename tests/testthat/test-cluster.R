
library(SummarizedExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)

set.seed(1000)
spe <- simulateDataset(rate=2)
spe <- computeBanksy(spe, assay_name='counts')
spe <- runBanksyPCA(spe, seed=1000)

test_that('clusterBanksy with invalid algo', {
    expect_error(clusterBanksy(spe, algo='?'))
})

test_that('clusterBanksy with leiden specifies correct args', {
    expect_error(clusterBanksy(spe, k_neighbors = NULL))
    expect_error(clusterBanksy(spe, resolution = NULL))
})

test_that('clusterBanksy gives message when seeded', {
    expect_message(clusterBanksy(spe, seed=1000))
})

test_that('clusterBanksy gives expected output', {
    expect_equal(ncol(colData(clusterBanksy(spe, use_pcs = TRUE))), 3)
    expect_equal(ncol(colData(clusterBanksy(spe, use_pcs = FALSE))), 3)
})

spe <- clusterBanksy(spe, resolution=1, seed=1000)

test_that('connectClusters with invalid map_to', {
    expect_error(connectClusters(spe, map_to = '?'))
})

test_that('connectClusters gives expected output', {
    
    cnm <- tail(colnames(colData(spe)), n=1)
    
    before <- sum(spe$cluster ==  colData(spe)[,cnm])

    spe_connect <- connectClusters(spe)
    after <- sum(spe_connect$cluster == colData(spe_connect)[,cnm])
    expect_gte(after, before)

    spe_connect <- connectClusters(spe, map_to = cnm)
    after <- sum(spe_connect$cluster == colData(spe_connect)[,cnm])
    expect_gte(after, before)
    
    
})

table(spe$cluster, spe$clust_M1_lam0.2_k50_res1)

test_that('smoothLabels with invalid cluster_names', {
    expect_error(smoothLabels(spe, cluster_names='?'))
})

test_that('smoothLabels with invalid cluster_names', {
    data(rings)
    rings$cluster <- NULL
    expect_error(smoothLabels(rings))
})

test_that('smoothLabels with gives expected output', {
    cnm <- tail(colnames(colData(spe)), n=1)
    spe <- smoothLabels(spe, cnm)
    expect_true(any(grepl('smooth', colnames(colData(spe)))))
})

test_that('compareClusters with no clusters', {
    spe$cluster <- NULL
    expect_error(compareClusters(spe))
})

test_that('compareClusters gives expected output', {
    cnm <- tail(colnames(colData(spe)), n=1)
    expect_true(is.numeric(
        compareClusters(spe, func=c('ARI', 'NMI'))
    ))
})
