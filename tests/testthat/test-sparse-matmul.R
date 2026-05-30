# Test equivalence between sparse matmul and legacy data.table implementations

library(SummarizedExperiment)
library(SpatialExperiment)

data(rings)
spe <- rings

# Precompute neighbors once for reuse across tests
knn_median <- Banksy:::computeNeighbors(
    spatialCoords(spe),
    spatial_mode = "kNN_median", k_geom = 15, verbose = FALSE
)
knn_unif <- Banksy:::computeNeighbors(
    spatialCoords(spe),
    spatial_mode = "kNN_unif", k_geom = 15, verbose = FALSE
)
knn_rank <- Banksy:::computeNeighbors(
    spatialCoords(spe),
    spatial_mode = "kNN_rank", k_geom = 15, verbose = FALSE
)

gcm <- as.matrix(assay(spe, "counts"))

test_that("sparse matmul matches legacy for M=0, center=FALSE", {
    for (knn in list(knn_median, knn_unif, knn_rank)) {
        res_new <- Banksy:::computeHarmonics(
            gcm, knn, M = 0, center = FALSE, verbose = FALSE
        )
        res_old <- Banksy:::.computeHarmonics_legacy(
            gcm, knn, M = 0, center = FALSE, verbose = FALSE
        )
        expect_equal(res_new, res_old,
                     tolerance = sqrt(.Machine$double.eps))
    }
})

test_that("sparse matmul matches legacy for M=1, center=FALSE", {
    for (knn in list(knn_median, knn_unif, knn_rank)) {
        res_new <- Banksy:::computeHarmonics(
            gcm, knn, M = 1, center = FALSE, verbose = FALSE
        )
        res_old <- Banksy:::.computeHarmonics_legacy(
            gcm, knn, M = 1, center = FALSE, verbose = FALSE
        )
        expect_equal(res_new, res_old,
                     tolerance = sqrt(.Machine$double.eps))
    }
})

test_that("sparse matmul matches legacy for M=1, center=TRUE", {
    for (knn in list(knn_median, knn_unif, knn_rank)) {
        res_new <- Banksy:::computeHarmonics(
            gcm, knn, M = 1, center = TRUE, verbose = FALSE
        )
        res_old <- Banksy:::.computeHarmonics_legacy(
            gcm, knn, M = 1, center = TRUE, verbose = FALSE
        )
        expect_equal(res_new, res_old,
                     tolerance = sqrt(.Machine$double.eps))
    }
})

test_that("sparse matmul matches legacy for M=2, center=TRUE", {
    knn_30 <- Banksy:::computeNeighbors(
        spatialCoords(spe),
        spatial_mode = "kNN_median", k_geom = 30, verbose = FALSE
    )
    res_new <- Banksy:::computeHarmonics(
        gcm, knn_30, M = 2, center = TRUE, verbose = FALSE
    )
    res_old <- Banksy:::.computeHarmonics_legacy(
        gcm, knn_30, M = 2, center = TRUE, verbose = FALSE
    )
    expect_equal(res_new, res_old,
                 tolerance = sqrt(.Machine$double.eps))
})

test_that("sparse matmul works with sparse input matrix", {
    gcm_sparse <- Matrix::Matrix(gcm, sparse = TRUE)
    res_sparse <- Banksy:::computeHarmonics(
        gcm_sparse, knn_median, M = 0, center = FALSE, verbose = FALSE
    )
    res_dense <- Banksy:::computeHarmonics(
        gcm, knn_median, M = 0, center = FALSE, verbose = FALSE
    )
    expect_equal(res_sparse, res_dense,
                 tolerance = sqrt(.Machine$double.eps))

    res_sparse <- Banksy:::computeHarmonics(
        gcm_sparse, knn_median, M = 1, center = TRUE, verbose = FALSE
    )
    res_dense <- Banksy:::computeHarmonics(
        gcm, knn_median, M = 1, center = TRUE, verbose = FALSE
    )
    expect_equal(res_sparse, res_dense,
                 tolerance = sqrt(.Machine$double.eps))
})

test_that("sparse matmul with chunking matches unchunked", {
    res_unchunked <- Banksy:::computeHarmonics(
        gcm, knn_median, M = 1, center = TRUE, verbose = FALSE
    )
    res_chunked <- Banksy:::computeHarmonics(
        gcm, knn_median, M = 1, center = TRUE, verbose = FALSE,
        chunk_size = 10
    )
    expect_equal(res_chunked, res_unchunked,
                 tolerance = sqrt(.Machine$double.eps))
})

test_that("full pipeline produces correct dimensions with sparse input", {
    spe_sparse <- spe
    assay(spe_sparse, "counts") <- Matrix::Matrix(
        assay(spe, "counts"), sparse = TRUE
    )
    spe_result <- computeBanksy(spe_sparse, assay_name = "counts",
                                M = 1, k_geom = c(15, 30))
    expect_equal(dim(assay(spe_result, "H0")), dim(spe))
    expect_equal(dim(assay(spe_result, "H1")), dim(spe))
})
