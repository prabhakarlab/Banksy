# Test funcs. in computation.R

library(SummarizedExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)

data(rings)
spe <- rings
sce <- SingleCellExperiment(spe)
assay(sce) <- NULL
assay(sce, "counts") <- assay(spe, "counts")
colData(sce) <- cbind(colData(spe), spatialCoords(spe))

test_that("computeBanksy without assay_name", {
    expect_error(
        computeBanksy(spe)
    )
})

test_that("computeBanksy with M and k_geom mismatch", {
    expect_error(
        computeBanksy(spe, assay_name = "counts", M = 1, k_geom = c(15, 30, 60))
    )
})

test_that("computeBanksy with invalid assay_name", {
    expect_error(
        computeBanksy(spe, assay_name = "hello")
    )
})

test_that("computeBanksy with invalid assay_name", {
    expect_error(
        computeBanksy(spe, assay_name = "?")
    )
})

test_that("computeBanksy on SPE without spatial coords", {
    spatialCoords(spe) <- NULL
    expect_error(
        computeBanksy(spe, assay_name = "counts")
    )
})

test_that("computeBanksy on SCE without coord_names", {
    expect_error(
        computeBanksy(sce, assay_name = "counts")
    )
})

test_that("computeBanksy on SCE with invalid coord_names", {
    expect_error(
        computeBanksy(sce, assay_name = "counts", coord_names = c("x", "?"))
    )
})

test_that("computeNeighbors with invalid spatial mode", {
    expect_error(
        computeBanksy(spe, assay_name = "counts", spatial_mode = "?")
    )
})

test_that("computeBanksy with spatial mode kNN_median", {
    spe <- computeBanksy(spe, assay_name = "counts", spatial_mode = "kNN_median")
    expect_equal(dim(assay(spe, "H0")), dim(spe))
    expect_equal(dim(assay(spe, "H1")), dim(spe))

    expect_equal(metadata(spe)$BANKSY_params$assay_name, "counts")
    expect_equal(metadata(spe)$BANKSY_params$M, c(0, 1))
    expect_equal(metadata(spe)$BANKSY_params$k_geom, c(15, 15))
    expect_equal(metadata(spe)$BANKSY_params$spatial_mode, "kNN_median")

    sce <- computeBanksy(sce,
        assay_name = "counts", spatial_mode = "kNN_median",
        coord_names = c("x", "y")
    )
    expect_equal(dim(assay(sce, "H0")), dim(sce))
    expect_equal(dim(assay(sce, "H1")), dim(sce))
})

test_that("computeBanksy with subsampling", {
    spe <- computeBanksy(spe,
        assay_name = "counts", spatial_mode = "kNN_median",
        k_geom = 30, sample_size = 15, sample_renorm = TRUE,
        seed = 1000
    )
    expect_equal(dim(assay(spe, "H0")), dim(spe))
    expect_equal(dim(assay(spe, "H1")), dim(spe))
})

test_that("computeBanksy with spatial mode kNN_median with XL k_geom", {
    expect_error(
        computeBanksy(spe,
            assay_name = "counts",
            spatial_mode = "kNN_median", k_geom = ncol(spe) + 1
        )
    )
})

test_that("computeBanksy with spatial_mode kNN_unif", {
    spe <- computeBanksy(spe, assay_name = "counts", spatial_mode = "kNN_unif")
    expect_equal(dim(assay(spe, "H0")), dim(spe))
    expect_equal(dim(assay(spe, "H1")), dim(spe))
})

test_that("computeBanksy with spatial mode kNN_unif with XL k_geom", {
    expect_error(
        computeBanksy(spe,
            assay_name = "counts",
            spatial_mode = "kNN_unif", k_geom = ncol(spe) + 1
        )
    )
})

test_that("computeBanksy with spatial_mode kNN_rn", {
    spe <- computeBanksy(spe, assay_name = "counts", spatial_mode = "kNN_rn")
    expect_equal(dim(assay(spe, "H0")), dim(spe))
    expect_equal(dim(assay(spe, "H1")), dim(spe))
})

test_that("computeBanksy with spatial mode kNN_rn with XL k_geom", {
    expect_error(
        computeBanksy(spe,
            assay_name = "counts",
            spatial_mode = "kNN_rn", k_geom = ncol(spe) + 1
        )
    )
})

test_that("computeBanksy with spatial_mode kNN_r", {
    spe <- computeBanksy(spe, assay_name = "counts", spatial_mode = "kNN_r")
    expect_equal(dim(assay(spe, "H0")), dim(spe))
    expect_equal(dim(assay(spe, "H1")), dim(spe))
})

test_that("computeBanksy with spatial mode kNN_r with XL k_geom", {
    expect_error(
        computeBanksy(spe,
            assay_name = "counts",
            spatial_mode = "kNN_r", k_geom = ncol(spe) + 1
        )
    )
})

test_that("computeBanksy with spatial_mode kNN_rank", {
    spe <- computeBanksy(spe, assay_name = "counts", spatial_mode = "kNN_rank")
    expect_equal(dim(assay(spe, "H0")), dim(spe))
    expect_equal(dim(assay(spe, "H1")), dim(spe))
})

test_that("computeBanksy with spatial mode kNN_rank with XL k_geom", {
    expect_error(
        computeBanksy(spe,
            assay_name = "counts",
            spatial_mode = "kNN_rank", k_geom = ncol(spe) + 1
        )
    )
})

test_that("computeBanksy with spatial_mode rNN_gauss", {
    spe <- computeBanksy(spe, assay_name = "counts", spatial_mode = "rNN_gauss")
    expect_equal(dim(assay(spe, "H0")), dim(spe))
    expect_equal(dim(assay(spe, "H1")), dim(spe))
})

test_that("computeBanksy with spatial mode rNN_gauss with XL k_spatial", {
    expect_error(
        computeBanksy(spe,
            assay_name = "counts",
            spatial_mode = "rNN_gauss", k_spatial = ncol(spe) + 1
        )
    )
})

test_that("getBanksyMatrix without computeBanksy", {
    expect_error(
        getBanksyMatrix(spe, M = 1, lambda = 0.2)
    )
})

test_that("getBanksyMatrix", {
    spe <- computeBanksy(spe, assay_name = "counts")

    dim_B0 <- c(nrow(spe) * 2, ncol(spe))
    expect_equal(dim(getBanksyMatrix(spe, M = 0, lambda = 0.2)), dim_B0)

    dim_B1 <- c(nrow(spe) * 3, ncol(spe))
    expect_equal(dim(getBanksyMatrix(spe, M = 1, lambda = 0.2)), dim_B1)

    dim_B0 <- c(nrow(spe) * 2, ncol(spe))
    expect_equal(dim(getBanksyMatrix(spe, M = 0, lambda = 0.2, scale = TRUE)), dim_B0)

    dim_B1 <- c(nrow(spe) * 3, ncol(spe))
    expect_equal(dim(getBanksyMatrix(spe, M = 1, lambda = 0.2, scale = TRUE)), dim_B1)

    expect_equal(
        as.numeric(rowMeans(getBanksyMatrix(
            spe, M = 1, lambda = 0.2, scale = TRUE))),
        rep(0, nrow(spe) * 3)
    )

    mat <- getBanksyMatrix(spe, M = 0, lambda = 0)
    expect_equal(unique(tail(rowMeans(mat), nrow(spe))), 0)

    mat <- getBanksyMatrix(spe, M = 0, lambda = 1)
    expect_equal(unique(head(rowMeans(mat), nrow(spe))), 0)
})

test_that("getBanksyMatrix multi-sample", {
    spe <- computeBanksy(spe, "counts")

    expect_error(getBanksyMatrix(
        spe,
        M = 1, lambda = 0.2, scale = TRUE, group = "?"
    ))

    expect_equal(
        as.numeric(rowMeans(getBanksyMatrix(
            spe,
            M = 1, lambda = 0.2, scale = TRUE, group = "cluster"
        ))),
        rep(0, nrow(spe) * 3)
    )

    rmeans <- as.numeric(rowMeans(getBanksyMatrix(
        spe,
        M = 1, lambda = 0.2, scale = TRUE, group = "cluster"
    )))

    expect_equal(
        rmeans[seq(nrow(spe))],
        rep(0, nrow(spe))
    )
    expect_equal(
        rmeans[seq(nrow(spe) + 1, 2 * nrow(spe))],
        rep(0, nrow(spe))
    )
})
