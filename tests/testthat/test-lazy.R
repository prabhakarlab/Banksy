# Test funcs. in lazy.R

library(SummarizedExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)

data(rings)
# Order cells by group so the per-group reference keeps the same column order
spe <- rings[, order(rings$cluster)]

npcs <- 10L
lambda <- 0.2
k_geom <- 15L
# Groups are a quarter of the object, so use a smaller neighborhood
k_geom_group <- 5L

# Smallest singular value of the leading-k bases' crossprod: 1 means the two
# embeddings span the same subspace. PC sign and the order of near-degenerate
# PCs are arbitrary, so compare subspaces, not embeddings.
subspace_overlap <- function(e1, e2, k = k_cmp) {
    q1 <- qr.Q(qr(scale(e1[, seq_len(k), drop = FALSE], scale = FALSE)))
    q2 <- qr.Q(qr(scale(e2[, seq_len(k), drop = FALSE], scale = FALSE)))
    min(svd(crossprod(q1, q2))$d)
}

# Drop the trailing PCs: they sit at the edge of the Krylov subspace and are
# least converged, and rings has signal in only ~n_rings directions, so beyond
# that the spectrum is near-degenerate and the span is ill-determined.
k_cmp <- npcs - 2L

# Dense reference: materialize the BANKSY matrix, then irlba
dense_pca <- function(se, ...) {
    suppressMessages(runBanksyPCA(
        se, assay_name = "counts", lazy = FALSE, lambda = lambda,
        npcs = npcs, seed = 1000, verbose = FALSE, ...
    ))
}

lazy_pca <- function(se, ...) {
    suppressMessages(runBanksyPCA(
        se, assay_name = "counts", lazy = TRUE, lambda = lambda,
        npcs = npcs, verbose = FALSE, ...
    ))
}

pca_name <- sprintf("PCA_M0_lam%s", lambda)

test_that("lazy PCA matches the dense path", {
    ref <- computeBanksy(spe, assay_name = "counts", compute_agf = FALSE,
                         k_geom = k_geom, verbose = FALSE)
    e_dense <- reducedDim(dense_pca(ref), pca_name)

    for (backend in c("cpp", "r")) {
        e_lazy <- reducedDim(
            lazy_pca(spe, k_geom = k_geom, pca_backend = backend), pca_name)
        expect_equal(dim(e_lazy), dim(e_dense))
        expect_gt(subspace_overlap(e_dense, e_lazy), 0.99)
    }
})

test_that("lazy PCA matches the dense path per group", {
    # The dense path computes one kNN over all cells, so build the per-group
    # reference group by group. Overwrite H0 in place rather than cbind-ing:
    # keeps the column order and avoids SpatialExperiment's cbind method.
    ref <- computeBanksy(spe, assay_name = "counts", compute_agf = FALSE,
                         k_geom = k_geom_group, verbose = FALSE)
    h0 <- assay(ref, "H0")
    for (cid in split(seq_len(ncol(spe)), spe$cluster)) {
        sub <- computeBanksy(spe[, cid], assay_name = "counts",
                             compute_agf = FALSE, k_geom = k_geom_group,
                             verbose = FALSE)
        h0[, cid] <- as.matrix(assay(sub, "H0"))
    }
    assay(ref, "H0") <- h0

    # split_scale = FALSE: per-group neighborhoods, one global z-scaling
    e_dense <- reducedDim(dense_pca(ref), pca_name)
    e_lazy <- reducedDim(
        lazy_pca(spe, group = "cluster", split_scale = FALSE,
                 k_geom = k_geom_group), pca_name)
    expect_gt(subspace_overlap(e_dense, e_lazy), 0.99)

    # split_scale = TRUE: features z-scaled within each group
    e_dense <- reducedDim(dense_pca(ref, group = "cluster"), pca_name)
    e_lazy <- reducedDim(
        lazy_pca(spe, group = "cluster", split_scale = TRUE,
                 k_geom = k_geom_group), pca_name)
    expect_gt(subspace_overlap(e_dense, e_lazy), 0.99)
})

test_that("lazy PCA gives expected output", {
    se <- lazy_pca(spe, k_geom = k_geom)
    expect_equal(dim(reducedDim(se, pca_name)), c(ncol(spe), npcs))
    expect_in("percentVar", names(attributes(reducedDim(se, pca_name))))
    expect_equal(metadata(se)$BANKSY_params$npcs, npcs)
})

test_that("lazy PCA does not require computeBanksy", {
    expect_false("H0" %in% assayNames(spe))
    expect_error(dense_pca(spe))
    se <- lazy_pca(spe, k_geom = k_geom)
    expect_true(pca_name %in% reducedDimNames(se))
})

test_that("lazy PCA is deterministic", {
    # Both backends start from a fixed vector, so repeated runs are identical
    # without a seed and the caller's RNG stream is untouched
    for (backend in c("cpp", "r")) {
        e1 <- reducedDim(lazy_pca(spe, k_geom = k_geom,
                                  pca_backend = backend), pca_name)
        e2 <- reducedDim(lazy_pca(spe, k_geom = k_geom,
                                  pca_backend = backend), pca_name)
        expect_identical(e1, e2)
    }
    set.seed(1); before <- .Random.seed
    invisible(lazy_pca(spe, k_geom = k_geom))
    expect_identical(before, .Random.seed)
})

test_that("lazy PCA rejects the AGF", {
    expect_error(lazy_pca(spe, use_agf = TRUE, k_geom = k_geom))
    expect_error(lazy_pca(spe, M = 1, k_geom = k_geom))
})

test_that("lazy PCA rejects empty groups", {
    # unique() keeps NA but no cell matches it, giving a zero-sized group
    spe_na <- spe
    spe_na$cluster[1] <- NA
    expect_error(
        lazy_pca(spe_na, group = "cluster", k_geom = k_geom_group),
        "no cells"
    )
})
