#' Run PCA on a BANKSY embedding.
#'
#' @param se A \code{SpatialExperiment},
#' \code{SingleCellExperiment} or \code{SummarizedExperiment}
#'   object with \code{computeBanksy} ran.
#' @param use_agf A logical vector specifying whether to use the AGF for
#'   computing principal components.
#' @param lambda A numeric vector in \eqn{\in [0,1]} specifying a spatial
#'   weighting parameter. Larger values incorporate more spatial neighborhood
#'   information.
#' @param npcs An integer scalar specifying the number of principal components
#'   to compute.
#' @param assay_name A string scalar specifying the name of the assay used in
#'   \code{computeBanksy}.
#' @param group A string scalar specifying a grouping variable for samples in
#'   \code{se}. This is used to scale the samples in each group separately.
#' @param M Advanced usage. An integer vector specifying the highest azimuthal
#'   Fourier harmonic to use. If specified, overwrites the \code{use_agf}
#'   argument.
#' @param seed Seed for PCA.
#'
#' @importFrom irlba prcomp_irlba
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom matrixStats rowSds
#' @importFrom S4Vectors metadata
#'
#' @return A SpatialExperiment / SingleCellExperiment / SummarizedExperiment
#'   object with PC coordinates in \code{reducedDims(se)}.
#'
#' @export
#'
#' @examples
#' data(rings)
#' spe <- computeBanksy(rings, assay_name = "counts", M = 1, k_geom = c(15, 30))
#' spe <- runBanksyPCA(spe, M = 1, lambda = 0.2, npcs = 20)
#'
runBanksyPCA <- function(se,
                         use_agf = TRUE,
                         lambda = 0.2,
                         npcs = 20,
                         assay_name = NULL,
                         group = NULL,
                         M = NULL,
                         seed = NULL) {
    # Get all combinations of M and lambdas
    param <- expand.grid(lambda, getM(use_agf, M))
    param_names <- sprintf("PCA_M%s_lam%s", param[, 2], param[, 1])

    # Compute PCs
    out <- mapply(function(m, lam) {
        joint <- getBanksyMatrix(
            se,
            M = m,
            lambda = lam,
            assay_name = assay_name,
            scale = TRUE,
            group = group
        )
        verbose.seed(seed)
        joint <- joint[rowSds(joint) != 0, ]
        pca <- prcomp_irlba(t(joint), n = npcs, scale. = FALSE)
        pca_x <- pca$x
        attr(pca_x, "percentVar") <- 100 * pca$sdev^2 / sum(pca$sdev^2)
        pca_x
    }, param[, 2], param[, 1], SIMPLIFY = FALSE)

    # Add PCs to Experiment
    for (i in seq(nrow(param))) reducedDim(se, param_names[i]) <- out[[i]]

    # Log
    metadata(se)$BANKSY_params$lambda <- lambda
    metadata(se)$BANKSY_params$npcs <- npcs
    metadata(se)$BANKSY_params$pca_seed <- seed

    se
}


#' Run UMAP on a BANKSY embedding.
#'
#' @param se A \code{SpatialExperiment},
#' \code{SingleCellExperiment} or \code{SummarizedExperiment}
#'   object with \code{computeBanksy} ran.
#' @param use_agf A logical vector specifying whether to use the AGF for
#'   computing UMAP.
#' @param lambda A numeric vector in \eqn{\in [0,1]} specifying a spatial
#'   weighting parameter. Larger values incorporate more spatial neighborhood
#'   information.
#' @param use_pcs A logical scalar specifying whether to run UMAP on PCs. If
#'   FALSE, runs on the BANKSY matrix.
#' @param npcs An integer scalar specifying the number of principal components
#'   to use if \code{use_pcs} is TRUE.
#' @param dimred A string scalar specifying the name of an existing
#'   dimensionality reduction result to use. Will overwrite \code{use_pcs} if
#'   supplied.
#' @param ndims An integer scalar specifying the number of dimensions to use if
#'   \code{dimred} is supplied.
#' @param assay_name A string scalar specifying the name of the assay used in
#'   \code{computeBanksy}.
#' @param group A string scalar specifying a grouping variable for samples in
#'   \code{se}. This is used to scale the samples in each group separately.
#' @param n_neighbors An integer scalar specifying the number of neighbors to
#'   use for UMAP.
#' @param spread A numeric scalar specifying the effective scale of embedded
#'   points.
#' @param min_dist A numeric scalar specifying the effective min. dist. between
#'   embedded points.
#' @param n_epochs An integer scalar specifying the number of epochs to run
#'   UMAP optimization.
#' @param M Advanced usage. An integer vector specifying the highest azimuthal
#'   Fourier harmonic to use. If specified, overwrites the \code{use_agf}
#'   argument.
#' @param seed Seed for UMAP.
#' @param ... parameters to pass to uwot::umap
#'
#' @importFrom uwot umap
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom S4Vectors metadata<-
#'
#' @return A SpatialExperiment / SingleCellExperiment / SummarizedExperiment
#'   object with UMAP coordinates in \code{reducedDims(se)}.
#'
#' @export
#'
#' @examples
#' data(rings)
#' spe <- computeBanksy(rings, assay_name = "counts", M = 1, k_geom = c(15, 30))
#' spe <- runBanksyPCA(spe, M = 1, lambda = 0.2, npcs = 20)
#' spe <- runBanksyUMAP(spe, M = 1, lambda = 0.2)
#'
runBanksyUMAP <- function(se,
                          use_agf = TRUE,
                          lambda = 0.2,
                          use_pcs = TRUE,
                          npcs = 20,
                          dimred = NULL,
                          ndims = NULL,
                          assay_name = NULL,
                          group = NULL,
                          n_neighbors = 30,
                          spread = 3,
                          min_dist = 0.1,
                          n_epochs = 300,
                          M = NULL,
                          seed = NULL,
                          ...) {
    if (!is.null(dimred)) {
        # Use a custom dimensionality reduction
        ndims <- checkDimred(se, dimred, ndims)
        umap_input <- reducedDim(se, dimred)[, seq(ndims)]
        verbose.seed(seed)
        out <- umap(
            umap_input,
            n_neighbors = n_neighbors,
            min_dist = min_dist,
            n_epochs = n_epochs,
            spread = spread,
            ...
        )
        reducedDim(se, sprintf("UMAP_%s", dimred)) <- out
    } else {
        # Get all combinations of M and lambdas
        param <- expand.grid(lambda, getM(use_agf, M))
        param_names <- sprintf("UMAP_M%s_lam%s", param[, 2], param[, 1])
        if (use_pcs) checkPCA(se, param[, 2], param[, 1], npcs)
        # Compute UMAPs
        out <- mapply(function(m, lam) {
            if (!use_pcs) {
                umap_input <- t(getBanksyMatrix(
                    se,
                    assay_name = assay_name,
                    lambda = lam,
                    M = m,
                    scale = TRUE,
                    group = group
                ))
            } else {
                umap_input <- reducedDim(se, sprintf("PCA_M%s_lam%s", m, lam))
                umap_input <- umap_input[, seq(npcs)]
            }
            verbose.seed(seed)
            umap(
                umap_input,
                n_neighbors = n_neighbors,
                min_dist = min_dist,
                n_epochs = n_epochs,
                spread = spread,
                ...
            )
        }, param[, 2], param[, 1], SIMPLIFY = FALSE)

        # Add UMAPs to Experiment
        for (i in seq(nrow(param))) reducedDim(se, param_names[i]) <- out[[i]]
    }

    # Log
    metadata(se)$BANKSY_params$n_neighbors <- n_neighbors
    metadata(se)$BANKSY_params$min_dist <- min_dist
    metadata(se)$BANKSY_params$spread <- spread
    metadata(se)$BANKSY_params$umap_seed <- seed

    se
}
