#' Run PCA on a BANKSY matrix.
#' 
#' @details
#' This function runs PCA on the BANKSY matrix 
#' (see \link[Banksy]{getBanksyMatrix}) with features scaled to zero mean and 
#' unit standard deviation.   
#'
#' @param se A \code{SpatialExperiment},
#' \code{SingleCellExperiment} or \code{SummarizedExperiment}
#'   object with \code{computeBanksy} ran.
#' @param use_agf A logical vector specifying whether to use the AGF for
#'   computing principal components.
#' @param lambda A numeric vector in \eqn{\in [0,1]} specifying a spatial
#'   weighting parameter. Larger values (e.g. 0.8) incorporate more spatial
#'   neighborhood and find spatial domains, while smaller values (e.g. 0.2)
#'   perform spatial cell-typing. 
#' @param npcs An integer scalar specifying the number of principal components
#'   to compute.
#' @param assay_name A string scalar specifying the name of the assay used in
#'   \code{computeBanksy}.
#' @param scale A logical scalar specifying whether to scale features before
#'   PCA. Defaults to TRUE.
#' @param group A string scalar specifying a grouping variable for samples in
#'   \code{se}. This is used to scale the samples in each group separately.
#' @param M Advanced usage. An integer vector specifying the highest azimuthal
#'   Fourier harmonic to use. If specified, overwrites the \code{use_agf}
#'   argument.
#' @param seed Seed for PCA. If not specified, no seed is set. 
#'
#' @importFrom irlba prcomp_irlba irlba
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom MatrixGenerics rowSds
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
                         use_agf = FALSE,
                         lambda = 0.2,
                         npcs = 20L,
                         assay_name = NULL,
                         scale = TRUE,
                         group = NULL,
                         M = NULL,
                         seed = NULL) {
    
    # Check parameters
    checkBanksyPCA(as.list(environment()))
    
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
            scale = scale,
            group = group
        )
        verbose.seed(seed)
        joint <- joint[MatrixGenerics::rowSds(joint) != 0, ]
        pca <- irlba::irlba(
            Matrix::t(joint), nv = npcs, nu = npcs, 
            scale. = FALSE, center = TRUE)
        pca_x <- pca$u %*% diag(pca$d)
        percentVar <- 100 * pca$d^2/sum(pca$d^2)
        colnames(pca_x) <- paste0("PC", seq_len(npcs))
        attr(pca_x, "percentVar") <- percentVar
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

# Argument checks for runBanksyPCA
checkBanksyPCA <- function(params) {
    stopifnot("use_agf should be a logical vector" = 
                  is.logical(as.logical(params$compute_agf)))
    stopifnot("lambda should be a numeric vector with each entry in [0,1]" = 
                  is.numeric(params$lambda) & 
                  max(params$lambda) <= 1 & 
                  min(params$lambda >= 0))
    stopifnot("npcs should be an integer scalar" = 
                  is.integer(as.integer(params$npcs)) & 
                  length(params$npcs) == 1)
}

#' Run UMAP on a BANKSY embedding.
#'
#' @details
#' This function runs UMAP on the principal components computed on the 
#' BANKSY matrix.   
#'
#' @param se A \code{SpatialExperiment},
#' \code{SingleCellExperiment} or \code{SummarizedExperiment}
#'   object with \code{computeBanksy} ran.
#' @param use_agf A logical vector specifying whether to use the AGF for
#'   computing UMAP.
#' @param lambda A numeric vector in \eqn{\in [0,1]} specifying a spatial
#'   weighting parameter. Larger values (e.g. 0.8) incorporate more spatial
#'   neighborhood and find spatial domains, while smaller values (e.g. 0.2)
#'   perform spatial cell-typing.
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
#' @param scale A logical scalar specifying whether to scale features before
#'   UMAP. Only used when use_pcs is FALSE. Defaults to TRUE.
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
#' @param seed Seed for UMAP. If not specified, no seed is set. 
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
                          use_agf = FALSE,
                          lambda = 0.2,
                          use_pcs = TRUE,
                          npcs = 20L,
                          dimred = NULL,
                          ndims = NULL,
                          assay_name = NULL,
                          scale = TRUE,
                          group = NULL,
                          n_neighbors = 30L,
                          spread = 3,
                          min_dist = 0.1,
                          n_epochs = 300L,
                          M = NULL,
                          seed = NULL,
                          ...) {
    # Check params
    checkBanksyUMAP(as.list(environment()))
    
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
                    scale = scale,
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

# Argument checks for runBanksyUMAP
checkBanksyUMAP <- function(params) {
    stopifnot("use_agf should be a logical vector" = 
                  is.logical(as.logical(params$compute_agf)))
    stopifnot("lambda should be a numeric vector with each entry in [0,1]" = 
                  is.numeric(params$lambda) & 
                  max(params$lambda) <= 1 & 
                  min(params$lambda >= 0))
    stopifnot("use_pcs should be an integer scalar" = 
                  is.logical(params$use_pcs) & length(params$use_pcs) == 1)
    stopifnot("npcs should be an integer scalar" = 
                  is.integer(as.integer(params$npcs)) & 
                  length(params$npcs) == 1)
    stopifnot("n_neighbors should be an integer scalar" = 
                  is.integer(as.integer(params$n_neighbors)) & 
                  length(params$n_neighbors) == 1)
    stopifnot("spread should be a numeric scalar" = 
                  is.numeric(params$spread) & length(params$spread) == 1)
    stopifnot("min_dist should be a numeric scalar" = 
                  is.numeric(params$min_dist) & length(params$min_dist) == 1)
    stopifnot("n_epochs should be an integer scalar" =
                  is.integer(as.integer(params$n_epochs)) &
                  length(params$n_epochs) == 1)
}

#' Run lazy PCA on a BANKSY matrix.
#'
#' @details
#' This function computes PCA on the BANKSY matrix without materializing the
#' full matrix in memory. Instead, it uses an implicit operator that applies
#' the BANKSY transform on the fly during the iterative PCA solver. This
#' enables analysis of very large datasets (millions of cells) that would
#' otherwise exceed available memory.
#'
#' Unlike \code{\link{runBanksyPCA}}, this function does not require
#' \code{\link{computeBanksy}} to be run first. It computes the kNN graph
#' and weight matrices internally.
#'
#' Currently only supported for M=0 (no AGF).
#'
#' @param se A \code{SpatialExperiment},
#'   \code{SingleCellExperiment} or \code{SummarizedExperiment}
#'   object.
#' @param lambda A numeric scalar in \eqn{[0,1]} specifying the spatial
#'   weighting parameter.
#' @param npcs An integer scalar specifying the number of principal components
#'   to compute (default 50).
#' @param assay_name A string scalar specifying the name of the assay to use.
#' @param coord_names A string vector specifying the names in \code{colData}
#'   corresponding to spatial coordinates.
#' @param k_geom An integer scalar specifying the number of neighbors to use.
#' @param spatial_mode A string scalar specifying the kernel for neighborhood
#'   computation (default: kNN_median).
#' @param group A string scalar specifying a grouping variable for samples in
#'   \code{se}. This is used for per-group kNN computation and optionally
#'   per-group scaling.
#' @param split_scale A logical scalar specifying whether to scale features
#'   per group. Only used when \code{group} is not NULL.
#' @param scale_max A numeric scalar specifying the maximum absolute z-score
#'   for clipping (default 10).
#' @param pca_backend A string scalar specifying the PCA backend.
#'   \code{"cpp"} (default) uses C++ irlba for lower memory and faster runtime.
#'   \code{"r"} uses R's irlba package.
#' @param name A string scalar specifying the name for the dimensionality
#'   reduction in \code{reducedDims(se)} (default: \code{"lazyPCA_M0_lamX"}).
#' @param verbose A logical scalar specifying verbosity.
#' @param ... Additional arguments passed to \code{computeNeighbors}.
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom S4Vectors metadata metadata<-
#'
#' @return A SpatialExperiment / SingleCellExperiment / SummarizedExperiment
#'   object with PC coordinates in \code{reducedDims(se)}.
#'
#' @export
#'
#' @examples
#' data(rings)
#' spe <- runBanksyLazyPCA(spe, assay_name = "counts", lambda = 0.2,
#'                          k_geom = 15, npcs = 20)
#'
runBanksyLazyPCA <- function(se,
                              lambda = 0.2,
                              npcs = 50L,
                              assay_name = NULL,
                              coord_names = NULL,
                              k_geom = 15L,
                              spatial_mode = c("kNN_median", "kNN_r",
                                               "kNN_rn", "kNN_rank",
                                               "kNN_unif", "rNN_gauss"),
                              group = NULL,
                              split_scale = TRUE,
                              scale_max = 10,
                              pca_backend = c("cpp", "r"),
                              name = NULL,
                              verbose = TRUE,
                              ...) {

    # Validate
    if (lambda < 0 || lambda > 1) stop('lambda must be between 0 and 1')
    spatial_mode <- match.arg(spatial_mode)
    pca_backend <- match.arg(pca_backend)
    if (is.null(assay_name)) {
        assay_name <- SummarizedExperiment::assayNames(se)[1]
        if (verbose) message('Using assay: ', assay_name)
    }

    # Default reduction name
    if (is.null(name)) name <- sprintf("lazyPCA_M0_lam%s", lambda)

    # Extract expression (genes x cells)
    data_own <- assay(se, assay_name)
    if (!inherits(data_own, 'dgCMatrix'))
        data_own <- as(data_own, 'dgCMatrix')

    # Get spatial coordinates
    locs <- getLocs(se, coord_names)

    do_split_scale <- !is.null(group) && split_scale

    if (!is.null(group)) {
        # Per-group path
        groups_vec <- colData(se)[, group]
        ugroups <- unique(groups_vec)
        group_idx <- lapply(ugroups, function(g) which(g == groups_vec))

        if (verbose) {
            message('Computing per-group neighbors')
            for (gr in seq_along(group_idx))
                message('  ', ugroups[gr], ': ', length(group_idx[[gr]]), ' cells')
        }

        # Per-group kNN on unstaggered coordinates
        group_locs <- lapply(group_idx, function(cid) locs[cid, , drop = FALSE])
        group_knn <- lapply(group_locs, function(loc_slice)
            lapply(k_geom, function(kg)
                computeNeighbors(loc_slice, spatial_mode = spatial_mode,
                                 k_geom = kg, verbose = FALSE, ...)))
        rm(group_locs)

        result <- .banksy_lazy_pca_core(
            data_own, knn_list = NULL,
            group_knn = group_knn, group_idx = group_idx,
            lambda = lambda, npcs = npcs,
            split_scale = do_split_scale,
            scale_max = scale_max, pca_backend = pca_backend,
            verbose = verbose)
        rm(data_own, group_knn)
    } else {
        # Single-matrix path
        knn_list <- lapply(k_geom, function(kg)
            computeNeighbors(locs, spatial_mode = spatial_mode,
                             k_geom = kg, verbose = verbose, ...))

        result <- .banksy_lazy_pca_core(
            data_own, knn_list = knn_list,
            lambda = lambda, npcs = npcs,
            split_scale = FALSE,
            scale_max = scale_max, pca_backend = pca_backend,
            verbose = verbose)
        rm(data_own, knn_list)
    }

    # Store in reducedDims: matrix with percentVar attribute
    pca_x <- result$embeddings
    percentVar <- 100 * result$stdev^2 / result$total_var *
                  (max(1, ncol(pca_x)) - 1)
    # Simpler: use stdev^2 directly since total_var = sum(d^2)
    percentVar <- 100 * result$stdev^2 * max(1, nrow(pca_x) - 1) /
                  result$total_var
    attr(pca_x, "percentVar") <- percentVar

    reducedDim(se, name) <- pca_x

    # Log
    metadata(se)$BANKSY_params$lambda <- lambda
    metadata(se)$BANKSY_params$npcs <- npcs
    metadata(se)$BANKSY_params$pca_backend <- pca_backend

    if (verbose) message('Done. Access reduction with reducedDim(se, "', name, '")')
    se
}
