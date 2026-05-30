#' Compute the component neighborhood matrices for the BANKSY matrix.
#' 
#' @details
#' Given an expression matrix (as specified by \code{assay_name}), this function
#' computes the mean neighborhood matrix (\code{H0}) and optionally, the 
#' azimuthal Gabor filter (AGF) matrix (\code{H1}). The number of neighbors 
#' used to define the spatial neighborhood is given by \code{k_geom}. 
#' Different kernels may be used to compute the neighborhood features, 
#' specified by \code{spatial_mode}.   
#'
#' @param se A \code{SpatialExperiment},
#' \code{SingleCellExperiment} or \code{SummarizedExperiment}
#'   object. If not a SpatialExperiment object, argument \code{coord_names}
#'   must be provided.
#' @param assay_name A string scalar specifying the name of the assay to use.
#' @param coord_names A string vector specifying the names in \code{colData}
#'   corresponding to spatial coordinates.
#' @param compute_agf A logical scalar specifying whether to compute the AGF.
#' @param k_geom An integer scalar specifying the number of neighbors to use.
#'   Values \eqn{\in [15,30]} work well.
#' @param spatial_mode A string scalar specifying the kernel for neighborhood
#'   computation (default: kNN_median).
#' \describe{
#'  \item{kNN_median}{k-nearest neighbors with median-scaled Gaussian kernel}
#'  \item{kNN_r}{k-nearest neighbors with $1/r$ kernel}
#'  \item{kNN_rn}{k-nearest neighbors with $1/r^n$ kernel}
#'  \item{kNN_rank}{k-nearest neighbors with rank Gaussian kernel}
#'  \item{kNN_unif}{k-nearest neighbors with uniform kernel}
#'  \item{rNN_gauss}{radial nearest neighbors with Gaussian kernel}
#' }
#' @param n A numeric scalar specifying the exponent of radius (for kNN_rn).
#' @param sigma A numeric scalar specifying the std. dev. of Gaussian kernel
#'   (for rNN_gauss).
#' @param alpha A numeric scalar specifying the radius used: larger alphas give
#'   smaller radii (for rNN_gauss).
#' @param k_spatial An integer scalar specifying the initial number of neighbors
#'   to use (for rNN_gauss)
#' @param M Advanced usage. A integer scalar specifying the highest azimuthal
#'   Fourier harmonic to compute. If specified, overwrites the \code{use_agf}
#'   argument.
#' @param sample_size An integer scalar number of neighbors to sample from the
#'   neighborhood.
#' @param sample_renorm A logical scalar specifying whether to renormalize the
#'   neighbor weights to 1.
#' @param seed An integer scalar specifying seed for sampling the neighborhood.
#' @param dimensions A character vector specifying the dimensions to use when
#'   computing neighborhood. One of:
#' \describe{
#'  \item{a subset of \code{spatialCoords} colnames}{use the specified dimensions}
#'  \item{all}{use all \code{spatialCoords} colnames (default)}
#' }
#' @param center A logical scalar specifying whether to center higher order
#'   harmonics in local neighborhoods.
#' @param chunk_size A integer scalar specifying the number of rows / genes of 
#'   the neighborhood cell matrix to compute. Must be strictly less than floor 
#'   of 2e31-1 / number of cells, though this may still give rise to negative
#'   length vector errors. For safety, set this to less than 
#'   2e31-1 * row_limit_factor / number of cells. This is automatically 
#'   computed based on the latter, but can be specified.
#' @param parallel A logical scalar specifying whether to compute chunks in 
#'   parallel using bplapply. Not implemented for Windows.
#' @param num_cores A integer scalar specifying the number of cores to use 
#'   if parallel is TRUE.
#' @param row_limit_factor A numeric scalar specifying the safety factor 
#'   applied to the maximum vector length (2^31-1) when computing chunk sizes. 
#'   Accounts for data.table overhead during grouping operations. 
#'   Default is 0.75.
#' @param verbose A logical scalar specifying verbosity. 
#'
#' @importFrom SummarizedExperiment assay assay<- assayNames
#' @importFrom S4Vectors metadata<-
#' @importFrom Matrix sparseMatrix
#'
#' @return A SpatialExperiment / SingleCellExperiment / SummarizedExperiment
#'   object with neighborhood matrices added.
#'
#' @export
#'
#' @examples
#' data(rings)
#' spe <- computeBanksy(rings, assay_name = "counts", M = 1, k_geom = c(15, 30))
#'
computeBanksy <- function(se,
                          assay_name,
                          coord_names = NULL,
                          compute_agf = FALSE,
                          k_geom = 15,
                          spatial_mode = c("kNN_median", "kNN_r", 
                                           "kNN_rn", "kNN_rank", "kNN_unif", 
                                           "rNN_gauss"),
                          n = 2,
                          sigma = 1.5,
                          alpha = 0.05,
                          k_spatial = 100L,
                          M = NULL,
                          sample_size = NULL,
                          sample_renorm = TRUE,
                          seed = NULL,
                          dimensions = "all",
                          center = TRUE,
                          chunk_size = NULL,
                          parallel = FALSE,
                          num_cores = NULL,
                          row_limit_factor = 0.75,
                          verbose = TRUE) {
    
    # Check args
    if (missing(assay_name)) {
        stop("Provide argument `assay_name`. One of ", assayNames(se))
    }
    spatial_mode <- match.arg(spatial_mode)
    checkComputeBanksy(as.list(environment()))

    # Compute Ms
    M <- seq(0, max(getM(compute_agf, M)))
    if (length(k_geom) == 1) k_geom <- rep_len(k_geom, max(M) + 1)
    # Compute kgeoms
    if (length(k_geom) != length(M)) {
        stop(
            "Specify either one k_geom only or sufficient k_geoms for each of ",
            length(M), " harmonics."
        )
    }

    # Extract expression and locations (keep sparse if input is sparse)
    expr <- assay(se, assay_name)
    locs <- getLocs(se, coord_names)
    if (ncol(locs) == 0) stop("No spatial coordinates found")

    # Compute neighbors
    knn_list <- lapply(k_geom, function(kg) {
        computeNeighbors(locs,
            spatial_mode = spatial_mode, k_geom = kg,
            n = n, sigma = sigma, alpha = alpha, k_spatial = k_spatial,
            sample_size = sample_size, sample_renorm = sample_renorm,
            seed = seed, dimensions = dimensions, verbose = verbose
        )
    })

    # Compute harmonics with different k_geoms
    center <- c(FALSE, rep(center, length(M) - 1))
    har <- Map(function(knn_df, M, center) {
        out <- computeHarmonics(expr, knn_df, M, center, 
                                verbose = verbose, chunk_size = chunk_size, 
                                parallel = parallel, num_cores = num_cores, 
                                row_limit_factor = row_limit_factor)
        rownames(out) <- rownames(expr)
        out
    }, knn_list, M, center)
    names(har) <- paste0("H", M)

    # Add harmonics to Experiment
    for (i in seq(length(M))) assay(se, paste0("H", i - 1)) <- har[[i]]

    # Log
    metadata(se)$BANKSY_params <- list(
        assay_name = assay_name,
        M = M,
        k_geom = k_geom,
        spatial_mode = spatial_mode
    )

    se
}

# Argument checks for computeBanksy
checkComputeBanksy <- function(params) {
    stopifnot("compute_agf should be a logical scalar" = 
                  is.logical(as.logical(params$compute_agf)) & 
                  length(params$compute_agf) == 1)
    stopifnot("n should be a numeric scalar" = 
                  is.numeric(params$n) & length(params$n) == 1)
    stopifnot("sigma should be a numeric scalar" = 
                  is.numeric(params$sigma) & length(params$sigma) == 1)
    stopifnot("k_spatial should be an integer scalar" = 
                  is.integer(as.integer(params$k_spatial)) & 
                  length(params$k_spatial) == 1)
    stopifnot("sample_size should be an integer scalar" = 
                  (is.integer(as.integer(params$sample_size)) & 
                       length(params$sample_size) == 1) | 
                  is.null(params$sample_size))
    stopifnot("sample_renorm should be a logical scalar" = 
                  is.logical(params$sample_renorm) & 
                  length(params$sample_renorm) == 1)
    stopifnot("center should be a logical scalar" = 
                  is.logical(params$center) & length(params$center) == 1)
    stopifnot("verbose should be a logical scalar" = 
                  is.logical(params$verbose) & length(params$verbose) == 1)
}

#' Builds the BANKSY matrix from neighborhood matrices. 
#' 
#' @details
#' After computation of the neighborhood matrices 
#' (see \link[Banksy]{computeBanksy}), this function builds the BANKSY matrix by 
#' concatenating the original expression matrix with the neighborhood matrices,
#' and scales each matrix by an appropriate weight as determined by 
#' \code{lambda}. The weights of the own expression matrix, mean neighborhood
#'  matrix and azimuthal Gabor filter are given by \eqn{\sqrt{1-\lambda}}, 
#'  \eqn{\sqrt{\lambda/\mu}} and \eqn{\sqrt{\lambda/2\mu}} respectively, where 
#'  \eqn{\mu=1.5}. In the case where the AGF is not computed, the weights for
#'  the own and mean neighborhood expression matrix simplify to 
#'  \eqn{\sqrt{1-\lambda}} and \eqn{\sqrt{\lambda}} respectively. 
#'
#' @param se A \code{SpatialExperiment},
#' \code{SingleCellExperiment} or \code{SummarizedExperiment}
#'   object with \code{computeBanksy} ran.
#' @param M A integer scalar specifying the highest azimuthal
#'   Fourier harmonic to compute.
#' @param lambda A numeric vector in \eqn{\in [0,1]} specifying a spatial
#'   weighting parameter. Larger values (e.g. 0.8) incorporate more spatial
#'   neighborhood and find spatial domains, while smaller values (e.g. 0.2)
#'   perform spatial cell-typing.
#' @param assay_name A string scalar specifying the name of the assay used in
#'   \code{computeBanksy}.
#' @param scale A logical scalar specifying whether to scale the features to
#'   zero mean and unit standard deviation. This is performed before
#'   multiplying the assays by their corresponding lambda weighting factors.
#' @param group A string scalar specifying a grouping variable for samples in
#'   \code{se}. This is used to scale the samples in each group separately.
#' @param verbose A logical scalar specifying verbosity.
#'
#' @importFrom SummarizedExperiment assays assayNames
#' @importFrom S4Vectors metadata
#'
#' @return BANKSY matrix.
#'
#' @export
#'
#' @examples
#' data(rings)
#' spe <- computeBanksy(rings, assay_name = "counts", M = 1, k_geom = c(15, 30))
#' banksyMatrix <- getBanksyMatrix(spe, M = 1, lambda = 0.2)
#'
getBanksyMatrix <- function(se,
                            M,
                            lambda,
                            assay_name = NULL,
                            scale = FALSE,
                            group = NULL,
                            verbose = TRUE) {
    M <- seq(0, M)
    anames <- assayNames(se)
    if (is.null(assay_name)) {
        assay_name <- metadata(se)$BANKSY_params$assay_name
    }
    banksy_names <- c(assay_name, paste0("H", M))
    not_found <- which(!(banksy_names %in% anames))

    if (length(not_found) > 0) {
        err_msg <- paste0(banksy_names[not_found], collaspe = "; ")
        stop(
            "The following assays are missing: ", err_msg, "Run computeBanksy
            (see ?computeBanksy for more details). 
            "
        )
    }

    banksy_assays <- assays(se)[banksy_names]
    lambdas <- getLambdas(lambda, n_harmonics = length(banksy_assays) - 1)

    # Scale features
    if (scale) {
        if (is.null(group)) {
            # Single dataset case
            banksy_assays <- lapply(banksy_assays, scaler)
        } else {
            # Multi dataset case
            if (!(group %in% colnames(colData(se)))) {
                stop(
                    "Invalid group variable ", group
                )
            }
            # Scale each group separately
            groups <- colData(se)[, group]
            ugroups <- unique(groups)
            for (i in seq(length(banksy_assays))) {
                if (verbose) {
                    if (i == 1) {
                        message("Scaling ", assay_name)
                    } else {
                        message("Scaling harmonic m = ", i - 2)
                    }
                }
                curr_assay <- as.matrix(banksy_assays[[i]])
                for (curr_group in ugroups) {
                    if (verbose) message("Group: ", curr_group)
                    curr_group_id <- which(curr_group == groups)
                    curr_assay[, curr_group_id] <-
                        scaler(curr_assay[, curr_group_id])
                }
                banksy_assays[[i]] <- curr_assay
            }
        }
    }

    # Multiply BANKSY matrices by their factors
    banksy_assays <- Map(function(lam, mat) lam * mat, lambdas, banksy_assays)

    # Rownames
    nfeat <- dim(se)[1]
    suffix <- rep(paste0("_H", M), each = nfeat)

    # Concat
    joint <- do.call(rbind, banksy_assays)
    rownames(joint)[seq(nfeat + 1, nrow(joint))] <- paste0(
        rownames(joint)[seq(nfeat + 1, nrow(joint))], suffix
    )

    joint
}

#' @importFrom MatrixGenerics rowSds
scaler <- function(x) {
    x <- as.matrix(x)
    rm <- rowMeans(x)
    rsd <- MatrixGenerics::rowSds(x)
    x <- (x - rm) / rsd
    x[is.nan(x)] <- 0
    return(x)
}

getLambdas <- function(lambda, n_harmonics) {
    weights <- lambda * (2^-seq(0, n_harmonics - 1))
    weights <- weights / sum(2^-seq(0, n_harmonics - 1))
    lam <- c(1 - sum(weights), weights)
    sqrt(lam)
}

getM <- function(use_agf, M) {
    if (is.null(M)) M <- as.numeric(use_agf)
    sort(M)
}

#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment colData
getLocs <- function(se, coord_names) {
    # Extract coordinates
    if (inherits(se, "SpatialExperiment")) {
        locs <- spatialCoords(se)
    } else {
        if (is.null(coord_names)) {
            stop(
                "Specify coord_names corresponding to spatial coordinates."
            )
        }
        if (!all(coord_names %in% colnames(colData(se)))) {
            stop(
                "Specify valid coord_names."
            )
        }
        locs <- colData(se)[, coord_names]
    }
    as.matrix(locs)
}

#' @importFrom data.table `:=`
computeNeighbors <- function(locs,
                             spatial_mode = "kNN_median", k_geom = 15, n = 2,
                             sigma = 1.5, alpha = 0.05, k_spatial = 100,
                             sample_size = NULL, sample_renorm = TRUE,
                             seed = NULL, dimensions = "all", verbose = FALSE) {
    from <- to <- phi <- NULL
    locs <- as.matrix(locs)
    kernelRadius <- sqrt(-ncol(locs) * log(alpha))

    if (verbose) message("Computing neighbors...")
    knnDF <- switch(spatial_mode, 
        'rNN_gauss' = withRNNgauss(
            locs = locs, sigma = sigma, k_spatial = k_spatial,
            kernelRadius = kernelRadius, verbose = verbose),
        'kNN_rank' = withKNNrank(
            locs = locs, k_geom = k_geom, verbose = verbose),
        'kNN_r' = withKNNr(locs = locs, k_geom = k_geom, verbose = verbose),
        'kNN_rn' = withKNNrn(locs, k_geom = k_geom, n = n, verbose = verbose),
        'kNN_unif' = withKNNunif(
            locs = locs, k_geom = k_geom, verbose = verbose),
        'kNN_median' = withKNNmedian(
            locs = locs, k_geom = k_geom, verbose = verbose)
                    )
    
    knnDF[, phi := getPhi(locs, from, to), by = from][]
    if (!is.null(sample_size)) {
        if (verbose) message("Subsampling to ", sample_size, " neighbors")
        knnDF <- subsampler(knnDF,
            sample_size = sample_size,
            sample_renorm = sample_renorm,
            seed = seed
        )
    }
    if (verbose) message("Done")
    return(knnDF)
}


computeHarmonics <- function(gcm, knn_df, M, center, verbose,
                             chunk_size = NULL, parallel = FALSE,
                             num_cores = NULL, row_limit_factor = 0.75) {

    if (parallel && verbose) {
        message("Note: parallel is not used by the sparse matmul implementation; ",
                "BLAS-level parallelism is used automatically")
    }

    n_genes <- nrow(gcm)
    n_cells <- ncol(gcm)

    # Extract knn_df columns. Filter out sentinel rows (to=0) that
    # rNN_gauss inserts for isolated cells — these have weight=0 and
    # would crash sparseMatrix which requires positive indices.
    from_idx <- knn_df[["from"]]
    to_idx <- knn_df[["to"]]
    w <- knn_df[["weight"]]
    p <- knn_df[["phi"]]
    valid <- to_idx > 0L & from_idx > 0L
    if (!all(valid)) {
        from_idx <- from_idx[valid]
        to_idx <- to_idx[valid]
        w <- w[valid]
        p <- p[valid]
    }

    if (verbose) {
        mean_k <- round(mean(tabulate(from_idx, nbins = n_cells)), 1)
        message("Computing harmonic m = ", M, " with ", mean_k, " neighbors")
    }

    # Determine chunk size: balance cache efficiency vs chunk overhead.
    # Target: gcm chunk fits in L3 cache (~30 MB) for good locality
    # during sparse matmul. Also cap total intermediates at ~2 GB.
    if (is.null(chunk_size)) {
        cache_chunk <- as.integer(floor(30e6 / (8 * as.double(n_cells))))
        mem_chunk <- as.integer(floor(2e9 / (8 * 4 * as.double(n_cells))))
        max_chunk <- max(100L, min(cache_chunk, mem_chunk, n_genes))
    } else {
        max_chunk <- min(as.integer(chunk_size), n_genes)
    }
    num_chunks <- ceiling(n_genes / max_chunk)
    if (verbose && num_chunks > 1) {
        message("Processing in ", num_chunks, " chunks of max ",
                max_chunk, " genes")
    }

    # Build sparse weight matrices (n_cells x n_cells, k_geom nnz per col)
    if (M == 0) {
        # Real case: W[to, from] = weight
        W <- sparseMatrix(i = to_idx, j = from_idx, x = w,
                          dims = c(n_cells, n_cells))

        ncm <- matrix(0, nrow = n_genes, ncol = n_cells)
        for (ch in seq_len(num_chunks)) {
            ri <- .chunkIdx(ch, max_chunk, n_genes)
            ncm[ri, ] <- as.matrix(abs(gcm[ri, , drop = FALSE] %*% W))
        }
    } else {
        # Complex case: split into real and imaginary parts
        W_re <- sparseMatrix(i = to_idx, j = from_idx,
                             x = w * cos(M * p),
                             dims = c(n_cells, n_cells))
        W_im <- sparseMatrix(i = to_idx, j = from_idx,
                             x = w * sin(M * p),
                             dims = c(n_cells, n_cells))

        # Centering: H_m = |gcm %*% W_m - (gcm %*% U) * s|
        # where U is uniform neighbor weights and s = colSums(W_m)
        if (center) {
            k_per <- tabulate(from_idx, nbins = n_cells)
            k_per[k_per == 0] <- 1L
            U <- sparseMatrix(i = to_idx, j = from_idx,
                              x = 1 / k_per[from_idx],
                              dims = c(n_cells, n_cells))
            s_re <- Matrix::colSums(W_re)
            s_im <- Matrix::colSums(W_im)
            if (verbose) message("Centering")
        }

        ncm <- matrix(0, nrow = n_genes, ncol = n_cells)
        for (ch in seq_len(num_chunks)) {
            ri <- .chunkIdx(ch, max_chunk, n_genes)
            gcm_chunk <- gcm[ri, , drop = FALSE]

            re <- as.matrix(gcm_chunk %*% W_re)
            im <- as.matrix(gcm_chunk %*% W_im)

            if (center) {
                mn <- as.matrix(gcm_chunk %*% U)
                re <- re - sweep(mn, 2, s_re, `*`)
                im <- im - sweep(mn, 2, s_im, `*`)
            }

            ncm[ri, ] <- sqrt(re^2 + im^2)
        }
    }

    rownames(ncm) <- rownames(gcm)
    colnames(ncm) <- colnames(gcm)

    if (verbose) message("Done")
    return(ncm)
}

# Legacy data.table implementation, kept for equivalence testing
.computeHarmonics_legacy <- function(gcm, knn_df, M, center, verbose,
                                     chunk_size = NULL, row_limit_factor = 0.75) {
    from <- to <- weight <- phi <- .N <- count <- . <- NULL
    j <- sqrt(as.complex(-1))

    total_rows <- as.double(nrow(gcm)) * ncol(gcm)
    max_rows <- (2^31 - 1) * row_limit_factor
    if (total_rows > max_rows || !is.null(chunk_size)) {
        max_chunk_size <- floor(max_rows / ncol(gcm))
        if (!is.null(chunk_size)) {
            if (chunk_size > max_chunk_size) {
                stop('Specified chunk_size too large. Must be smaller than ',
                     floor(max_rows / ncol(gcm)))
            }
            max_chunk_size <- chunk_size
        }
    } else {
        max_chunk_size <- nrow(gcm)
    }

    num_chunks <- ceiling(nrow(gcm) / max_chunk_size)

    process_chunk <- function(chunk) {
        start_idx <- (chunk - 1) * max_chunk_size + 1
        end_idx <- min(chunk * max_chunk_size, nrow(gcm))
        gcm_chunk <- gcm[start_idx:end_idx, , drop = FALSE]

        if (center) {
            chunk_aggr <- knn_df[, abs(
                fscale(gcm_chunk[, to, drop = FALSE]) %*%
                    (weight * exp(j * M * phi))
            ), by = from]
        } else {
            chunk_aggr <- knn_df[, abs(
                gcm_chunk[, to, drop = FALSE] %*%
                    (weight * exp(j * M * phi))
            ), by = from]
        }
        chunk_aggr$V1
    }

    chunk_results <- lapply(seq_len(num_chunks), process_chunk)

    ncm <- matrix(0, nrow = nrow(gcm), ncol = ncol(gcm))
    for (i in seq_along(chunk_results)) {
        start_idx <- (i - 1) * max_chunk_size + 1
        end_idx <- min(i * max_chunk_size, nrow(gcm))
        ncm[start_idx:end_idx, ] <- matrix(chunk_results[[i]],
                                           nrow = end_idx - start_idx + 1,
                                           ncol = ncol(gcm))
    }

    rownames(ncm) <- rownames(gcm)
    colnames(ncm) <- colnames(gcm)
    return(ncm)
}

.chunkIdx <- function(chunk, max_chunk, n_total) {
    start <- (chunk - 1L) * max_chunk + 1L
    end <- min(chunk * max_chunk, n_total)
    start:end
}

# Build sparse weight matrix W from knn_df (M=0 only).
# W[j, i] = w_ij (weight of neighbor j for cell i).
.buildWeightMatrix <- function(knn_df, n_cells) {
    sparseMatrix(
        i = knn_df[["to"]], j = knn_df[["from"]],
        x = knn_df[["weight"]], dims = c(n_cells, n_cells)
    )
}

# Compute row-wise mean and sd of (gcm %*% W) without forming the full product.
# Processes in gene-row chunks to control memory.
.computeH0ScalingParams <- function(gcm, W, chunk_size = 100L) {
    n_genes <- nrow(gcm)
    n_cells <- ncol(gcm)
    mu <- numeric(n_genes)
    ss <- numeric(n_genes)

    for (ch_start in seq(1L, n_genes, by = chunk_size)) {
        ch_end <- min(ch_start + chunk_size - 1L, n_genes)
        ri <- ch_start:ch_end
        chunk <- as.matrix(gcm[ri, , drop = FALSE] %*% W)
        mu[ri] <- rowMeans(chunk)
        ss[ri] <- rowSums(chunk * chunk)
    }

    # Sample sd (n-1 denominator) to match Seurat::FastRowScale
    sd <- sqrt(pmax(n_cells / (n_cells - 1) * (ss / n_cells - mu * mu), 0))
    sd[sd == 0] <- 1
    list(mu = mu, sd = sd)
}


getPhi <- function(locs, from, to) {
    out <- sweep(locs[to, , drop = FALSE], 2, locs[from, , drop = FALSE], "-")
    phi <- atan2(out[, 2, drop = FALSE], out[, 1, drop = FALSE])
    phi + as.integer(phi < 0) * 2 * pi
}


fscale <- function(x) {
    rm <- rowMeans(x)
    x <- (x - rm)
    return(x)
}


#' @importFrom dbscan kNN
#' @importFrom data.table data.table setDT  setnames `:=` rbindlist
#' @importFrom stats dnorm median
withRNNgauss <- function(locs, sigma, k_spatial, kernelRadius, verbose) {
    if (verbose) message("Spatial mode is rNN gaussian")
    if (verbose) message("Parameters: sigma=", sigma, ", k_spatial=", k_spatial)

    tryCatch(
        {
            knn <- dbscan::kNN(x = locs, k = k_spatial)
        },
        error = function(cond) {
            stop("Not enough neighbours at k_spatial=", k_spatial)
        }
    )

    medianDist <- median(as.vector(knn$dist[, 1]))
    knnDF <- data.table(
        from = rep(seq_len(nrow(knn$id)), k_spatial),
        to = as.vector(knn$id),
        weight = dnorm(as.vector(knn$dist),
            mean = 0,
            sd = medianDist * sigma
        ),
        distance = as.vector(knn$dist)
    )

    distance <- norm.weight <- weight <- from <- to <- NULL
    knnDF <- knnDF[distance < sigma * kernelRadius * medianDist, ]
    setDT(knnDF)[, norm.weight := weight / sum(weight), by = from]
    knnDF <- knnDF[, -3, with = FALSE]

    ## Create dummy entries for filtered out cells
    iso <- setdiff(seq_len(nrow(locs)), unique(knnDF$from))
    isomat <- c(rep(iso, 2), rep(0, 2 * length(iso)))
    isomat <- data.table(matrix(isomat, ncol = ncol(knnDF)))
    knnDF <- rbindlist(list(knnDF, isomat), use.names = FALSE)
    knnDF <- knnDF[order(from, to)]
    setnames(knnDF, "norm.weight", "weight")

    return(knnDF)
}


#' @importFrom dbscan kNN
#' @importFrom data.table data.table setnames `:=`
withKNNrank <- function(locs, k_geom, verbose) {
    if (verbose) message("Spatial mode is kNN_rank")
    if (verbose) message("Parameters: k_geom=", k_geom)

    tryCatch(
        {
            knn <- dbscan::kNN(x = locs, k = k_geom)
            unnormWt <- exp(-seq(1, k_geom, 1)^2 / (2 * (k_geom / 1.5)^2))
            normWt <- unnormWt / sum(unnormWt)

            weightMatrix <-
                t(matrix(normWt, nrow = k_geom, ncol = nrow(knn$id)))
            knnDF <- data.table(
                from = rep(seq_len(nrow(knn$id)), k_geom),
                to = as.vector(knn$id),
                weight = as.vector(weightMatrix),
                distance = as.vector(knn$dist)
            )
        },
        error = function(cond) {
            stop("Not enough neighbours at k_geom=", k_geom)
        }
    )

    return(knnDF)
}


#' @importFrom dbscan kNN
#' @importFrom data.table data.table setnames `:=`
withKNNr <- function(locs, k_geom, verbose) {
    if (verbose) message("Spatial mode is kNN_r")
    if (verbose) message("Parameters: k_geom=", k_geom)

    tryCatch(
        {
            knn <- dbscan::kNN(x = locs, k = k_geom)
        },
        error = function(cond) {
            stop("Not enough neighbours at k_geom=", k_geom)
        }
    )

    norm.weight <- weight <- from <- NULL
    knnDF <- data.table(
        from = rep(seq_len(nrow(knn$id)), k_geom),
        to = as.vector(knn$id),
        weight = 1 / as.vector(knn$dist),
        distance = as.vector(knn$dist)
    )
    knnDF[, norm.weight := weight / sum(weight), by = from]
    knnDF <- knnDF[, -3, with = FALSE]
    setnames(knnDF, "norm.weight", "weight")

    return(knnDF)
}


#' @importFrom dbscan kNN
#' @importFrom data.table data.table setnames `:=`
withKNNrn <- function(locs, k_geom, n, verbose) {
    if (verbose) message("Spatial mode is kNN_rn")
    if (verbose) message("Parameters: k_geom=", k_geom, ", n=", n)

    tryCatch(
        {
            knn <- dbscan::kNN(x = locs, k = k_geom)
        },
        error = function(cond) {
            stop("Not enough neighbours at k_geom=", k_geom)
        }
    )

    norm.weight <- weight <- from <- NULL
    knnDF <- data.table(
        from = rep(seq_len(nrow(knn$id)), k_geom),
        to = as.vector(knn$id),
        weight = 1 / (as.vector(knn$dist)^n),
        distance = as.vector(knn$dist)
    )
    knnDF[, norm.weight := weight / sum(weight), by = from]
    knnDF <- knnDF[, -3, with = FALSE]
    setnames(knnDF, "norm.weight", "weight")

    return(knnDF)
}


#' @importFrom dbscan kNN
#' @importFrom data.table data.table setnames `:=`
withKNNunif <- function(locs, k_geom, verbose) {
    if (verbose) message("Spatial mode is kNN_unif")
    if (verbose) message("Parameters: k_geom = ", k_geom)

    tryCatch(
        {
            knn <- dbscan::kNN(x = locs, k = k_geom)
        },
        error = function(cond) {
            stop("Not enough neighbours at k_geom=", k_geom)
        }
    )

    norm.weight <- weight <- from <- NULL
    knnDF <- data.table(
        from = rep(seq_len(nrow(knn$id)), k_geom),
        to = as.vector(knn$id),
        weight = 1,
        distance = as.vector(knn$dist)
    )
    knnDF[, norm.weight := weight / sum(weight), by = from]
    knnDF <- knnDF[, -3, with = FALSE]
    setnames(knnDF, "norm.weight", "weight")

    return(knnDF)
}


#' @importFrom dbscan kNN
#' @importFrom data.table data.table setnames `:=`
#' @importFrom stats median
withKNNmedian <- function(locs, k_geom, verbose) {
    if (verbose) message("Spatial mode is kNN_median")
    if (verbose) message("Parameters: k_geom=", k_geom)

    tryCatch(
        {
            knn <- dbscan::kNN(x = locs, k = k_geom)
        },
        error = function(cond) {
            stop("Not enough neighbours at k_geom=", k_geom)
        }
    )

    norm.weight <- weight <- from <- distance <- NULL
    knnDF <- data.table(
        from = rep(seq_len(nrow(knn$id)), k_geom),
        to = as.vector(knn$id),
        distance = as.vector(knn$dist)
    )
    knnDF[, weight := exp(-distance^2 / median(distance)^2), by = from]
    knnDF[, norm.weight := weight / sum(weight), by = from]
    knnDF <- knnDF[, -c(3, 4), with = FALSE]
    setnames(knnDF, "norm.weight", "weight")

    return(knnDF)
}

#' @importFrom data.table data.table `:=` .SD .N
subsampler <- function(knnDF,
                       sample_size = NULL,
                       sample_renorm = TRUE,
                       seed = NULL) {
    verbose.seed(seed)
    from <- weight <- NULL
    x <- knnDF[,
        .SD[sample(.N, min(sample_size, .N), replace = FALSE)],
        by = from
    ]
    if (sample_renorm) x[, weight := weight / sum(weight), by = from]
    data.table(x)
}
