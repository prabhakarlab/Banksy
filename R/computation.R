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
#' \itemize{
#'  \item{kNN_median: k-nearest neighbors with median-scaled Gaussian kernel}
#'  \item{kNN_r: k-nearest neighbors with $1/r$ kernel}
#'  \item{kNN_rn: k-nearest neighbors with $1/r^n$ kernel}
#'  \item{kNN_rank: k-nearest neighbors with rank Gaussian kernel}
#'  \item{kNN_unif: k-nearest neighbors wth uniform kernel}
#'  \item{rNN_gauss: radial nearest neighbors with Gaussian kernel}
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
#'   computing neighborhood.
#' \itemize{
#'  \item{subset of colnames of cell.locs}
#'  \item{all}{Uses all colnames of spatialCoords to compute (default)}
#' }
#' @param center A logical scalar specifying whether to center higher order
#'   harmonics in local neighborhoods.
#' @param verbose A logical scalar specifying verbosity. 
#'
#' @importFrom SummarizedExperiment assay assay<- assayNames
#' @importFrom S4Vectors metadata<-
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

    # Extract expression and locations
    expr <- as.matrix(assay(se, assay_name))
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
        out <- computeHarmonics(expr, knn_df, M, center, verbose = verbose)
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

#' @importFrom matrixStats rowSds
scaler <- function(x) {
    x <- as.matrix(x)
    rm <- rowMeans(x)
    rsd <- rowSds(x)
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


computeHarmonics <- function(gcm, knn_df, M, center, verbose) {
    from <- to <- weight <- phi <- .N <- count <- . <- NULL
    j <- sqrt(as.complex(-1))

    mean_k <- round(mean(knn_df[, .(count = .N), by = from]$count), 1)

    if (verbose) message("Computing harmonic m = ", M)
    if (verbose) message("Using ", mean_k, " neighbors")
    if (center) {
        if (verbose) message("Centering")
        aggr <- knn_df[, abs(
            fscale(gcm[, to, drop = FALSE]) %*% (weight * exp(j * M * phi))
        ), by = from]
    } else {
        aggr <- knn_df[, abs(
            gcm[, to, drop = FALSE] %*% (weight * exp(j * M * phi))
        ), by = from]
    }
    ncm <- matrix(aggr$V1, nrow = nrow(gcm), ncol = ncol(gcm))
    rownames(ncm) <- rownames(gcm)
    colnames(ncm) <- colnames(gcm)
    if (verbose) message("Done")

    return(ncm)
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
