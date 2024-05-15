#' Perform clustering in BANKSY's neighborhood-augmented feature space.
#' 
#' @details
#' This function performs clustering on the principal components computed on 
#' the BANKSY matrix, i.e., the BANKSY embedding. The PCA corresponding to the 
#' parameters \code{use_agf} and \code{lambda} must have been computed with  
#' \link[Banksy]{runBanksyPCA}. Clustering may also be performed directly on the 
#' BANKSY matrix with \code{use_pcs} set to \code{FALSE} (this is not 
#' recommended). 
#' 
#' Four clustering algorithms are implemented. 
#' \itemize{
#'  \item{leiden: Leiden graph-based clustering. The arguments 
#'  \code{k_neighbors} and \code{resolution} should be specified.}
#'  \item{louvain: Louvain graph-based clustering. The arguments 
#'  \code{k_neighbors} and \code{resolution} should be specified.}
#'  \item{kmeans: kmeans clustering. The argument \code{kmeans.centers} should 
#'  be specified.}
#'  \item{mclust: Gaussian mixture model-based clustering. The argument 
#'  \code{mclust.G} should be specified.}
#' }
#' 
#' By default, no seed is set for clustering. If a seed is specified, the same
#' seed is used for clustering across the input parameters. 
#' 
#'
#' @param se A \code{SpatialExperiment},
#' \code{SingleCellExperiment} or \code{SummarizedExperiment}
#'   object with \code{computeBanksy} ran.
#' @param use_agf A logical vector specifying whether to use the AGF for
#'   clustering.
#' @param lambda A numeric vector in \eqn{\in [0,1]} specifying a spatial
#'   weighting parameter. Larger values (e.g. 0.8) incorporate more spatial
#'   neighborhood and find spatial domains, while smaller values (e.g. 0.2)
#'   perform spatial cell-typing.
#' @param use_pcs A logical scalar specifying whether to cluster on PCs. If
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
#' @param algo A string scalar specifying the clustering algorithm to use; one
#'   of leiden, louvain, mclust, kmeans.
#' @param k_neighbors An integer vector specifying number of neighbors for
#'   constructing sNN (for louvain / leiden).
#' @param resolution A numeric vector specifying resolution used for clustering
#'   (louvain / leiden).
#' @param leiden.iter An integer scalar specifying the number of leiden
#'   iterations. For running till convergence, set to -1 (leiden).
#' @param kmeans.centers An integer vector specifying the number of kmeans 
#'   clusters (kmeans). 
#' @param mclust.G An integer vector specifying the number of mixture 
#' components (Mclust). 
#' @param M Advanced usage. An integer vector specifying the highest azimuthal
#'   Fourier harmonic to cluster with. If specified, overwrites the
#'   \code{use_agf} argument.
#' @param seed Random seed for clustering. If not specified, no seed is set. 
#' @param ... to pass to methods
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom utils tail
#' @importFrom S4Vectors metadata<-
#'
#' @return A SpatialExperiment / SingleCellExperiment / SummarizedExperiment
#'   object with cluster labels in \code{colData(se)}.
#'
#' @export
#'
#' @examples
#' data(rings)
#' spe <- computeBanksy(rings, assay_name = "counts", M = 1, k_geom = c(15, 30))
#' spe <- runBanksyPCA(spe, M = 1, lambda = c(0, 0.2), npcs = 20)
#' spe <- clusterBanksy(spe, M = 1, lambda = c(0, 0.2), resolution = 1)
#' 
clusterBanksy <-
    function(se,
             use_agf = FALSE,
             lambda = 0.2,
             use_pcs = TRUE, 
             npcs = 20L,
             dimred = NULL, 
             ndims = NULL,
             assay_name = NULL, 
             group = NULL, 
             algo = c("leiden", "louvain", "kmeans", "mclust"),
             k_neighbors = 50,
             resolution = 1,
             leiden.iter = -1,
             kmeans.centers = 5,
             mclust.G = 5, 
             M = NULL,
             seed = NULL,
             ...) {
        
        # Check args
        M <- getM(use_agf, M)
        all_params <- c(as.list(environment(), list(...)))
        all_params$se <- NULL
        algo <- match.arg(algo)
        checkClusterArgs(se = se, all_params = all_params)
        
        # Compute cluster matrices
        params <- expand.grid(lambda = lambda, M = M)
        colnames(params) <- c("lam", "har")
        if (is.null(dimred)) {
            cluster_matrices <- lapply(seq(nrow(params)), function(i) {
                getClusterMatrix(se,
                                 assay_name = assay_name,
                                 M = params[i,]$har,
                                 lambda = params[i,]$lam, 
                                 use_pcs = use_pcs,
                                 npcs = npcs, 
                                 group = group)})
        } else {
            # Dim. reduction provided
            ndims <- checkDimred(se, dimred, ndims)
            cluster_matrices <- list(
                reducedDim(se, dimred)[, seq(ndims), drop=FALSE])
        }
        
        # Run clustering
        labels <- switch(
            algo, 
            leiden = runGraphBased(cluster_matrices, 
                                   algo, 
                                   k_neighbors, 
                                   resolution, 
                                   leiden.iter, 
                                   seed),
            louvain = runGraphBased(cluster_matrices, 
                                    algo, 
                                    k_neighbors, 
                                    resolution, 
                                    leiden.iter, 
                                    seed),
            kmeans = runKmeans(cluster_matrices,
                                  kmeans.centers,
                                  seed, ...),
            mclust = runMclust(cluster_matrices,
                                  mclust.G,
                                  seed, ...)
        )
        
        # Formatting output
        clust_param <- switch(algo,
            leiden = expand.grid(res=resolution, k=k_neighbors, 
                                 lam=lambda, M=M), 
            louvain = expand.grid(res=resolution, k=k_neighbors, 
                                  lam=lambda, M=M), 
            kmeans = expand.grid(kmeans=kmeans.centers, lam=lambda, M=M),
            mclust = expand.grid(mclust=mclust.G, lam=lambda, M=M))
        labels <- do.call(cbind.data.frame, unlistNested(labels))
        cluster_names <- generateClusterNames(rev(clust_param))
        if (!is.null(dimred)) {
            cluster_names <- unique(
                gsub("M\\d+_lam\\d+(\\.\\d+)?", dimred, cluster_names)
            )
        }
        colnames(labels) <- cluster_names
        colData(se) <- cbind(colData(se), labels)
        
        # Log
        metadata(se)$BANKSY_params$algo <- algo
        metadata(se)$BANKSY_params$k_neighbors <- k_neighbors
        metadata(se)$BANKSY_params$resolution <- resolution
        metadata(se)$BANKSY_params$kmeans.centers <- kmeans.centers
        metadata(se)$BANKSY_params$mclust.G <- mclust.G
        metadata(se)$BANKSY_params$cluster_seed <- seed
        
        se
    }


#' @importFrom leidenAlg leiden.community
#' @importFrom igraph as.undirected cluster_louvain
runGraphBased <- function(cluster_matrices, algo, 
                         k_neighbors, resolution, leiden.iter, seed) {
    lapply(cluster_matrices, function(cmat) {
        lapply(k_neighbors, function(k) {
            graph <- as.undirected(getGraph(cmat, k))
            lapply(resolution, function(res){
                verbose.seed(seed)
                memb <- switch(
                    algo,
                    leiden = {
                        cout <- leiden.community(graph, resolution = res,
                                                 n.iterations = leiden.iter)
                        memb <- cout$membership
                        factor(memb, labels = seq(
                            as.numeric(tail(levels(memb), n = 1)) + 1))
                    },
                    louvain = {
                        cout <- cluster_louvain(graph, resolution = res)
                        factor(cout$membership) 
                    })
            })
        })
    })
}

#' @importFrom stats kmeans
runKmeans <- function(cluster_matrices, kmeans.centers, seed, ...) {
    lapply(cluster_matrices, function(cmat) {
        lapply(kmeans.centers, function(k) {
            verbose.seed(seed)
            cout <- kmeans(cmat, centers = k, ...)
            factor(cout$cluster)
        })
    })
}

#' @importFrom mclust Mclust mclustBIC
runMclust <- function(cluster_matrices, mclust.G, seed, ...) {
    lapply(cluster_matrices, function(cmat) {
        lapply(mclust.G, function(G) {
            verbose.seed(seed)
            cout <- Mclust(cmat, G = G)
            factor(cout$classification)
        })
    })
}

# Generate cluster labels
generateClusterNames <- function(df) {
    processRow <- function(row) {
        paste0("clust_", paste0(names(df), row, sep = "", collapse = "_"))
    }
    unlist(lapply(seq_len(nrow(df)), function(i) processRow(df[i, ])))
}

# Unlist a nested list
unlistNested <- function(nested){
    lapply(rapply(nested, enquote, how = "unlist"), eval)
}

# Check arguments to clusterBansky
checkClusterArgs <- function(se, all_params) {
    all_params <- all_params[!vapply(all_params, is.null, logical(1))]
    param_nms <- names(all_params)
    
    stopifnot("use_agf should be a logical vector" = 
                  is.logical(as.logical(all_params$compute_agf)))
    stopifnot("lambda should be a numeric vector with each entry in [0,1]" = 
                  is.numeric(all_params$lambda) & 
                  max(all_params$lambda) <= 1 & 
                  min(all_params$lambda >= 0))
    stopifnot("use_pcs should be an integer scalar" = 
                  is.logical(all_params$use_pcs) & 
                  length(all_params$use_pcs) == 1)
    stopifnot("npcs should be an integer scalar" = 
                  is.integer(as.integer(all_params$npcs)) & 
                  length(all_params$npcs) == 1)

    # Check PCs
    if (all_params$use_pcs & is.null(all_params$dimred)) {
        har <- all_params$M
        lam <- all_params$lambda
        npcs <- all_params$npcs
        checkPCA(se, har, lam, npcs)
    }
}

#' @importFrom SingleCellExperiment reducedDimNames reducedDim
checkPCA <- function(se, har, lam, npcs) {
    pgrid <- expand.grid(lam, har)
    har <- pgrid[, 2]
    lam <- pgrid[, 1]
    param_names <- sprintf("PCA_M%s_lam%s", har, lam)


    pca.names <- sprintf("PCA_M%s_lam%s", har, lam)
    check <- pca.names %in% reducedDimNames(se)

    # Check if PCs in object
    if (any(check == FALSE)) {
        id <- which(check == FALSE)
        err_msg <- paste0(sprintf("(%s,%s)", har[id], lam[id]), collapse = "; ")
        stop("Run PCA for (use_agf,lambda)=", err_msg)
    }

    # Check if sufficient PCs
    pca.mats <- lapply(pca.names, function(nm) reducedDim(se, nm))
    pca.ncols <- vapply(pca.mats, ncol, FUN.VALUE = numeric(1))
    if (any(pca.ncols < npcs)) {
        id <- which(pca.ncols < npcs)
        err_msg <- paste0(sprintf("(%s,%s)", har[id], lam[id]),
            collapse = "; "
        )
        stop(
            "Not enough PCs for (use_agf,lambda)=", err_msg,
            "\nCall runBanksyPCA and increase npcs for these (use_agf,lambdas)."
        )
    }
}

checkDimred <- function(se, dimred, ndims) {
    check <- dimred %in% reducedDimNames(se)
    # Check if reduction in object
    if (check == FALSE) {
        stop("Following dim. reductions not found: ", dimred)
    } else {
        mat <- reducedDim(se, dimred)
        if (is.null(ndims)) ndims <- ncol(mat)
        if (ncol(mat) < ndims) {
            err_msg <- sprintf(
                "Not enough dimensions for %s. Requested %s but found %s.",
                dimred, ndims, ncol(mat)
            )
            stop(err_msg)
        }
    }
    ndims
}


#' @importFrom SingleCellExperiment reducedDim
getClusterMatrix <- function(se, assay_name, M, lambda, use_pcs, npcs, group) {
    if (use_pcs) {
        pca.name <- sprintf("PCA_M%s_lam%s", M, lambda)
        x <- reducedDim(se, pca.name)[, seq(npcs)]
    } else {
        x <-
            t(getBanksyMatrix(
                se,
                assay_name = assay_name,
                M = M,
                lambda = lambda,
                scale = TRUE,
                group = group
            ))
    }
    x
}


#' @importFrom dbscan kNN sNN
#' @importFrom igraph graph_from_data_frame
#' @importFrom data.table setorder data.table
getGraph <- function(x, k) {
    snet <- sNN(kNN(x, k = k), k = k)
    from <- shared <- .N <- NULL
    snet.dt <- data.table(
        from = rep(seq_len(nrow(snet$id)), k),
        to = as.vector(snet$id),
        weight = 1 / (1 + as.vector(snet$dist)),
        distance = as.vector(snet$dist),
        shared = as.vector(snet$shared)
    )
    data.table::setorder(snet.dt, from, -shared)
    snet.dt[, rank := seq_len(.N), by = from]
    snet.dt <- snet.dt[rank <= 3 | shared >= 5]
    graph <- graph_from_data_frame(snet.dt)
    return(graph)
}

#' Relabel cluster labels across parameter runs to maximise their similarity.
#'
#' @param se A \code{SpatialExperiment},
#' \code{SingleCellExperiment} or \code{SummarizedExperiment}
#'   object with cluster labels in \code{colData(se)}.
#' @param map_to A string scalar specify a cluster to map to.
#' @param verbose A logical scalar specifying verbosity.
#'
#' @return A SpatialExperiment / SingleCellExperiment / SummarizedExperiment
#'   object with 'connected' cluster labels in \code{colData(se)}.
#'
#' @export
#'
#' @examples
#' data(rings)
#' spe <- computeBanksy(rings, assay_name = "counts", M = 1, k_geom = c(15, 30))
#' spe <- runBanksyPCA(spe, M = 1, lambda = c(0, 0.2), npcs = 20)
#' spe <- clusterBanksy(spe, M = 1, lambda = c(0, 0.2), resolution = 1)
#' spe <- connectClusters(spe)
#'
connectClusters <- function(se, map_to = NULL, verbose = TRUE) {
    clust_nm <- colnames(colData(se))
    clust_nm <- clust_nm[grep("^clust", clust_nm)]
    clust_df <- data.frame(colData(se)[, clust_nm])

    if (is.null(map_to)) {
        # Perform a cascading map based on number of levels
        n_clust <- apply(clust_df, 2, function(x) length(unique(x)))
        map_order <- order(n_clust)

        for (i in seq(length(map_order) - 1)) {
            curr_many <- map_order[i + 1]
            curr_few <- map_order[i]
            if (verbose) {
                message(clust_nm[curr_many], " --> ", clust_nm[curr_few])
            }
            clust_df[, curr_many] <-
                mapToSeed(clust_df[, curr_few], clust_df[, curr_many])
        }
    } else {
        found <- match(map_to, clust_nm)
        if (is.na(found)) stop(map_to, " not found in cluster names /^clust/.")
        seed <- clust_df[, found]
        for (i in setdiff(seq(length(clust_nm)), found)) {
            clust_df[, i] <- mapToSeed(seed, clust_df[, i])
        }
    }

    colData(se)[, clust_nm] <- clust_df
    se
}

#' @importFrom RcppHungarian HungarianSolver
mapToSeed <- function(seed, val) {
    max_seed <- max(as.numeric(seed))
    max_val <- max(as.numeric(val))
    max_both <- max(max_seed, max_val)

    # Get min cost matching
    con.mat <- table(seed, val)
    cost.mat <- max(con.mat) - con.mat
    matching <- HungarianSolver(cost.mat)$pairs

    # Remove unmapped labels if val has fewer labels than seed
    if (any(matching[, 2, drop = FALSE] == 0)) {
        matching <- matching[matching[, 2, drop = FALSE] != 0, ]
    }
    from <- matching[, 2, drop = FALSE]
    to <- matching[, 1, drop = FALSE]

    # Deal with the extra unmapped labels in val
    unmapped_from <- setdiff(seq(max_val), from)
    unmapped_to <- seq(length(from) + 1, length.out = length(unmapped_from))
    from <- c(from, unmapped_from)
    to <- c(to, unmapped_to)

    val <- factor(val, labels = to[order(from)])
    val <- factor(val, levels = seq_len(max_both))

    table(seed, val)
    val
}


#' k-Nearest neighbor cluster label smoothing.
#'
#' @details 
#' As described in SpiceMix (https://doi.org/10.1038/s41588-022-01256-z). 
#' Implemented for labels that can be coerced to numeric only.
#'
#' @param se A \code{SpatialExperiment},
#' \code{SingleCellExperiment} or \code{SummarizedExperiment}
#'   object with cluster labels in \code{colData(se)}.
#' @param cluster_names A string vector of label names to smooth. If NULL,
#'   smooths labels in colData(se) matching /^clust/
#' @param coord_names A string vector specifying the names in \code{colData}
#'   corresponding to spatial coordinates.
#' @param k An integer scalar specifying number of neighbors for smooething.
#' @param prop_thres A numeric scalar \eqn{\in [0,1]} specifying a label
#'   proportions threshold If the fraction of neighbors with a certain label
#'   exceeds this proportion, change the label of the current sample
#'   (default: 0.5).
#' @param max_iter An integer scalar specifying the max number of smoothing
#'   iterations. Set to -1 for smoothing to convergence.
#' @param verbose A logical scalar specifying verbosity.

#' @return A SpatialExperiment / SingleCellExperiment / SummarizedExperiment
#'   object with smoothed cluster labels in \code{colData(se)} suffixed with
#'   '_smooth'.
#'
#' @export
#'
#' @examples
#' data(rings)
#' spe <- computeBanksy(rings, assay_name = "counts", M = 1, k_geom = c(15, 30))
#' spe <- runBanksyPCA(spe, M = 1, lambda = 0.2, npcs = 20)
#' spe <- clusterBanksy(spe, M = 1, lambda = 0.2, resolution = 1)
#' spe <- smoothLabels(spe, cluster_names = "clust_M1_lam0.2_k50_res1")
#'
smoothLabels <- function(se,
                         cluster_names = NULL,
                         coord_names = NULL,
                         k = 15L,
                         prop_thres = 0.5,
                         max_iter = 10,
                         verbose = TRUE) {
    # Parse args
    clust_nm <- colnames(colData(se))
    if (is.null(cluster_names)) {
        clust_nm <- clust_nm[grep("^clust", clust_nm)]
    } else {
        found <- cluster_names %in% clust_nm
        if (any(!found)) {
            err_msg <- cluster_names[!found]
            stop(
                "The following cluster names were not found in colData: ",
                err_msg
            )
        }
        clust_nm <- cluster_names[cluster_names %in% clust_nm]
    }
    if (length(clust_nm) == 0) stop("No clustering runs found.")

    # Extract clusters
    clust_df <- colData(se)[, clust_nm, drop = FALSE]

    # Smooth for each cluster
    for (i in seq(length(clust_nm))) {
        new_name <- paste0(clust_nm[i], "_smooth")
        colData(se)[, new_name] <- smoother(
            clust_df[, clust_nm[i]],
            getLocs(se, coord_names),
            k = k, prop_thres = prop_thres,
            max_iter = max_iter, verbose = verbose
        )
    }
    se
}


#' @importFrom dbscan kNN
smoother <-
    function(labels_curr,
             locs,
             k = 10L,
             prop_thres = 0.5,
             max_iter = 10,
             verbose = TRUE) {
        # Get neighbors
        knn <- kNN(locs, k = k)$id
        N <- nrow(knn)

        labels_raw <- labels_curr
        labels_update <- labels_curr

        if (max_iter == -1) max_iter <- Inf

        iter <- 1
        while (iter < max_iter) {
            if (verbose) {
                message("Iteration ", iter)
            }

            # Iterate across cells
            for (i in seq(N)) {
                # Get neighbors
                neighbor_labels <- labels_curr[c(i, knn[i, ])]
                neighbor_props <- table(neighbor_labels) /
                    length(neighbor_labels)
                # Change label based on condition
                if (any(neighbor_props > prop_thres)) {
                    labels_update[i] <- as.numeric(
                        names(which.max(neighbor_props))
                    )
                } else {
                    labels_update[i] <- labels_curr[i]
                }
            }

            change <- sum(labels_update != labels_curr)
            if (verbose) {
                message("Change: ", change)
            }

            if (change == 0) {
                break
            }

            labels_curr <- labels_update
            iter <- iter + 1
        }
        labels_update
    }


#' Compare cluster outputs based on various clustering comparison measures.
#'
#' @param se A \code{SpatialExperiment},
#' \code{SingleCellExperiment} or \code{SummarizedExperiment}
#'   object with cluster labels in \code{colData(se)}.
#' @param func A string scalar specifying what clustering comparison measure to
#'   compute. See \code{?aricode} for more information.
#' @param digits An integer scalar specifying the number of digits to round to.
#'
#' @importFrom aricode AMI ARI MARI MARIraw RI NID NMI NVI
#' @importFrom utils combn
#'
#' @return A matrix of cluster comparison measures.
#'
#' @export
#'
#' @examples
#' data(rings)
#' spe <- computeBanksy(rings, assay_name = "counts", M = 1, k_geom = c(15, 30))
#' spe <- runBanksyPCA(spe, M = 1, lambda = 0.2, npcs = 20)
#' spe <- clusterBanksy(spe, M = 1, lambda = 0.2, resolution = c(0.1, 1))
#' spe <- connectClusters(spe)
#' compareClusters(spe)
#'
compareClusters <-
    function(se,
             func = c(
                 "ARI", "AMI", "MARI", "MARIraw", "RI",
                 "NID", "NMI", "NVI"
             ),
             digits = 3) {
        if (length(func) > 1) func <- func[1]
        func <- match.arg(arg = func)
        measure_func <- get(func, mode = "function")

        clust_nm <- colnames(colData(se))
        clust_nm <- clust_nm[grep("^clust", clust_nm)]
        clust_df <- data.frame(colData(se)[, clust_nm, drop = FALSE])
        n_clust <- ncol(clust_df)

        if (n_clust < 2) {
            stop("ARI will only be calculated for at least 2 clustering runs.")
        }

        comb <- combn(names(clust_df), 2)
        n_comb <- ncol(comb)
        res <- numeric(length = ncol(comb))
        for (i in seq_len(n_comb)) {
            res[i] <- measure_func(
                clust_df[[comb[1, i]]],
                clust_df[[comb[2, i]]]
            )
        }

        res_mat <- diag(nrow = n_clust, ncol = n_clust)
        res_mat[lower.tri(res_mat)] <- res
        rownames(res_mat) <- colnames(res_mat) <- colnames(clust_df)
        res_mat <- res_mat + t(res_mat) - diag(nrow = n_clust, ncol = n_clust)
        res_mat <- round(res_mat, digits)
        res_mat
    }

#' Get names of clustering runs.
#'
#' @param se A \code{SpatialExperiment},
#' \code{SingleCellExperiment} or \code{SummarizedExperiment}
#'   object with \code{clusterBanksy} ran.
#'
#' @importFrom SummarizedExperiment colData
#'
#' @return A character vector of names of clustering runs.
#'
#' @export
#'
#' @examples
#' data(rings)
#' spe <- computeBanksy(rings, assay_name = "counts", M = 1, k_geom = c(15, 30))
#' spe <- runBanksyPCA(spe, M = 1, lambda = c(0, 0.2), npcs = 20)
#' spe <- clusterBanksy(spe, M = 1, lambda = c(0, 0.2), resolution = 1)
#' clusterNames(spe)
#'
clusterNames <- function(se) {
    grep("^clust_", colnames(colData(se)), value = TRUE)
}
