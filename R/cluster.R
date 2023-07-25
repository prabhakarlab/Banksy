#' Cluster based on joint expression matrix
#'
#' @param se A \code{SpatialExperiment},
#' \code{SingleCellExperiment} or \code{SummarizedExperiment}
#'   object with \code{computeBanksy} ran.
#' @param use_agf A logical vector specifying whether to use the AGF for 
#'   clustering.
#' @param lambda A numeric scalar in \eqn{\in [0,1]} specifying a spatial
#'   weighting parameter. Larger values incorporate more spatial neighborhood
#'   information.
#' @param use_pcs A logical scalar specifying whether to cluster on PCs. If
#'   FALSE, runs on the BANKSY matrix.
#' @param npcs An integer scalar specifying the number of principal components
#'   to use if \code{use_pcs} is TRUE.
#' @param assay_name A string scalar specifying the name of the assay used in
#'   \code{computeBanksy}.
#' @param group A string scalar specifying a grouping variable for samples in
#'   \code{se}. This is used to scale the samples in each group separately.
#' @param algo A string scalar specifying the clustering algorithm to use; one
#'   of leiden, louvain, mclust, kmeans.
#' @param k_neighbors An integer scalar specifying number of neighbors for
#'   constructing sNN (for louvain / leiden).
#' @param resolution A numeric scalar specifying resolution used for clustering
#'   (leiden).
#' @param leiden.iter An integer scalar specifying the number of leiden
#'   iterations. For running till convergence, set to -1 (leiden).
#' @param M Advanced usage. An integer vector specifying the highest azimuthal
#'   Fourier harmonic to cluster with. If specified, overwrites the 
#'   \code{use_agf} argument.
#' @param seed Random seed.
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
clusterBanksy <- function(se,
                          use_agf = TRUE,
                          lambda = 0.2,
                          algo = "leiden",
                          use_pcs = TRUE,
                          npcs = 20,
                          assay_name = NULL,
                          group = NULL,
                          k_neighbors = 50,
                          resolution = 1,
                          leiden.iter = -1,
                          M = NULL,
                          seed = NULL,
                          ...) {
    M <- getM(use_agf, M)
    all_params <- c(as.list(environment(), list(...)))
    all_params$se <- NULL
    checkClusterArgs(se = se, all_params = all_params)

    if (algo == "leiden") {
        cout <- runLeiden(
            se,
            assay_name = assay_name,
            M = M,
            lambda = lambda,
            use_pcs = use_pcs,
            npcs = npcs,
            group = group,
            k_neighbors = k_neighbors,
            resolution = resolution,
            leiden.iter = leiden.iter,
            seed = seed
        )
    }

    colData(se) <- cbind(colData(se), cout)

    # Log
    metadata(se)$BANKSY_params$algo <- algo
    metadata(se)$BANKSY_params$k_neighbors <- k_neighbors
    metadata(se)$BANKSY_params$resolution <- resolution
    metadata(se)$BANKSY_params$cluster_seed <- seed

    se
}

#' @importFrom leidenAlg leiden.community
#' @importFrom igraph as.undirected
runLeiden <- function(se,
                      assay_name,
                      M,
                      lambda,
                      use_pcs,
                      npcs,
                      group = group,
                      k_neighbors,
                      resolution,
                      leiden.iter,
                      seed) {
    max.iters <- prod(
        length(M),
        length(lambda),
        length(k_neighbors),
        length(resolution)
    )
    params <- expand.grid(resolution, k_neighbors, lambda, M)
    colnames(params) <- c("res", "k", "lam", "har")
    param_names <- sprintf(
        "clust_M%s_lam%s_k%s_res%s",
        params$har, params$lam, params$k, params$res
    )

    out <- lapply(M, function(har) {
        lapply(lambda, function(lam) {
            x <- getClusterMatrix(se,
                assay_name = assay_name,
                M = har,
                lambda = lam,
                use_pcs = use_pcs,
                npcs = npcs,
                group = group
            )
            lapply(k_neighbors, function(k) {
                graph <- as.undirected(getGraph(x, k))
                lapply(resolution, function(res) {
                    if (!is.null(seed)) {
                        set.seed(seed)
                        message("Using seed=", seed)
                    }
                    cout <- leiden.community(
                        graph,
                        resolution = res, n.iterations = leiden.iter
                    )
                    memb <- cout$membership
                    factor(memb, labels = seq(
                        as.numeric(tail(levels(memb), n = 1)) + 1
                    ))
                })
            })
        })
    })
    out <- do.call(cbind.data.frame, out)
    colnames(out) <- param_names
    out
}


checkClusterArgs <- function(se, all_params) {
    
    all_params <- all_params[!vapply(all_params, is.null, logical(1))]
    param_nms <- names(all_params)

    # Check params for clustering algos.
    if (all_params$algo == "leiden") {
        if (!("k_neighbors" %in% param_nms)) {
            stop("Specify k_neighbors")
        }
        if (!("resolution" %in% param_nms)) {
            stop("Specify resolution")
        }
    } else {
        stop('Invalid clustering algo.')
    }

    # Check PCs
    if (all_params$use_pcs) {
        
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
        stop("Run PCA for (M,lambda)=", err_msg)
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
            "Not enough PCs for (M,lambda)=", err_msg,
            "\nCall runBanksyPCA and increase npcs for these (M,lambdas)."
        )
    }
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
    if (any(matching[,2,drop=FALSE] == 0)) {
        matching <- matching[matching[,2,drop=FALSE] != 0, ]
    }
    from <- matching[,2,drop=FALSE]
    to <- matching[,1,drop=FALSE]

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


#' Label smoothing as described in SpiceMix 
#' (https://doi.org/10.1038/s41588-022-01256-z).
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
                if (any(neighbor_props > 0.5)) {
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


#' Compare cluster outputs.
#'
#' @param se A \code{SpatialExperiment},
#' \code{SingleCellExperiment} or \code{SummarizedExperiment}
#'   object with cluster labels in \code{colData(se)}.
#' @param func A string scalar specifying what clustering comparison measure to
#'   compute. See `?aricode` for more information.
#' @param digits An integer scalar specifying the number of digits to round to.
#'
#' @importFrom aricode AMI ARI MARI MARIraw RI NID NMI NVI
#' @importFrom utils combn
#'
#' @return A matirx of cluster comparison measures.
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
        func <- match.arg(
            arg = func,
            choices = c(
                "ARI", "AMI", "MARI", "MARIraw", "RI",
                "NID", "NMI", "NVI"
            )
        )
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
