# ClusterBanksy
# ConnectClusters

#' Cluster based on joint expression matrix
#'
#' @param bank BanksyObject
#' @param lambda (numeric) weighting parameter
#' @param M (numeric) compute up to the k-th azimuthal fourier harmonic (default: 1)
#' @param pca (logical) runs clustering on PCs (TRUE) or BANKSY matrix (FALSE)
#' @param npcs (numeric) number of pcs to use for clustering
#' @param method (character) one of leiden, louvain, mclust, kmeans
#' @param k.neighbors (numeric) parameter for constructing sNN (for louvain / leiden)
#' @param resolution (numeric) parameter used for clustering (leiden)
#' @param leiden.iter (numeric) number of leiden iterations (leiden)
#' @param mclust.G (numeric) number of mixture components (mclust)
#' @param kmeans.centers (numeric) number of clusters (kmeans)
#' @param kmeans.iter.max (numeric) max number of iterations (kmeans)
#' @param seed seed
#' @param ... to pass to methods
#'
#' @return BanksyObject with cluster labels in meta.data
#'
#' @export
#'
#' @examples
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' bank <- NormalizeBanksy(bank)
#' bank <- ComputeBanksy(bank)
#' bank <- ScaleBanksy(bank)
#' bank <- RunBanksyPCA(bank, lambda = 0.3)
#' bank <- ClusterBanksy(bank, lambda = 0.3, npcs = 20, k.neighbors = 50, resolution = 0.5)
#'
ClusterBanksy <-
    function(bank,
             lambda = 0.2,
             M = 1,
             pca = TRUE,
             npcs = 20,
             method = c('leiden', 'louvain', 'mclust', 'kmeans'),
             k.neighbors = 50,
             resolution = 1,
             leiden.iter = -1,
             mclust.G = NULL,
             kmeans.centers = NULL,
             kmeans.iter.max = 10,
             seed = 42,
             verbose = TRUE,
             ...) {
        method <- checkMethod(method)
        params <- c(as.list(environment(), list(...)))
        checkArgs(bank, method, lambda, M, pca, npcs, params)
        
        if (method == 'kmeans') {
            bank <- runKmeans(bank, lambda, M, pca, npcs, 
                              kmeans.centers, kmeans.iter.max, ...)
        }
        
        if (method == 'mclust') {
            bank <- runMclust(bank, lambda, M, pca, npcs, mclust.G, ...)
        }
        
        if (method == 'leiden') {
            bank <- runLeiden(bank, lambda, M, pca, npcs, 
                              k.neighbors, resolution, leiden.iter, verbose, ...)
        }
        
        if (method == 'louvain') {
            bank <- runLouvain(bank, lambda, M, pca, npcs, k.neighbors, resolution)
        }
        
        return(bank)
    }

checkMethod <- function(method) {
    if (length(method) > 1) {
        return('leiden')
    }
    
    methods <- c('kmeans', 'mclust', 'louvain', 'leiden')
    valid <- pmatch(method, methods)
    if (is.na(valid))
        stop('Specify a valid method. One of ',
             paste(methods, collapse = ', '))
    
    return(method)
}

checkArgs <- function(bank, method, lambda, M, pca, npcs, call) {
    call <- call[!sapply(call, is.null)]
    args <- names(call)
    
    if (pca) {
        # Check if pca has been run for lambda
        params <- expand.grid(M, lambda)
        pca.names <- paste0('pca_M', params[, 1], '_lam', params[, 2])
        check <- pca.names %in% names(bank@reduction)
        # Check if enough pcs
        if (any(check == FALSE)) {
            id <- which(check == FALSE)
            stop(
                'Run PCA for M=',
                paste(params[id, ][, 1], collapse = ','),
                ' lambda=',
                paste(params[id, ][, 2], collapse = ',')
            )
        }
        pca.mats <-
            lapply(pca.names, function(nm)
                bank@reduction[[nm]]$x)
        pca.ncols <- vapply(pca.mats, ncol, FUN.VALUE = numeric(1))
        if (any(pca.ncols < npcs)) {
            id <- which(pca.ncols < npcs)
            stop(
                'Not enough PCs for M=',
                paste(params[id, ][, 1], collapse = ','),
                ' lambda=',
                paste(params[id, ][, 2], collapse = ','),
                '\nCall RunBanksyPCA and increase npcs for these (M,lambdas).'
            )
        }
    }
    
    if (method == 'kmeans') {
        if (!('kmeans.centers' %in% args))
            stop('Specify kmeans.centers')
    }
    
    if (method == 'mclust') {
        if (!('mclust.G' %in% args))
            stop('Specify mclust.G')
    }
    
    if (method == 'leiden') {
        if (!('k.neighbors' %in% args))
            stop('Specify k.neighbors')
        if (!('resolution' %in% args))
            stop('Specify resolution')
    }
    
    if (method == 'louvain') {
        if (!('k.neighbors' %in% args))
            stop('Specify k.neighbors')
        if (!('resolution' %in% args))
            stop('Specify resolution')
    }
    
}

getClusterMatrix <- function(bank, lambda, M, pca, npcs) {
    if (pca) {
        pca.name <- paste0('pca_M', M, '_lam', lambda)
        x <- bank@reduction[[pca.name]]$x[, seq_len(npcs)]
    } else {
        x <- t(getBanksyMatrix(bank, lambda = lambda, M = M)$expr)
    }
    return(x)
}

#' @importFrom dbscan kNN sNN
#' @importFrom igraph graph_from_data_frame
#' @importFrom data.table setorder data.table
getGraph <- function(x, k) {
    snet <- sNN(kNN(x, k = k), k = k)
    from <- shared <- .N <- NULL
    snet.dt = data.table(
        from = rep(seq_len(nrow(snet$id)), k),
        to = as.vector(snet$id),
        weight = 1 / (1 + as.vector(snet$dist)),
        distance = as.vector(snet$dist),
        shared = as.vector(snet$shared)
    )
    data.table::setorder(snet.dt, from,-shared)
    snet.dt[, rank := seq_len(.N), by = from]
    snet.dt <- snet.dt[rank <= 3 | shared >= 5]
    graph <- graph_from_data_frame(snet.dt)
    return(graph)
    
}

#' @importFrom stats kmeans
runKmeans <- function(bank,
                      lambda,
                      M,
                      pca,
                      npcs,
                      kmeans.centers,
                      kmeans.iter.max,
                      ...) {
    max.iters <- prod(length(lambda), length(kmeans.centers), length(M))
    iter <- 1
    message('Iteration ', iter, ' out of ', max.iters)
    for (har in M) {
        for (lam in lambda) {
            x <- getClusterMatrix(bank, lam, har, pca, npcs)
            for (k in kmeans.centers) {
                out <- kmeans(x,
                              centers = k,
                              iter.max = kmeans.iter.max,
                              ...)
                clust.name <- paste0('clust_M', har, '_lam', lam, '_kmeans', k)
                bank@meta.data[[clust.name]] <- out$cluster
                iter <- iter + 1
                if (iter <= max.iters)
                    message('Iteration ', iter, ' out of ', max.iters)
            }
        }
    }
    return(bank)
}

#' @importFrom mclust Mclust mclustBIC
runMclust <- function(bank, lambda, M, pca, npcs, mclust.G, ...) {
    max.iters <- prod(length(lambda), length(mclust.G), length(M))
    iter <- 1
    message('Iteration ', iter, ' out of ', max.iters)
    for (har in M) {
        for (lam in lambda) {
            x <- getClusterMatrix(bank, lam, har, pca, npcs)
            for (G in mclust.G) {
                out <- Mclust(x, G = G, ...)
                clust.name <- paste0('clust_M', har, '_lam', lam, '_mclust', G)
                bank@meta.data[[clust.name]] <- out$classification
                iter <- iter + 1
                if (iter <= max.iters)
                    message('Iteration ', iter, ' out of ', max.iters)
            }
        }
    }
    return(bank)
}

#' @importFrom leidenAlg leiden.community
#' @importFrom progress progress_bar
runLeiden <- function(bank,
                      lambda,
                      M,
                      pca,
                      npcs,
                      k.neighbors,
                      resolution,
                      leiden.iter,
                      verbose) {
    max.iters <-
        prod(length(lambda),
             length(k.neighbors),
             length(resolution),
             length(M))
    pb <- progress_bar$new(
        format = " [:bar] :percent eta: :eta",
        total = max.iters, clear = FALSE, width = 60)
    for (har in M) {
        for (lam in lambda) {
            x <- getClusterMatrix(bank, lam, har, pca, npcs)
            for (k in k.neighbors) {
                graph <- getGraph(x, k)
                for (res in resolution) {
                    if (verbose) pb$tick()
                    out <-
                        leiden.community(graph,
                                         resolution = res,
                                         n.iterations = leiden.iter)
                    clust.name <- paste0('clust_M', har, '_lam', lam, '_k', k, '_res', res)
                    bank@meta.data[[clust.name]] <-
                        as.numeric(out$membership)
                    # iter <- iter + 1
                }
            }
        }
    }
    return(bank)
}

#' @importFrom igraph cluster_louvain membership as.undirected
runLouvain <-
    function(bank,
             lambda,
             M,
             pca,
             npcs,
             k.neighbors,
             resolution) {
        max.iters <-
            prod(length(lambda),
                 length(k.neighbors),
                 length(resolution),
                 length(M))
        iter <- 1
        message('Iteration ', iter, ' out of ', max.iters)
        for (har in M) {
            for (lam in lambda) {
                x <- getClusterMatrix(bank, lam, har, pca, npcs)
                for (k in k.neighbors) {
                    graph <- as.undirected(getGraph(x, k))
                    for (res in resolution) {
                        out <- membership(cluster_louvain(graph, resolution = res))
                        clust.name <- paste0('clust_M', har, '_lam', lam, '_k', k, '_res', res, '_louvain')
                        bank@meta.data[[clust.name]] <- as.numeric(out)
                        iter <- iter + 1
                        if (iter <= max.iters)
                            message('Iteration ', iter, ' out of ', max.iters)
                    }
                }
            }
        }
        return(bank)
    }


#' Harmonize cluster labels among parameter runs
#'
#' @param bank Banksy Object
#' @param map.to (character) specify a cluster to map to
#'
#' @return BanksyObject with harmonized cluster labels
#'
#' @export
#'
#' @examples
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' bank <- NormalizeBanksy(bank)
#' bank <- ScaleBanksy(bank)
#' bank <- ComputeBanksy(bank)
#' bank <- RunBanksyPCA(bank, lambda = 0.2)
#' bank <- ClusterBanksy(bank, lambda = 0.2, npcs = 20, k.neighbors = 50, resolution = c(0.5,1.5))
#' bank <- ConnectClusters(bank)
#'
ConnectClusters <- function(bank, map.to = NULL) {
    clusters <- bank@meta.data[, clust.names(bank)]
    order.name <- names(clusters)
    seed.name <- getSeed(clusters, map.to)
    seed <- clusters[[seed.name]]
    values <-
        clusters[,-match(seed.name, names(clusters)), drop = FALSE]
    new.clusters <-
        data.frame(apply(values, 2, function(x)
            mapToSeed(x, seed)))
    new.clusters[[seed.name]] <- seed
    new.clusters <- new.clusters[, order.name]
    bank@meta.data[, clust.names(bank)] <- new.clusters
    
    return(bank)
}

#' @importFrom stats median
getSeed <- function(clusters, map.to) {
    if (!is.null(map.to)) {
        found <- map.to %in% names(clusters)
        if (!found)
            stop(map.to, ' not found in cluster names.')
        seed <- map.to
    } else {
        n.clust <- apply(clusters, 2, function(x)
            length(unique(x)))
        med <- median(n.clust)
        seed <- names(which.min(abs(n.clust - med)))
    }
    return(seed)
}

#' @importFrom plyr mapvalues
#' @importFrom RcppHungarian HungarianSolver
mapToSeed <- function(val, seed) {
 
    # seed and val must be numeric for this to work (because max(seed) is used 
    # for assigning unmatched values)
    con.mat <- table(val, seed)
    cost.mat <- max(con.mat) - con.mat
    matching <- HungarianSolver(cost.mat)$pairs
    
    # split into matched and unmatched
    matched <- matching[!(matching[, 2] == 0), , drop = FALSE]
    unmatched <- matching[matching[, 2] == 0, , drop = FALSE]
    ## Mapping fr more clusters to fewer
    unmapped <- which(matching[, 2] == 0)
    impute <- max(seed) + seq_len(length(unmapped))
    unmatched[, 2] <- impute
    
    # matched is currently the indices of the rows and cols of the contingency matrix
    matched.names.from <-
        c(as.numeric(rownames(cost.mat)[matched[, 1]]), unmatched[, 1])
    matched.names.to <-
        c(as.numeric(colnames(cost.mat)[matched[, 2]]), unmatched[, 2])
    
    new.val <-
        mapvalues(x = val, from = matched.names.from, to = matched.names.to)
    return(new.val)
}


#' Calculate adjusted rand index for harmonized clusters
#'
#' @param bank BanksyObject with harmonized clusters
#' @param digits (numeric) number of digits to round ARI to
#'
#' @return matrix of ARI
#'
#' @importFrom mclust adjustedRandIndex
#' @importFrom utils combn
#'
#' @export
#'
#' @examples
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' bank <- NormalizeBanksy(bank)
#' bank <- ScaleBanksy(bank)
#' bank <- ComputeBanksy(bank)
#' bank <- RunBanksyPCA(bank, lambda = 0.2)
#' bank <- ClusterBanksy(bank, lambda = 0.2, npcs = 20, k.neighbors = 50, resolution = c(0.5,1.5))
#' ari <- getARI(bank)
#' ari
#'
getARI <- function(bank, digits = 3) {
    clust <- bank@meta.data[, clust.names(bank), drop = FALSE]
    n.clust <- ncol(clust)
    if (n.clust < 2) {
        stop('ARI will only be calculated for at least 2 clustering runs.')
    }
    comb <- combn(names(clust), 2)
    n.comb <- ncol(comb)
    ari <- numeric(length = ncol(comb))
    for (i in seq_len(n.comb)) {
        ari[i] <- adjustedRandIndex(clust[[comb[1, i]]], clust[[comb[2, i]]])
    }
    
    ari.mat <- diag(nrow = n.clust, ncol = n.clust)
    ari.mat[lower.tri(ari.mat)] <- ari
    rownames(ari.mat) <- colnames(ari.mat) <- colnames(clust)
    ari.mat <-
        ari.mat + t(ari.mat) - diag(nrow = n.clust, ncol = n.clust)
    ari.mat <- round(ari.mat, digits)
    return(ari.mat)
}
