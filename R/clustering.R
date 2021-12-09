# ClusterBanksy
# ConnectClusters

#' Cluster based on joint expression matrix
#'
#' @param bank BanksyObject
#' @param lambda weighting parameter
#' @param pca if TRUE, runs clustering on PCA, else runs on Banksy matrix
#' @param npcs number of pcs to use for clustering
#' @param method one of leiden, mclust, kmeans
#' @param k.neighbors leiden - parameter for constructing shared nearest neighbor network
#' @param resolution leiden - parameter used for clustering
#' @param leiden.iter leiden - number of leiden iterations
#' @param mclust.G mclust - number of mixture components G
#' @param kmeans.centers kmeans - number of clusters
#' @param kmeans.iter.max kmeans - max number of iterations
#' @param ... to pass to methods
#'
#' @return BanksyObject with cluster labels in meta.data
#'
#' @export
ClusterBanksy <- function(bank, lambda = 0.25, pca = TRUE, npcs = 30,
                          method = c('leiden', 'mclust', 'kmeans'),
                          k.neighbors = NULL, resolution = NULL,
                          leiden.iter = -1, mclust.G = NULL,
                          kmeans.centers = NULL, kmeans.iter.max = 10, ...) {

  method <- checkMethod(method)
  checkArgs(bank, method, lambda, pca, match.call())

  if (method == 'kmeans') {
    bank <- runKmeans(bank, lambda, pca, npcs,
                      kmeans.centers, kmeans.iter.max, ...)
  }

  if (method == 'mclust') {
    bank <- runMclust(bank, lambda, pca, npcs, mclust.G, ...)
  }

  if (method == 'leiden') {
    bank <- runLeiden(bank, lambda, pca, npcs,
                      k.neighbors, resolution, leiden.iter, ...)
  }

  return(bank)
}

checkMethod <- function(method) {

  if (length(method) > 1) {
    return('leiden')
  }

  methods <- c('kmeans', 'mclust', 'louvain', 'leiden')
  valid <- pmatch(method, methods)
  if (is.na(valid)) stop('Specify a valid method. One of ',
                         paste(methods, collapse = ', '))

  return(method)
}

checkArgs <- function(bank, method, lambda, pca, call) {

  args <- names(call)

  if (pca) {
    pca.names <- paste0('pca_', lambda)
    check <- pca.names %in% names(bank@reduction)
    if (any(check == FALSE)) {
      id <- which(check == FALSE)
      stop('Run PCA for lambda=', paste(lambda[id], collaspe = ''))
    }
  }

  if (method == 'kmeans') {
    if (!('kmeans.centers' %in% args)) stop('Specify kmeans.centers')
  }

  if (method == 'mclust') {
    if (!('mclust.G' %in% args)) stop('Specify mclust.G')
  }

  if (method == 'leiden') {
    if (!('k.neighbors' %in% args)) stop('Specify k.neighbors')
    if (!('resolution' %in% args)) stop('Specify resolution')
  }

}

getClusterMatrix <- function(bank, lambda, pca, npcs) {

  if (pca) {
    pca.name <- paste0('pca_', lambda)
    x <- bank@reduction[[pca.name]]$x[,seq_len(npcs)]
  } else {
    x <- getBanksyMatrix(bank, lambda = lambda)$expr
  }
  return(x)
}

#' @importFrom dbscan kNN sNN
#' @importFrom igraph graph_from_data_frame
#' @importFrom data.table setorder data.table
getGraph <- function(x, k) {

  snet <- sNN(kNN(x, k = k), k = k)
  from <- shared <- .N <- NULL
  snet.dt = data.table(from = rep(1:nrow(snet$id), k),
                       to = as.vector(snet$id),
                       weight = 1/(1 + as.vector(snet$dist)),
                       distance = as.vector(snet$dist),
                       shared = as.vector(snet$shared))
  data.table::setorder(snet.dt, from, -shared)
  snet.dt[, rank := 1:.N, by = from]
  snet.dt <- snet.dt[rank <= 3 | shared >= 5]
  graph <- graph_from_data_frame(snet.dt)
  return(graph)

}

#' @importFrom stats kmeans
runKmeans <- function(bank, lambda, pca, npcs,
                      kmeans.centers, kmeans.iter.max, ...) {

  max.iters <- prod(length(lambda), length(kmeans.centers))
  iter <- 1
  message('Iteration ', iter, ' out of ', max.iters)
  for (lam in lambda) {
    x <- getClusterMatrix(bank, lam, pca, npcs)
    for (k in kmeans.centers) {
      out <- kmeans(x, centers = k, iter.max = kmeans.iter.max, ...)
      clust.name <- paste0('clust_lam', lam, '_kmeans', k)
      bank@meta.data[[clust.name]] <- out$cluster
      iter <- iter + 1
      if (iter <= max.iters) message('Iteration ', iter, ' out of ', max.iters)
    }
  }
  return(bank)
}

#' @importFrom mclust Mclust mclustBIC
runMclust <- function(bank, lambda, pca, npcs, mclust.G, ...) {

  max.iters <- prod(length(lambda), length(mclust.G))
  iter <- 1
  message('Iteration ', iter, ' out of ', max.iters)
  for (lam in lambda) {
    x <- getClusterMatrix(bank, lam, pca, npcs)
    for (G in mclust.G) {
      out <- Mclust(x, G = G, ...)
      clust.name <- paste0('clust_lam', lam, '_mclust', G)
      bank@meta.data[[clust.name]] <- out$classification
      iter <- iter + 1
      if (iter <= max.iters) message('Iteration ', iter, ' out of ', max.iters)
    }
  }
  return(bank)
}

#' @importFrom leidenAlg leiden.community
runLeiden <- function(bank, lambda, pca, npcs,
                      k.neighbors, resolution, leiden.iter) {

  max.iters <- prod(length(lambda), length(k.neighbors), length(resolution))
  iter <- 1
  message('Iteration ', iter, ' out of ', max.iters)
  for (lam in lambda) {
    x <- getClusterMatrix(bank, lam, pca, npcs)
    for (k in k.neighbors) {
      graph <- getGraph(x, k)
      for (res in resolution) {
        out <- leiden.community(graph, resolution = res, n.iterations = leiden.iter)
        clust.name <- paste0('clust_lam', lam, '_k', k, '_res', res)
        bank@meta.data[[clust.name]] <- as.numeric(out$membership)
        iter <- iter + 1
        if (iter <= max.iters) message('Iteration ', iter, ' out of ', max.iters)
      }
    }
  }
  return(bank)
}

#' Harmonize cluster labels among parameter runs
#'
#' @param bank Banksy Object
#' @param map.to specify a cluster to map to
#'
#' @return BanksyObject with harmonized cluster labels
#'
#' @export
ConnectClusters <- function(bank, map.to = NULL) {

  clusters <- bank@meta.data[, grepl('^clust', names(bank@meta.data))]
  order.name <- names(clusters)
  seed.name <- getSeed(clusters, map.to)
  seed <- clusters[[seed.name]]
  values <- clusters[, -match(seed.name, names(clusters)), drop = FALSE]
  new.clusters <- data.frame(apply(values, 2, function(x) mapToSeed(x, seed)))
  new.clusters[[seed.name]] <- seed
  new.clusters <- new.clusters[, order.name]
  bank@meta.data[, grepl('^clust', names(bank@meta.data))] <- new.clusters

  return(bank)
}

#' @importFrom stats median
getSeed <- function(clusters, map.to) {
  if (!is.null(map.to)) {
    found <- map.to %in% names(clusters)
    if (!found) stop(map.to, ' not found in cluster names.')
    seed <- map.to
  } else {
    n.clust <- apply(clusters, 2, function(x) length(unique(x)))
    med <- median(n.clust)
    seed <- names(which.min(abs(n.clust - med)))
  }
  return(seed)
}

#' @importFrom plyr mapvalues
#' @importFrom RcppHungarian HungarianSolver
mapToSeed <- function(val, seed) {
  # ----- commenting out old version to test new. ------
#   con.mat <- table(val, seed)
#   cost.mat <- max(con.mat) - con.mat
#   matching <- HungarianSolver(cost.mat)$pairs

#   ## Mapping fr more clusters to fewer
#   unmapped <- which(matching[,2] == 0)
#   impute <- max(seed) + seq_len(length(unmapped))
#   matching[,2][matching[,2] == 0] <- impute

#   new.val <- mapvalues(x = val, from = matching[,1], to = matching[,2])
#   return(new.val)
  
  # ------ new version --------
   # seed and val must be numeric for this to work (because max(seed) is used for assigning unmatched values)
  con.mat <- table(val, seed)
  cost.mat <- max(con.mat) - con.mat
  matching <- HungarianSolver(cost.mat)$pairs
  
  # split into matched and unmatched
  matched <- matching[!(matching[,2] == 0),]
  unmateched <- matching[matching[,2] == 0,]
  ## Mapping fr more clusters to fewer
  unmapped <- which(matching[,2] == 0)
  impute <- max(seed) + seq_len(length(unmapped))
  unmateched[,2] <- impute
  
  # matched is currently the indices of the rows and cols of the contingency matrix
  matched.names.from <- c(as.numeric(rownames(cost.mat)[matched[,1]]), unmateched[,1])
  matched.names.to <- c(as.numeric(colnames(cost.mat)[matched[,2]]), unmateched[,2])
  
  new.val <- mapvalues(x = val, from = matched.names.from, to = matched.names.to)
  return(new.val)
}
                     

#' Calculate adjusted rand index for harmonized clusters
#'
#' @param bank BanksyObject with harmonized clusters
#' @param digits number of digits to round ARI to
#'
#' @return matrix of ARI
#'
#' @importFrom mclust adjustedRandIndex
#' @importFrom utils combn
#'
#' @export
getARI <- function(bank, digits = 3) {
  clust <- bank@meta.data[, grepl('^clust', names(bank@meta.data))]
  n.clust <- ncol(clust)
  comb <- combn(names(clust), 2)
  n.comb <- ncol(comb)
  ari <- numeric(length = ncol(comb))
  for (i in seq_len(n.comb)) {
    ari[i] <- adjustedRandIndex(clust[[comb[1,i]]], clust[[comb[2,i]]])
  }

  ari.mat <- diag(nrow = n.clust, ncol = n.clust)
  ari.mat[lower.tri(ari.mat)] <- ari
  rownames(ari.mat) <- colnames(ari.mat) <- colnames(clust)
  ari.mat <- ari.mat + t(ari.mat) - diag(nrow = n.clust, ncol = n.clust)
  ari.mat <- round(ari.mat, digits)
  return(ari.mat)
}

