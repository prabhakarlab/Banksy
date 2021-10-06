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

#' Harmonise cluster labels among parameter runs
#'
#' @param bank Banksy Object
#' @param mapto specify a cluster to map to
#' @param reverse reverse the map
#'
#' @importFrom mclust adjustedRandIndex
#' @importFrom plyr mapvalues
#' @importFrom RcppHungarian HungarianSolver
#' @importFrom stats median
#'
#' @return A list with BanksyObject with harmonised cluster labels and
#'   adjusted Rand indices for each cluster mapping
#'
#' @export
ConnectClusters <- function(bank, mapto = NULL, reverse = FALSE) {

  mdata <- bank@meta.data
  mnames <- names(mdata)
  if (is.null(mapto)) labels <- mnames[grepl('^res', mnames)]
  else labels <- mnames[grepl('^res', mnames) | grepl(mapto, mnames)]

  if (length(labels) == 1) {
    message('Only one cluster')
    return(bank)
  }

  clusters <- mdata[,labels]
  clustnames <- names(clusters)
  if (is.null(mapto)) {
    numClust <- apply(clusters, 2, function(x) length(unique(x)))
    medClust <- median(numClust)
    parent <- which.min(abs(numClust - medClust))
  } else {
    if (!(mapto %in% mnames)) stop('Specify a valid metadata column to map to.')
    parent <- which(labels == mapto)
  }

  message(paste0('Mapping to ', clustnames[parent]))
  children <- setdiff(seq_len(ncol(clusters)), parent)
  newLabels <- lapply(children, function(child, reverse) {

    contingency <- table(clusters[, parent], clusters[, child])
    if (reverse) {
      contingency <- t(contingency)
      mat <- max(contingency) - contingency
      matching <- HungarianSolver(mat)
      map <- matching$pairs[,2]
      map[which(map == 0)] <- apply(contingency[which(map == 0),],
                                    1, which.max)
      newChild <- mapvalues(clusters[, child],
                            from = matching$pairs[,1],
                            to = map,
                            warn_missing = FALSE)
    } else {
      mat <- max(contingency) - contingency
      matching <- HungarianSolver(mat)
      newChild <- mapvalues(clusters[, child],
                            from = matching$pairs[,2],
                            to = matching$pairs[,1],
                            warn_missing = FALSE)
    }
    adjRI <- adjustedRandIndex(clusters[, parent], clusters[, child])
    adjRI <- round(adjRI, 3)
    message(paste0('Mapped ', clustnames[child],
                   ' with adjusted Rand index ', adjRI))
    return(list(newChild, adjRI))
  }, reverse = reverse)

  result <- do.call(cbind.data.frame, lapply(newLabels, `[[`, 1))
  rand <- sapply(newLabels, `[[`, 2)
  names(result) <- names(clusters)[children]
  names(rand) <- names(clusters)[children]
  bank@meta.data <- cbind(bank@meta.data[,
                                         -which(names(bank@meta.data) %in% names(result))], result)
  return(list(BanksyObject = bank, rand.index = rand))
}
