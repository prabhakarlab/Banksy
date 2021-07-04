# NormalizeBanksy
# ScaleBanksy
# SubsetBanksy
# ComputeBanksy
# ClusterBanksy
# ConnectClusters

#' Normalize columns to a normFactor
#'
#' @param bank BanksyObject
#' @param normFactor normalization factor
#'
#' @return normalized Banksy object
#'
#' @export
NormalizeBanksy <- function(bank, normFactor = 100) {
  ## Own expression
  if (!is.null(bank@own.expr)) {
    if (is.list(bank@own.expr)) {
      bank@own.expr <- lapply(bank@own.expr, normalizer, normFactor)
    } else {
      bank@own.expr <- normalizer(bank@own.expr, normFactor)
    }
  }
  ## Nbr expression
  if (!is.null(bank@nbr.expr)) {
    if (is.list(bank@nbr.expr)) {
      bank@nbr.expr <- lapply(bank@nbr.expr, normalizer, normFactor)
    } else {
      bank@nbr.expr <- normalizer(bank@nbr.expr, normFactor)
    }
  }
  return(bank)
}

#' Scale rows
#'
#' @param bank BanksyObject
#' @param separate scale datasets separately
#'
#' @return scaled BanksyObject
#'
#' @export
ScaleBanksy <- function(bank, separate = TRUE) {
  if (!is.null(bank@own.expr)) {
    if (is.list(bank@own.expr)) {
      if (separate) bank@own.expr <- lapply(bank@own.expr, scaler)
      else {
        bank@own.expr <- scalerAll(bank@own.expr)
      }
    } else {
      bank@own.expr <- scaler(bank@own.expr)
    }
  }
  if (!is.null(bank@nbr.expr)) {
    if (is.list(bank@nbr.expr)) {
      bank@nbr.expr <- lapply(bank@nbr.expr, scaler)
    } else {
      bank@nbr.expr <- scaler(bank@nbr.expr)
    }
  }
  return(bank)
}

#' Subset a BanksyObject
#' @param x a BanksyObject
#' @param cells cells to filter by
#' @param dims dimensions to filter by - must correspond to valid columns
#' @param features genes to filter by
#' @param metadata metadata to filter by - must be valid colulmn
#' @param dataset dataset to subset dimensions by
#'
#' @importFrom rlang enquo quo_get_expr
#'
#' @return subset BanksyObject
#'
#' @export
SubsetBanksy <- function(
  x, cells=NULL, features=NULL, dims=TRUE, metadata=TRUE, dataset=NULL) {

  nfeaturesBef <- nrow(x@own.expr)
  ncellsBef <- ncol(x@own.expr)

  ## Enquote all input conditions
  dims <- rlang::enquo(dims)
  metadata <- rlang::enquo(metadata)

  ## Subset dimensions
  x <- subsetLocations(x, dims, dataset)

  ## Subset cells
  x <- subsetCells(x, cells)

  ## Subset metadata
  x@meta.data <- subset(x@meta.data, subset = eval(rlang::quo_get_expr(metadata)))
  ## Subset features
  x <- subsetFeatures(x, features)

  ## Consistent Cell filtering
  survivingCells <- subsetConsistent(x)
  x <- subsetCells(x, survivingCells)

  ## Filter cell locations, metadata and dim reductions based on surviving cells
  x@meta.data <- x@meta.data[x@meta.data$cell_ID %in% survivingCells,,drop=FALSE]
  x@dim.reduction <- lapply(x@dim.reduction, function(dim.red) {
    dim.red <- dim.red[rownames(dim.red) %in% survivingCells,,drop=FALSE]
    dim.red
  })
  if (is.list(x@own.expr)) {
    x@cell.locs <- lapply(x@cell.locs, function(x) {
      keep <- which(rownames(x) %in% survivingCells)
      x <- x[keep,,drop=FALSE]
    })
  } else {
    keep <- which(rownames(x@cell.locs) %in% survivingCells)
    x@cell.locs <- x@cell.locs[keep,,drop=FALSE]
  }

  return(x)
}

#' Compute Banksy Matrices
#' @param bank BanksyObject
#' @param sigma sigma
#' @param alpha alpha
#' @param kspatial kspatial
#' @param dimensions dims
#' @param spatialMode spatialMode
#' @param k_geom k_geom
#'
#' @importFrom data.table data.table setnames
#'
#' @return BanksyObject
#'
#' @export
ComputeBanksy <- function(bank,
                          ## For computing nbr matrix
                          sigma = 1.5, alpha = 0.05,
                          kspatial = 1000, dimensions = 'all',
                          spatialMode = 'kNN_r', k_geom = 10) {

  if (is.list(bank@own.expr)) {
    locs <- lapply(bank@cell.locs, function(x) {
      x <- data.table(x, keep.rownames = TRUE)
      setnames(x, 'rn', 'cell_ID')
    })
    nbr <- Map(function(expr, locs) {
      compute.banksyMatrices(expr, locs,
                             sigma = sigma,
                             alpha = alpha,
                             kspatial = kspatial,
                             dimensions = dimensions,
                             spatialMode = spatialMode,
                             k_geom = k_geom)},
      bank@own.expr, locs)
    bank@nbr.expr <- nbr
    names(bank@nbr.expr) <- names(bank@own.expr)

  } else {
    locs <- data.table(bank@cell.locs, keep.rownames = TRUE)
    setnames(locs, 'rn', 'cell_ID')
    nbr <- compute.banksyMatrices(bank@own.expr, locs,
                                  sigma=sigma,
                                  alpha=alpha,
                                  kspatial=kspatial,
                                  dimensions = dimensions,
                                  spatialMode = spatialMode ,
                                  k_geom = k_geom)
    bank@nbr.expr <- nbr
  }
  return(bank)
}


#' Cluster based on joint expression matrix
#'
#' @param bank BanksyObject
#' @param lambda lambda mixing
#' @param resolution resolution
#' @param kneighbours kneighbours (sNN)
#' @param npcs num princ. comp.
#' @param nneighbors umap
#' @param spread umap
#' @param nepochs umap
#' @param mindist umap
#' @param leideniters leiden iters
#' @param verbose messages
#'
#' @importFrom irlba prcomp_irlba
#' @importFrom uwot umap
#' @importFrom dbscan kNN sNN
#' @importFrom leidenAlg leiden.community
#' @importFrom igraph graph_from_data_frame
#' @importFrom data.table setorder data.table
#'
#' @return BanksyObject with cluster labels in meta.data
#'
#' @export
ClusterBanksy <- function(bank,
                          ## Grid parameters
                          lambda = 0.25, resolution = 0.8, kneighbours = 40,
                          ## PCA
                          npcs = 50,
                          ## UMAP parameters
                          nneighbors = 30, spread = 3,
                          nepochs = 300, mindist = 0.3,
                          ## Leiden parameters
                          leideniters = -1,
                          verbose = FALSE) {

  max_iters <- prod(length(lambda), length(resolution), length(kneighbours))
  iter = 1
  message(paste0('Iteration ', iter, ' out of ', max_iters) )
  for (lam in lambda) {

    joint <- getBanksyMatrix(bank)$expr

    if (verbose) message('Running PCA')
    pca <- prcomp_irlba(t(joint), n = npcs)$x
    rownames(pca) <- bank@meta.data$cell_ID
    bank@dim.reduction[[paste0('pca_', lam)]] <- pca

    if (verbose) message('Computing UMAP')
    umap <- umap(pca, n_neighbors = nneighbors, spread = spread,
                 n_epochs = nepochs, min_dist = mindist)
    rownames(umap) <- bank@meta.data$cell_ID
    bank@dim.reduction[[paste0('umap_', lam)]] <- umap

    for (res in resolution) {
      for (knb in kneighbours) {

        ## Giotto implementation
        if (verbose) message(paste0('Computing sNN'))
        knet <- kNN(pca, k = knb)
        snet <- sNN(knet, k = knb)
        from <- shared <- .N <- NULL
        snet_dt = data.table(from = rep(1:nrow(snet$id), knb),
                             to = as.vector(snet$id),
                             weight = 1/(1 + as.vector(snet$dist)),
                             distance = as.vector(snet$dist),
                             shared = as.vector(snet$shared))
        data.table::setorder(snet_dt, from, -shared)
        snet_dt[, rank := 1:.N, by = from]
        snet_dt <- snet_dt[rank <= 3 | shared >= 5]
        graph <- graph_from_data_frame(snet_dt)

        if (verbose) message('Running Leiden clustering')
        clusLeiden <- as.numeric(leiden.community(graph,
            resolution = res, n.iterations = leideniters)$membership)
        clusName <- paste0('res',res,'_lam',lam, '_k', knb)
        bank@meta.data[[clusName]] <- clusLeiden


        message(paste0('Finished clustering for Lambda=', lam,
                       ', Resolution=', res, ', K Neighbours=', knb))
        iter <- iter + 1
        if (iter <= max_iters) message(paste0('Iteration ', iter,
                                              ' out of ', max_iters) )

      }
    }
  }
  return(bank)

}


#' Harmonise cluster labels among parameter runs
#'
#' @param bank Banksy Object
#' @param verbose verbose or not
#' @param optim optimise or not (KS-test)
#'
#' @importFrom data.table copy
#' @importFrom stats median ks.test
#' @importFrom plyr mapvalues
#'
#' @return BanksyObject with harmonised cluster labels
#'
#' @export
ConnectClusters <- function(bank, verbose=FALSE, optim=TRUE) {

  d <- bank@meta.data
  clust <- d[,grep('^res', colnames(d)),drop=FALSE]
  clustNames <- names(clust)

  if (length(clustNames) == 1) {
    message("Only one cluster.")
    return(bank)
  }

  ## Init the new clustering output
  newClust <- copy(clust)
  ## Use median clusters as seed
  numClust <- apply(clust, 2, function(x) length(unique(x)))
  medClust <- median(numClust)
  ## Get the clusters for each parameter run
  allClust <- apply(clust, 2, unique)
  ## Settle the parent cluster labels
  #parent <- max(which(numClust == medClust))
  parent <- which.min(abs(numClust - medClust))
  parentClusters <- allClust[[parent]]
  newClust[,parent] <- plyr::mapvalues(clust[,parent],
                                  from = parentClusters,
                                  to = seq_len(length(unique(parentClusters))),
                                  warn_missing = FALSE)
  parentDist <- as.numeric(table(newClust[,parent]))
  message(paste0('Mapping clusterings to ', clustNames[parent]))
  ## The rest will be children
  children <- setdiff(seq_len(ncol(clust)), parent)
  ## Iterate over the clusters of the parent
  for (child in children) {

    message(paste0('Processing ', clustNames[child]))
    ## Child-centered mapping approach
    childClust <- sort(allClust[[child]], decreasing = FALSE)
    childDist <- rep(0, medClust)

    for (i in childClust) {

      map <- factor(newClust[,parent][clust[,child]==i])
      tab <- tabulate(map)
      hit <- as.numeric(levels(map)[tab == max(tab)])

      ## Optimize
      if (optim & i > length(childClust) / 2 & length(tab) > 1) {
        topn <- sort(tab, decreasing = TRUE)[seq_len(2)]
        ambig <- min(topn)/max(topn) > 0.9
        if (ambig) {
          if (verbose) message(paste0('Ambiguous mapping for cluster ', i,
                                      ' - using KS-test'))
          hits <- as.numeric(levels(map)[tab %in% topn])
          test1 <- childDist + init(medClust, hits[1],
                                     tab[which(levels(map)==hits[1])])
          test2 <- childDist + init(medClust, hits[2],
                                     tab[which(levels(map)==hits[2])])
          stat1 <- ks.test(parentDist, jitter(test1), exact=FALSE)$stat
          stat2 <- ks.test(parentDist, jitter(test2), exact=FALSE)$stat
          hit <- ifelse(stat1 < stat2, hits[1], hits[2])
        }
      } else {
        hit <- max(hit)
      }
      newClust[,child][clust[,child]==i] <- hit
      childDist[hit] <- childDist[hit] + tab[hit]
      childDist[is.na(childDist)] <- 0
      if(verbose) message(paste0('Mapping child cluster ',
                                 i , ' to parent cluster ', hit))
    }
  }
  maxClust <- max(apply(newClust, 2, max))
  fromMap <- seq_len(maxClust)
  toMap <- getPalette(maxClust)
  bank@meta.data <- cbind(cell_ID = bank@meta.data$cell_ID, newClust)
  return(bank)
}

#' Returns the Banksy matrix (own + nbr)
#'
#' @param bank Banksy Object
#' @param lambda lambda
#'
#' @return BanksyMatrix
#'
#' @export
getBanksyMatrix <- function(bank, lambda = 0.25) {

  # Optimize this
  if (is.list(bank@own.expr)) {
    own <- do.call(cbind, bank@own.expr)
    nbr <- do.call(cbind, bank@nbr.expr)
    joint <- rbind(sqrt(1-lambda)*own, sqrt(lambda)*nbr)

    locs <- do.call(rbind, bank@cell.locs)
    rownames(locs) <- colnames(joint)

  } else {
    joint <- rbind(sqrt(1-lambda)*bank@own.expr,
                   sqrt(lambda)*bank@nbr.expr)
    locs <- bank@cell.locs
  }
  return(list(expr = joint,
              locs = locs))
}
