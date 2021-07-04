# Index of helper functions and where they are called
# normalizer ---------------- NormalizeBanksy
# scaler -------------------- ScaleBanksy
# scalerAll ----------------- ScaleBanksy
# geneFilter ---------------- BanksyObject
# init ---------------------- ConnectClusters
# compute.banksyMatrices ---- ComputeBanksy
# getPalette ---------------- Plotting
# subsetDims ---------------- SubsetBanksy
# subsetLocations ----------- SubsetBanksy
# subsetCells --------------- SubsetBanksy
# subsetFeatures ------------ SubsetBanksy
# subsetConsistent ---------- SubsetBanksy
# getAssay ------------------ Plotting
# getHeatmapPalette --------- Plotting

#' @importFrom collapse `%c/%`
#' @importFrom matrixStats colSums2
normalizer <- function(x, normFactor) {
    x <- t(t(x) %c/% colSums2(x)) * normFactor
    return(x)
}

#' @importFrom matrixStats rowMeans2 rowSds
scaler <- function(x) {
  rm <- rowMeans2(x)
  rsd <- rowSds(x)
  x <- (x - rm) / rsd
  x[is.nan(x)] <- 0
  return(x)
}

# Scale rows of the concatenate of matrices without concatenating them
#' @importFrom matrixStats rowSums2 rowVars rowMeans2
scalerAll <- function(x) {

  sumIndiv <- lapply(x, rowSums2)
  nIndiv <- lapply(x, ncol)
  nAll <- sum(unlist(nIndiv))
  sumAll <- Reduce(`+`, sumIndiv)
  meanAll <- sumAll / nAll

  varIndiv <- lapply(x, rowVars)
  meanIndiv <- lapply(x, rowMeans2)
  zeroIndiv <- lapply(meanIndiv, `-`, meanAll)
  a <- Map(`*`, varIndiv, nIndiv)
  b <- Map(`*`, zeroIndiv, nIndiv)
  numIndiv <- Map(`+`, a,b)
  numAll <- Reduce(`+`, numIndiv)
  sdAll <- sqrt(numAll / nAll)

  x <- lapply(x, function(d) {
    d <- (d - meanAll) / sdAll
    d[is.nan(d)] <- 0
    d
  })
  return(x)
}

geneFilter <- function(x, genes.filter, min.cells.expressed) {

    ngenesBef <- sapply(x, function(x) dim(x)[1])

    if (genes.filter == 'union') {
        all.genes <- Reduce(union, lapply(x, rownames))
        x <- lapply(x, function(x) {
        genes <- setdiff(all.genes, rownames(x))
          if (length(genes) == 0) {
              return(x)
          } else {
          append <- matrix(0, nrow = length(genes), ncol = ncol(x))
          colnames(append) <- colnames(x)
          rownames(append) <- genes
          x <- rbind(x, append)
          return(x)
          }
        })
    } else if (genes.filter == 'intersect') {
        common.genes <- Reduce(intersect, lapply(x, rownames))
        x <- lapply(x, function(x) x[rownames(x) %in% common.genes,])

        if (min.cells.expressed > 0) {
          message(paste0('Filering genes expressed in less than ', min.cells.expressed, ' cells'))
          pass.genes <- lapply(x, function(x) rownames(x)[rowSums(x > 0) >= min.cells.expressed])
          pass.genes <- Reduce(intersect, pass.genes)
          x <- lapply(x, function(x) x[rownames(x) %in% pass.genes,])
          ngenesAft <- sapply(x, function(x) dim(x)[1])
          filt <- ngenesBef - ngenesAft
          for (i in seq_len(length(x))) {
            message(paste0('Filtered ', filt[i], ' genes from dataset ',
                           names(x)[i]))
          }

        }
    }

    ## Harmonise gene name orderings
    gene.names <- sort(rownames(x[[1]]))
    x <- lapply(x, function(x) {
      x[match(gene.names, rownames(x)),]
    })

    return(x)
}

init <- function(n, idx, val) {
  x <- rep(0, n)
  x[idx] <- val
  x
}


#' @importFrom pals kelly glasbey polychrome
getPalette <- function(n) {
  all.cols <- as.character(pals::kelly()[-c(1,3)],
                           pals::glasbey(),
                           pals::polychrome())
  return(all.cols[seq_len(n)])
}

#' @importFrom rlang quo_get_expr
subsetDims <- function(locs, dims) {
  locs <- subset(data.frame(locs), subset = eval(rlang::quo_get_expr(dims)))
  return(locs)
}

subsetLocations <- function(x, dims=TRUE, dataset=NULL) {

  if (is.list(x@own.expr)) {

    if (is.null(dataset)) {
      message('Dataset not specified. Subsetting all data sets by dimension.')
      x@cell.locs <- lapply(x@cell.locs, function(x) {
        subsetDims(x, dims)
      })
    } else {
      if (!(dataset %in% names(x@cell.locs))) {
        stop(paste0('Dataset ', dataset, ' not found. One of ',
                    paste(names(x@cell.locs), collapse = ' ')))
      }
      assayid <- which(names(x@cell.locs) == dataset)
      assay <- x@cell.locs[[assayid]]
      x@cell.locs[[assayid]] <- subsetDims(assay, dims)
    }

  } else {
    x@cell.locs <- subsetDims(x@cell.locs, dims)
  }

  return(x)
}

subsetCells <- function(x, cells) {

  if (is.null(cells)) {
    return(x)
  }

  if (is.list(x@own.expr)) {
    x@own.expr <- lapply(x@own.expr, function(x) {
      keep <- which(colnames(x) %in% cells)
      x <- x[,keep, drop = FALSE]
      x
    })
    if (!is.null(x@nbr.expr)) {
      x@nbr.expr <- lapply(x@nbr.expr, function(x) {
        keep <- which(colnames(x) %in% cells)
        x <- x[,keep, drop = FALSE]
        x
      })
    }
  } else {
    keep <- which(colnames(x@own.expr) %in% cells)
    x@own.expr <- x@own.expr[, keep, drop = FALSE]
    if (!is.null(x@nbr.expr)) {
      x@nbr.expr <- x@nbr.expr[, keep, drop = FALSE]
    }
  }
  return(x)
}

subsetFeatures <- function(x, features) {

  if (is.null(features)) {
    return(x)
  }

  if (is.list(x@own.expr)) {
    keep <- which(rownames(x@own.expr[[1]]) %in% features )
    x@own.expr <- lapply(x@own.expr, function(x) {
      x <- x[keep,,drop = FALSE]
      x
    })
    if (!is.null(x@nbr.expr)) {
      x@nbr.expr <- lapply(x@nbr.expr, function(x) {
        x <- x[keep,,drop=FALSE]
        x
      })
    }

  } else {
    keep <- which(rownames(x@own.expr) %in% features)
    x@own.expr <- x@own.expr[keep,,drop=FALSE]
    if (!is.null(x@nbr.expr)) {
      x@nbr.expr <- x@nbr.expr[keep,,drop=FALSE]
    }
  }
  return(x)
}

subsetConsistent <- function(x) {

  cellMeta <- x@meta.data$cell_ID

  if (is.list(x@own.expr)) {
    cellCell <- as.character(unlist(lapply(x@own.expr, colnames)))
    cellLocs <- as.character(unlist(lapply(x@cell.locs, rownames)))
  } else {
    cellCell <- colnames(x@own.expr)
    cellLocs <- rownames(x@cell.locs)
  }

  cellOut <- intersect(cellMeta, intersect(cellCell, cellLocs))
  return(cellOut)
}

getAssay <- function(bank, assay, dataset) {

  if (!is.na(pmatch(assay, c('own.expr', 'nbr.expr', 'custom.expr')))) {
    slot <- get(assay, mode = 'function')
    mat <- slot(bank)
    if (is.list(mat)) {
      if (is.null(dataset)) {
        message('Dataset not specified. Choosing first dataset.')
        mat <- mat[[1]]
      } else if (!(dataset %in% names(mat))) {
        stop(paste0('Dataset ', dataset, ' not found.'))
      } else {
        mat <- mat[[dataset]]
      }
    }
    if (is.null(mat)) stop('Assay chosen is empty.')
  } else {
    stop('Specify a valid assay.')
  }
  return(mat)
}

#' @importFrom circlize colorRamp2
#' @importFrom stats quantile
getHeatmapPaletteold <- function(mat, col, col.breaks) {

  if (is.null(col)) {
    lim.high <- quantile(mat, 0.9)
    lim.low <- quantile(mat, 0.1)
    col.fun <- colorRamp2(c(lim.low, 0, lim.high),
                          c('blue', 'white', 'red'))
  } else {
    if (is.null(col.breaks)) {
      lim.high <- quantile(mat, 0.9)
      lim.low <- quantile(mat, 0.1)
      col.breaks <- seq(lim.low, lim.high, length.out = length(col))
    } else if (length(col.breaks) != length(col)) {
      stop ('Unequal number of breaks and colors.')
    }
    col.fun <- colorRamp2(col.breaks, col)
  }
  return(col.fun)
}

#' @importFrom circlize colorRamp2
#' @importFrom stats quantile
getHeatmapPalette <- function(mat, col, col.breaks) {

  if (!is.null(col) &
      !is.null(col.breaks) &
      (length(col) != length(col.breaks))) {
    stop('Unequal number of colours and breaks')
  }
  if (is.null(col)) {
    col <- c('blue', 'white', 'red')
  }
  if (is.null(col.breaks)) {
    lim.high <- quantile(mat, 0.9)
    lim.low <- quantile(mat, 0.1)
    col.breaks <- seq(lim.low, lim.high, length.out = length(col))
  }

  col.fun <- colorRamp2(col.breaks, col)

  return(col.fun)
}

getSpatialDims <- function(locs, dimensions, alpha) {

  cellID <- seq_len(nrow(locs))
  names(cellID) <- locs$cell_ID
  locs <- locs[, grepl('sdim', colnames(locs)), with = FALSE]

  if (dimensions != "all") {
    if (all(dimensions %in% colnames(locs))) stop('Invalid dimensions specified')
    locs <- locs[, dimensions]
  }

  # Set kernel radius based on number of dimensions
  ndims <- ncol(locs)
  nz <- length(unique(locs$dimz))
  if (ndims < 2) {
    stop('Must have at least 2 spatial locations.')
  } else if (nz == 1) {
    ndims <- ndims - 1
  }
  kernelRadius <- sqrt(-ndims * log(alpha))

  locs <- as.matrix(locs)
  rownames(locs) <- names(cellID)

  return(list(locs = locs, kr = kernelRadius, cellID = cellID))
}

#' @importFrom dbscan kNN
#' @importFrom data.table data.table setDT  setnames `:=` rbindlist
#' @importFrom stats dnorm
withRNNgauss <- function(locs, sigma, kspatial, kernelRadius, verbose) {

  if (verbose) message('Computing Banksy matrix')
  if (verbose) message('Spatial mode is rNN gaussian')
  if (verbose) message(paste0('Parameters: sigma = ', sigma, ', kspatial = ', kspatial))

  tryCatch({
    knn <- dbscan::kNN(x = locs, k = kspatial)
  },
  error=function(cond) {
    message(paste("Not enough neighbours at kspatial = ", kspatial, " level."))
    message(cond)
  })

  medianDist = median(as.vector(knn$dist[,1]))
  knnDF <- data.table(from = rep(seq_len(nrow(knn$id)), kspatial),
                      to = as.vector(knn$id),
                      weight = dnorm(as.vector(knn$dist), mean = 0,
                                     sd = medianDist * sigma),
                      distance = as.vector(knn$dist))

  distance <- norm.weight <- weight <- from <- to <- NULL
  knnDF <- knnDF[distance < sigma * kernelRadius * medianDist, ]
  setDT(knnDF)[, norm.weight := weight / sum(weight), by = from]
  knnDF <- knnDF[, -3, with = FALSE]

  ## Create dummy entries for filtered out cells
  iso <- setdiff(seq_len(nrow(locs)), unique(knnDF$from))
  isomat <- c(rep(iso, 2), rep(0, 2*length(iso)))
  isomat <- data.table(matrix(isomat, ncol = ncol(knnDF)))
  knnDF <- rbindlist(list(knnDF, isomat), use.names = FALSE)
  knnDF <- knnDF[order(from, to)]

  setnames(knnDF, 'norm.weight', 'weight')

  return(knnDF)
}

#' @importFrom dbscan kNN
#' @importFrom data.table data.table setnames `:=`
withKNNrank <- function(locs, k_geom, verbose) {

  if (verbose) message('Computing Banksy matrix')
  if (verbose) message('Spatial mode is kNNrank')
  if (verbose) message(paste0('Parameters: k_geom = ', k_geom))

  tryCatch({
    knn <- dbscan::kNN(x = locs, k = k_geom)
    unnormWt <- exp(-seq(1, k_geom, 1)^2 / (2*(k_geom/1.5)^2))
    normWt <- unnormWt / sum(unnormWt)

    weightMatrix = t(matrix(normWt, nrow = k_geom, ncol = nrow(knn$id)))
    knnDF = data.table(from = rep(seq_len(nrow(knn$id)) , k_geom),
                       to = as.vector(knn$id),
                       weight = as.vector(weightMatrix),
                       distance = as.vector(knn$dist))
  },
  error=function(cond) {
    message(paste("Not enough neighbours at kspatial = ", k_geom, " level."))
    message(cond)
  })

  return(knnDF)
}

#' @importFrom dbscan kNN
#' @importFrom data.table data.table setnames `:=`
withKNNr <- function(locs, k_geom, verbose) {

  if (verbose) message('Computing Banksy matrix')
  if (verbose) message('Spatial mode is kNNr')
  if (verbose) message(paste0('Parameters: k_geom = ', k_geom))

  tryCatch({
    knn <- dbscan::kNN(x = locs, k = k_geom)
  },
  error=function(cond) {
    message(paste("Not enough neighbours at k_geom = ", k_geom, " level."))
    message(cond)
  })

  norm.weight <- weight <- from <- NULL
  knnDF <- data.table(from = rep(seq_len(nrow(knn$id)), k_geom),
                      to = as.vector(knn$id),
                      weight = 1/as.vector(knn$dist),
                      distance = as.vector(knn$dist))
  knnDF[, norm.weight := weight / sum(weight), by = from]
  knnDF <- knnDF[,-3, with=FALSE]
  setnames(knnDF, 'norm.weight', 'weight')

  return(knnDF)
}

#' @importFrom dbscan kNN
#' @importFrom data.table data.table setnames `:=`
withKNNunif <- function(locs, k_geom, verbose) {

  if (verbose) message('Computing Banksy matrix')
  if (verbose) message('Spatial mode is kNNunif')
  if (verbose) message(paste0('Parameters: k_geom = ', k_geom))

  tryCatch({
    knn <- dbscan::kNN(x = locs, k = k_geom)
  },
  error=function(cond) {
    message(paste("Not enough neighbours at kspatial = ", k_geom, " level."))
    message(cond)
  })

  norm.weight <- weight <- from <- NULL
  knnDF <- data.table(from = rep(seq_len(nrow(knn$id)), k_geom),
                      to = as.vector(knn$id),
                      weight = 1,
                      distance = as.vector(knn$dist))
  knnDF[, norm.weight := weight / sum(weight), by = from]
  knnDF <- knnDF[,-3,with=FALSE]
  setnames(knnDF, 'norm.weight', 'weight')

  return(knnDF)
}

#' @importFrom zeallot `%<-%`
compute.banksyMatrices <- function(gcm, locs,
                                    sigma=3, alpha=0.05,
                                    kspatial=1000, k_geom = 10,
                                    dimensions = 'all',
                                    spatialMode = 'kNN_r',
                                    verbose = FALSE){

  # Global binding
  weight <- from <- to <- kernelRadius <- cellID <- NULL
  c(locs, kernelRadius, cellID) %<-% getSpatialDims(locs, dimensions, alpha)

  message('Computing neighbors...')
  if (spatialMode == 'rNN_gauss'){
    knnDF <- withRNNgauss(locs = locs, sigma = sigma, kspatial = kspatial,
                          kernelRadius = kernelRadius, verbose = verbose)

  } else if (spatialMode == 'kNN_rank' ){
    knnDF <- withKNNrank(locs = locs, k_geom = k_geom, verbose = verbose)

  } else if (spatialMode =='kNN_r'){
    knnDF <- withKNNr(locs = locs, k_geom = k_geom, verbose = verbose)

  } else if (spatialMode == 'kNN_unif') {
    knnDF <- withKNNunif(locs = locs, k_geom = k_geom, verbose = verbose)

  } else {
    stop('Invalid spatialMode. One of rNN_gauss, kNN_rank, kNN_r, kNN_unif')
  }

  message('Computing neighbor matrix...')
  aggr <- knnDF[, gcm[,to, drop = FALSE] %*% weight, by = from]
  ncm <- matrix(aggr$V1, nrow = nrow(gcm), ncol = ncol(gcm))
  rownames(ncm) <- paste0(rownames(gcm), '.nbr')
  colnames(ncm) <- colnames(gcm)
  message('Done')

  return(ncm)
}





















