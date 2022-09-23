# NormalizeBanksy
# ScaleBanksy
# SubsetBanksy
# ComputeBanksy

#' Normalize columns to a normFactor
#'
#' @param bank BanksyObject
#' @param assay (character) assay to scale
#' \itemize{
#'  \item{own: Scales own.expr}
#'  \item{nbr: Scales nbr.expr}
#'  \item{both: Scales own.expr and nbr.expr (default)}
#' }
#' @param norm_factor (numeric) normalization factor (default: 100)
#' @param log_norm (logical) log transforms data (default: FALSE)
#' @param pseudocount (numeric) pseudocount for log transform (default: 0.1)
#' @param base (numeric) base for log transform (default: 10)
#'
#' @return normalized Banksy object
#'
#' @export
#' 
#' @examples
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' # Normalize the own.expr matrix
#' bank <- NormalizeBanksy(bank)
#' 
NormalizeBanksy <- function(bank, assay = 'both', norm_factor = 100,
                            log_norm = FALSE, pseudocount = 0.1, base = 10) {

  scaleOwn <- TRUE
  scaleNbr <- TRUE
  if (assay == 'own') {
    scaleNbr <- FALSE
  } else if (assay == 'nbr') {
    scaleOwn <- FALSE
  } else if (assay != 'both'){
    stop('Specify a valid assay. One of both, own, or nbr.')
  }

  if (!is.null(bank@own.expr) & scaleOwn) {
    if (is.list(bank@own.expr)) {
      bank@own.expr <- lapply(bank@own.expr, normalizer, norm_factor,
                              log_norm, pseudocount, base)
    } else {
      bank@own.expr <- normalizer(bank@own.expr, norm_factor,
                                  log_norm, pseudocount, base)
    }
  }
  
  if (!is.null(bank@nbr.expr) & scaleNbr) {
    if (is.list(bank@nbr.expr)) {
      bank@nbr.expr <- lapply(bank@nbr.expr, normalizer, norm_factor,
                              log_norm, pseudocount, base)
      bank@harmonics <- lapply(names(bank@harmonics), function(x) {
          lapply(bank@harmonics[[x]], function(x) {
              normalizer(x, norm_factor, log_norm, pseudocount, base)
          })
      })
      names(bank@harmonics) <- names(bank@nbr.expr)
    } else {
      bank@nbr.expr <- normalizer(bank@nbr.expr, norm_factor,
                                  log_norm, pseudocount, base)
      bank@harmonics <- lapply(bank@harmonics, function(x) {
          normalizer(x, norm_factor, log_norm, pseudocount, base)
      })
    }
  }
  return(bank)
}

#' Scale rows
#'
#' @param bank BanksyObject
#' @param assay (character) assay to scale
#' \itemize{
#'  \item{own: Scales own.expr}
#'  \item{nbr: Scales nbr.expr}
#'  \item{both: Scales own.expr and nbr.expr (default)}
#' }
#' @param separate (logical) scale datasets separately
#'
#' @return scaled BanksyObject
#'
#' @export
#' 
#' @examples 
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' # Normalize the own.expr matrix
#' bank <- NormalizeBanksy(bank)
#' # Scale the own.expr matrix
#' bank <- ScaleBanksy(bank)
#' 
ScaleBanksy <- function(bank, assay = 'both', separate = TRUE) {

  scaleOwn <- TRUE
  scaleNbr <- TRUE
  if (assay == 'own') {
    scaleNbr <- FALSE
  } else if (assay == 'nbr') {
    scaleOwn <- FALSE
  } else if (assay != 'both'){
    stop('Specify a valid assay. One of both, own, or nbr.')
  }

  if (!is.null(bank@own.expr) & scaleOwn) {
    if (is.list(bank@own.expr)) {
      if (separate) bank@own.expr <- lapply(bank@own.expr, scaler)
      else {
        bank@own.expr <- scalerAll(bank@own.expr)
      }
    } else {
      bank@own.expr <- scaler(bank@own.expr)
    }
  }
  if (!is.null(bank@nbr.expr) & scaleNbr) {
    if (is.list(bank@nbr.expr)) {
      bank@nbr.expr <- lapply(bank@nbr.expr, scaler)
      bank@harmonics <- lapply(names(bank@harmonics), function(x) {
          lapply(bank@harmonics[[x]], scaler)
      })
      names(bank@harmonics) <- names(bank@nbr.expr)
    } else {
      bank@nbr.expr <- scaler(bank@nbr.expr)
      bank@harmonics <- lapply(bank@harmonics, scaler)
    }
  }
  return(bank)
}

#' Compute Banksy Matrices
#' @param bank BanksyObject
#' @param K (numeric) compute up to the k-th azimuthal fourier harmonic (default: 1) 
#' @param spatial_mode (character) 
#' \itemize{
#'  \item{kNN_r: k-nearest neighbors with $1/r$ kernel}
#'  \item{kNN_rn: k-nearest neighbors with $1/r^n$ kernel}
#'  \item{kNN_rank: k-nearest neighbors with rank Gaussian kernel}
#'  \item{kNN_unif: k-nearest neighbors wth uniform kernel}
#'  \item{kNN_median: k-nearest neighbors with median-scaled Gaussian kernel (default)}
#'  \item{rNN_gauss: radial nearest neighbors with Gaussian kernel}
#' }
#' @param k_geom (numeric) number of neighbors to use (for kNN)
#' @param n (numeric) exponent of radius (for kNN_rn)
#' @param sigma (numeric) std. dev. of Gaussian kernel (for rNN_gauss)
#' @param alpha (numeric) determines radius used: larger alphas give
#'   smaller radii (for rNN_gauss)
#' @param k_spatial (numeric) initial number of neighbors to use (for rNN_gauss)
#' @param dimensions (character vector) dimensions to use when computing neighborhood
#' \itemize{
#'  \item{subset of colnames of cell.locs}
#'  \item{all}{Uses all colnames of cell.locs to compute (default)}
#' }
#' @param center (logical) center higher order harmonics in local neighborhoods
#' @param verbose messages
#'
#' @importFrom data.table data.table setnames
#'
#' @return BanksyObject
#'
#' @export
#' 
#' @examples 
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' bank <- NormalizeBanksy(bank)
#' # Compute neighbors 
#' bank <- ComputeBanksy(bank)
#' 
ComputeBanksy <- function(bank, K = 1,
                          spatial_mode = 'kNN_median', k_geom = 10, n = 2,
                          sigma = 1.5, alpha = 0.05, k_spatial = 100,
                          dimensions = 'all', center = TRUE, verbose=TRUE) {
    
    K <- seq(0, K)
    if (is.list(bank@own.expr)) {
        # Multi-dataset case
        locs <- lapply(bank@cell.locs, function(x) {
            x <- data.table(x, keep.rownames = TRUE)
            setnames(x, 'rn', 'cell_ID')
        })
        knn_lst <- lapply(locs, function(dlocs) {
            computeNeighbors(dlocs,
                             spatial_mode = spatial_mode, k_geom = k_geom, n = n,
                             sigma=sigma, alpha=alpha, k_spatial=k_spatial,
                             dimensions = dimensions, verbose = verbose)
        })
        nbr_lst <- Map(function(expr, knn_df) {
            computeHarmonics(expr, knn_df, K = 0, center = FALSE)
        }, bank@own.expr, knn_lst)
        
        har_lst <- Map(function(expr, knn_df) {
            har <- lapply(setdiff(K, 0), function(k) {
                computeHarmonics(expr, knn_df, K = k, center = center)
            })
            if (length(har) > 0) names(har) <- paste0('k', setdiff(K, 0))
            har
        }, bank@own.expr, knn_lst)
        
        bank@nbr.expr <- nbr_lst
        bank@harmonics <- har_lst
        names(bank@nbr.expr) <- names(bank@own.expr)
        
    } else {
        # Single dataset case
        locs <- data.table(bank@cell.locs, keep.rownames = TRUE)
        setnames(locs, 'rn', 'cell_ID')
        knn_df <- computeNeighbors(locs,
                                   spatial_mode = spatial_mode, k_geom = k_geom, n = n,
                                   sigma=sigma, alpha=alpha, k_spatial=k_spatial,
                                   dimensions = dimensions, verbose = verbose)
        nbr <- computeHarmonics(bank@own.expr, knn_df, K = 0, center = FALSE)
        har <- lapply(setdiff(K, 0), function(k) {
            computeHarmonics(bank@own.expr, knn_df, K = k, center = center)
        })
        if (length(har) > 0) names(har) <- paste0('k', setdiff(K, 0))
        bank@nbr.expr <- nbr
        bank@harmonics <- har
        
    }
    
    return(bank)
}

# Index of helper functions and where they are called
# geneFilter ---------------- BanksyObject
# normalizer ---------------- NormalizeBanksy
# scaler -------------------- ScaleBanksy
# scalerAll ----------------- ScaleBanksy
# computeHarmonics ---------- ComputeBanksy
# computeNeighbors ---------- ComputeBanksy
# getSpatialDims ------------ ComputeBanksy
# with<Method> -------------- ComputeBanksy
# getPhi -------------------- computeHarmonics
# mmult --------------------- ComputeBanksy
# getLambdas ---------------- getBanksyMatrix

geneFilter <- function(x, genes.filter, min.cells.expressed) {

  ngenesBef <- vapply(x, function(x) dim(x)[1], FUN.VALUE = numeric(1))

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
      message('Filering genes expressed in less than ', min.cells.expressed, ' cells')
      pass.genes <- lapply(x, function(x) rownames(x)[rowSums(x > 0) >= min.cells.expressed])
      pass.genes <- Reduce(intersect, pass.genes)
      x <- lapply(x, function(x) x[rownames(x) %in% pass.genes,])
      ngenesAft <- vapply(x, function(x) dim(x)[1], FUN.VALUE = numeric(1))
      filt <- ngenesBef - ngenesAft
      for (i in seq_len(length(x))) {
        message('Filtered ', filt[i], ' genes from dataset ', names(x)[i])
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

normalizer <- function(x, norm_factor, log_norm, pseudocount, base) {
  x <- t(t(x) / colSums(x)) * norm_factor
  if (log_norm) x <- log(x + pseudocount, base = base)
  return(x)
}

#' @importFrom matrixStats rowSds
scaler <- function(x) {
  rm <- rowMeans(x)
  rsd <- rowSds(x)
  x <- (x - rm) / rsd
  x[is.nan(x)] <- 0
  return(x)
}

fscale <- function(x) {
    rm <- rowMeans(x)
    x <- (x - rm) 
    return(x)
}

#' @importFrom matrixStats rowVars
scalerAll <- function(x) {

  sumIndiv <- lapply(x, rowSums)
  nIndiv <- lapply(x, ncol)
  nAll <- sum(unlist(nIndiv))
  sumAll <- Reduce(`+`, sumIndiv)
  meanAll <- sumAll / nAll

  varIndiv <- lapply(x, rowVars)
  meanIndiv <- lapply(x, rowMeans)
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

getLambdasDeprecate <- function(lambda, n_harmonics) {
    lam = c(lambda, lambda * (2^-seq_len(n_harmonics-1)))
    lam = c(1 - sum(lam), lam)
    message('Squared lambdas: ', paste0(lam, collapse = ', '))
    sqrt(lam)
}

getLambdas <- function(lambda, n_harmonics) {
    weights = lambda * (2^-seq(0, n_harmonics-1))
    weights = weights / sum(2^-seq(0, n_harmonics-1))
    lam = c(1 - sum(weights), weights)
    message('Squared lambdas: ', paste0(round(lam,4), collapse = ', '))
    sqrt(lam)
}

#' Returns the Banksy matrix (own + nbr)
#'
#' @param bank Banksy Object
#' @param lambda (numeric) spatial weighting parameter
#' @param K (numeric) compute up to the k-th azimuthal fourier harmonic (default: 1) 
#'
#' @return BanksyMatrix
#'
#' @export
#' 
#' @examples 
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' # Normalize the own.expr matrix
#' bank <- NormalizeBanksy(bank)
#' # Compute BANKSY matrix
#' bank <- ComputeBanksy(bank)
#' bm <- getBanksyMatrix(bank)
getBanksyMatrix <- function(bank, lambda = 0.15, K = 1) {
    
    # Optimize this
    if (is.list(bank@own.expr)) {
        # Multiple dataset case
        own <- do.call(cbind, bank@own.expr)
        nbr <- do.call(cbind, bank@nbr.expr)
        harmonics <- names(bank@harmonics[[1]])
        out <- lapply(harmonics, function(k) {
            do.call(cbind, lapply(bank@harmonics, function(x) x[[k]]))
        })
        assays <- c(list(own, nbr), out)
        assays <- assays[seq_len(min(K + 2, length(assays)))]
        message('BANKSY matrix with own.expr, ', 
                paste0('F', seq(0, length(assays)-2), collapse = ', '))
        lambdas <- getLambdas(lambda, n_harmonics = length(assays)-1)
        assays <- Map(function(lam, mat) lam * mat, lambdas, assays)
        joint <- do.call(rbind, assays)
        
        locs <- do.call(rbind, bank@cell.locs)
        rownames(locs) <- colnames(joint)
        
    } else {
        # Single dataset case
        # Consolidate own, F0, and higher harmonics
        assays <- c(list(bank@own.expr, bank@nbr.expr), bank@harmonics)
        assays <- assays[seq_len(min(K + 2, length(assays)))]
        message('BANKSY matrix with own.expr, ', 
                paste0('F', seq(0, length(assays)-2), collapse = ', '))
        # Compute lambdas
        lambdas <- getLambdas(lambda, n_harmonics = length(assays)-1)
        # Multiple by lambdas
        assays <- Map(function(lam, mat) lam * mat, lambdas, assays)
        # Row concatenate
        joint <- do.call(rbind, assays)
        locs <- bank@cell.locs
    }
    
    return(list(expr = joint, locs = locs))
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

#' @importFrom data.table `:=`
computeNeighbors <- function(locs,
                             spatial_mode = 'kNN_r', k_geom = 10, n = 2,
                             sigma = 1.5, alpha = 0.05, k_spatial = 100, 
                             dimensions = 'all', verbose = FALSE){
    from <- to <- phi <- NULL
    out <- getSpatialDims(locs, dimensions, alpha)
    locs <- out[[1]]
    kernelRadius <- out[[2]]
    cellID <- out[[3]]
    
    if (verbose) message('Computing neighbors...')
    if (spatial_mode == 'rNN_gauss'){
        knnDF <- withRNNgauss(locs = locs, sigma = sigma, kspatial = k_spatial,
                              kernelRadius = kernelRadius, verbose = verbose)
        
    } else if (spatial_mode == 'kNN_rank' ){
        knnDF <- withKNNrank(locs = locs, k_geom = k_geom, verbose = verbose)
        
    } else if (spatial_mode =='kNN_r'){
        knnDF <- withKNNr(locs = locs, k_geom = k_geom, verbose = verbose)
        
    } else if (spatial_mode == 'kNN_rn') {
        knnDF <- withKNNrn(locs, k_geom = k_geom, n = n, verbose = verbose)
        
    } else if (spatial_mode == 'kNN_unif') {
        knnDF <- withKNNunif(locs = locs, k_geom = k_geom, verbose = verbose)
        
    } else if (spatial_mode == 'kNN_median') {
        knnDF <- withKNNmedian(locs = locs, k_geom = k_geom, verbose = verbose)
        
    } else {
        stop('Invalid spatial_mode. 
             One of rNN_gauss, kNN_rank, kNN_r, kNN_unif, kNN_median')
    }
    knnDF[, phi := getPhi(locs, from, to), by=from][]
    if (verbose) message('Done')
    return(knnDF)
}


computeHarmonics <- function(gcm, knn_df, K, center){
    from <- to <- weight <- phi <- NULL 
    j = sqrt(as.complex(-1))
    
    if (any(dim(gcm)==0)) return(NULL)
    
    suffix <- ifelse(K == 0, '.nbr', paste0('.k', K))
    
    message('Computing harmonic k = ', K)
    if (center) {
        aggr <- knn_df[, abs(
            fscale(gcm[, to, drop=FALSE]) %*% (weight * exp(j*K*phi))
        ), by = from]
    } else {
        aggr <- knn_df[, abs(
            gcm[, to, drop=FALSE] %*% (weight * exp(j*K*phi))
        ), by = from]
    }
    ncm <- matrix(aggr$V1, nrow = nrow(gcm), ncol = ncol(gcm))
    rownames(ncm) <- paste0(rownames(gcm), suffix)
    colnames(ncm) <- colnames(gcm)
    message('Done')
    
    return(ncm)
}


getPhi <- function(locs, from, to) {
    out = sweep(locs[to,], 2, locs[from,], '-')
    phi = atan2(out[,2], out[,1]) 
    phi + as.integer(phi < 0) * 2*pi
}

mmult <- function(mat, vec) {
    t(t(mat)*vec)
}

#' @importFrom dbscan kNN
#' @importFrom data.table data.table setDT  setnames `:=` rbindlist
#' @importFrom stats dnorm
withRNNgauss <- function(locs, sigma, kspatial, kernelRadius, verbose) {
    
    if (verbose) message('Computing Banksy matrix')
    if (verbose) message('Spatial mode is rNN gaussian')
    if (verbose) message('Parameters: sigma = ', sigma, ', kspatial = ', kspatial)
    
    tryCatch({
        knn <- dbscan::kNN(x = locs, k = kspatial)
    },
    error=function(cond) {
        message("Not enough neighbours at kspatial = ", kspatial, " level.")
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
    if (verbose) message('Parameters: k_geom = ', k_geom)
    
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
        message("Not enough neighbours at kspatial = ", k_geom, " level.")
        message(cond)
    })
    
    return(knnDF)
}

#' @importFrom dbscan kNN
#' @importFrom data.table data.table setnames `:=`
withKNNr <- function(locs, k_geom, verbose) {
    
    if (verbose) message('Computing Banksy matrix')
    if (verbose) message('Spatial mode is kNNr')
    if (verbose) message('Parameters: k_geom = ', k_geom)
    
    tryCatch({
        knn <- dbscan::kNN(x = locs, k = k_geom)
    },
    error=function(cond) {
        message("Not enough neighbours at k_geom = ", k_geom, " level.")
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
withKNNrn <- function(locs, k_geom, n, verbose) {
    
    if (verbose) message('Computing Banksy matrix')
    if (verbose) message('Spatial mode is kNNrn')
    if (verbose) message('Parameters: k_geom = ', k_geom, ', n = ', n)
    
    tryCatch({
        knn <- dbscan::kNN(x = locs, k = k_geom)
    },
    error=function(cond) {
        message("Not enough neighbours at k_geom = ", k_geom, " level.")
        message(cond)
    })
    
    norm.weight <- weight <- from <- NULL
    knnDF <- data.table(from = rep(seq_len(nrow(knn$id)), k_geom),
                        to = as.vector(knn$id),
                        weight = 1/(as.vector(knn$dist)^n),
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
    if (verbose) message('Parameters: k_geom = ', k_geom)
    
    tryCatch({
        knn <- dbscan::kNN(x = locs, k = k_geom)
    },
    error=function(cond) {
        message("Not enough neighbours at kspatial = ", k_geom, " level.")
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

#' @importFrom dbscan kNN
#' @importFrom data.table data.table setnames `:=`
withKNNmedian <- function(locs, k_geom, verbose) {
    
    if (verbose) message('Computing Banksy matrix')
    if (verbose) message('Spatial mode is kNN median')
    if (verbose) message('Parameters: k_geom = ', k_geom)
    
    tryCatch({
        knn <- dbscan::kNN(x = locs, k = k_geom)
    },
    error=function(cond) {
        message("Not enough neighbours at kspatial = ", k_geom, " level.")
        message(cond)
    })
    
    norm.weight <- weight <- from <- distance <- NULL
    knnDF <- data.table(from = rep(seq_len(nrow(knn$id)), k_geom),
                        to = as.vector(knn$id),
                        distance = as.vector(knn$dist))
    knnDF[, weight := exp(-distance^2 / median(distance)^2), by=from]
    knnDF[, norm.weight := weight / sum(weight), by = from]
    knnDF <- knnDF[,-c(3,4),with=FALSE]
    setnames(knnDF, 'norm.weight', 'weight')
    
    return(knnDF)
}
