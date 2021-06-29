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

#' @importFrom data.table as.data.table setDT `:=`
#' @importFrom dbscan kNN
#' @importFrom stats median dnorm
#' @importFrom tictoc tic toc
compute.banksyMatrices <- function(gcm, locs,
                                   sigma=3, alpha=0.05,
                                   kspatial=1000, k_geom = 10,
                                   dimensions = 'all',
                                   spatialMode = c('kNN_r', 'kNearestNeighbours',
                                                   'kNN_unif','kNN_rank', 'rNN_gauss'),
                                   multiple_datasets = FALSE){

    # Global binding
    norm_wt <- weight <- from <- to <- NULL
    # spatial_locations = as.data.table(cbind(locs , cell_ID = rownames(locs)))
    spatial_locations = as.data.table(locs )
    cell_ID_vec <- c(1:nrow(spatial_locations))
    names(cell_ID_vec) <- locs$cell_ID
    # cell_ID_vec = locs$cell_ID #!TODO what was this?
    spatial_locations =
      spatial_locations[, grepl("sdim",
                                colnames(spatial_locations)),
                        with = FALSE]
    if (dimensions != "all") {
      spatial_locations = spatial_locations[, dimensions]
    }
    numDims = ncol(spatial_locations)
    nx = length(unique(spatial_locations$sdimx)) # !!!! TODO FIX THIS.
    ny = length(unique(spatial_locations$sdimy))
    if (numDims==3){
        nz = length(unique(spatial_locations$sdimz))
        if (nz > 1){
            dataDimensions = 3
            kernel_radius = sqrt(-3*log(alpha)) # !TODO figure out the math here.
        } else {
            dataDimensions = 2
            kernel_radius = sqrt(-2*log(alpha))
        }
    } else if (numDims==2){
          dataDimensions = 2
          kernel_radius = sqrt(-2*log(alpha))
    } else {
          stop('The location data must have 2 or 3 dimensions.')
    }
    spatial_locations <- as.matrix(spatial_locations)
    rownames(spatial_locations) <- locs$cell_ID

    # different spatial neighbourhood weighting functions
    if (spatialMode == 'rNN_gauss'){
          message(paste0('Computing BANKSY Matrix, Sigma = ', sigma, ',
                  Alpha = ', alpha, ', kspatial = ', kspatial))
          message(paste0('Spatial mode is ' , spatialMode))
      tryCatch({
          knn_spatial <- dbscan::kNN(x = spatial_locations, k = kspatial)
      },
      error=function(cond) {
            message(paste("Not enough neighbours at kspatial = ", kspatial, " level."))
            message("Here's the original error message:")
            message(cond)
      })
      median_dist_to_NN = median( as.vector(knn_spatial$dist[,1]))
      SIGMA = sigma
      knn_sptial.norm = data.frame(from = rep(1:nrow(knn_spatial$id),kspatial),
                                   to = as.vector(knn_spatial$id),
                                   weight = dnorm(as.vector(knn_spatial$dist),
                                                  mean = 0,
                                                  sd=median_dist_to_NN*SIGMA),
                                   distance = as.vector(knn_spatial$dist))
      knn_sptial.norm2 = knn_sptial.norm
      knn_sptial.norm2 =
          knn_sptial.norm2[as.vector(knn_sptial.norm2[,'distance']<SIGMA*kernel_radius*median_dist_to_NN),]
      setDT(knn_sptial.norm2)[, norm_wt := weight/sum(weight), by = from]
      knn_sptial.norm2 = knn_sptial.norm2[,c(1,2,5,4)]
      colnames(knn_sptial.norm2)[3] = 'weight'
    } else if (spatialMode == 'kNearestNeighbours' || spatialMode == 'kNN_rank' ){
          print(paste0('spatial mode is kNN_rank, k_geom = ', k_geom))
          tryCatch({
              knn_spatial <- dbscan::kNN(x = spatial_locations,
                                   k = k_geom)
              unnormWt = exp(-seq(1, k_geom, 1)^2 / (2*(k_geom/1.5)^2))
              normWt = unnormWt/sum(unnormWt)

          weightMatrix = t(matrix(normWt, nrow = k_geom, ncol = nrow(knn_spatial$id)))
          knn_sptial.norm = data.frame(from = rep(1:nrow(knn_spatial$id),k_geom),
                                     to = as.vector(knn_spatial$id),
                                     weight = as.vector(weightMatrix),
                                     distance = as.vector(knn_spatial$dist))
          knn_sptial.norm2 = setDT(knn_sptial.norm)
          },
          error=function(cond) {
              message(paste("Not enough neighbours at kspatial = ", k_geom, " level."))
              message("Here's the original error message:")
              message(cond)
          })
    } else if (spatialMode=='kNN_r'){
          message(paste0('Spatial mode is ' , spatialMode, ', k_geom = ', k_geom))
          tryCatch({
              knn_spatial <- dbscan::kNN(x = spatial_locations,
                                   k = k_geom)
          },
          error=function(cond) {
              message(paste("Not enough neighbours at kspatial = ", k_geom, " level."))
              message("Here's the original error message:")
              message(cond)
          })
      knn_sptial.norm = data.frame(from = rep(1:nrow(knn_spatial$id),k_geom),
                                   to = as.vector(knn_spatial$id),
                                   weight = 1/as.vector(knn_spatial$dist),
                                   distance = as.vector(knn_spatial$dist))
      knn_sptial.norm2 = knn_sptial.norm
      # knn_sptial.norm2 =
      #   knn_sptial.norm2[as.vector(knn_sptial.norm2[,'distance']<SIGMA*kernel_radius*median_dist_to_NN),]
      setDT(knn_sptial.norm2)[, norm_wt := weight/sum(weight), by = from]
      knn_sptial.norm2 = knn_sptial.norm2[,c(1,2,5,4)]
      colnames(knn_sptial.norm2)[3] = 'weight'

    } else {
          # just kNN with uniform weights.
          message(paste0('Spatial mode is kNN_unif, k_geom = ', k_geom))
          tryCatch({
              knn_spatial <- dbscan::kNN(x = spatial_locations,
                                       k = k_geom)
          },
          error=function(cond) {
              message(paste("Not enough neighbours at kspatial = ", kspatial, " level."))
              message("Here's the original error message:")
              message(cond)
          })

          knn_sptial.norm = data.frame(from = rep(1:nrow(knn_spatial$id),k_geom),
                                       to = as.vector(knn_spatial$id),
                                       weight = 1,
                                       distance = as.vector(knn_spatial$dist))
          knn_sptial.norm2 = knn_sptial.norm
          # knn_sptial.norm2 =
          #   knn_sptial.norm2[as.vector(knn_sptial.norm2[,'distance']<SIGMA*kernel_radius*median_dist_to_NN),]
          setDT(knn_sptial.norm2)[, norm_wt := weight/sum(weight), by = from]
          knn_sptial.norm2 = knn_sptial.norm2[,c(1,2,5,4)]
          colnames(knn_sptial.norm2)[3] = 'weight'
    }


    fromcellids = unique(knn_sptial.norm2$from)
    fromcells = names(cell_ID_vec[fromcellids])
    num_from_cell = length(fromcells)
    cell_expr_matrix = gcm
    nbr_expr_matrix <- cell_expr_matrix
    nbr_expr_matrix[] <- 0L
    # can we speed up this for loop? yes. I think we can do the
    # data table method of the old method,
    # 100 genes at a time. Nice idea, since the gene computations are independent.
    # This will not take up too much memory, and will be **much** faster.
    tic('Banksy matrix')
    for (fromcellid in fromcellids){
        cellNames = names(cell_ID_vec)
        fromcell = cellNames[fromcellid]
        tocells = names(cell_ID_vec[knn_sptial.norm2[from==fromcellid,to]])
        weights = knn_sptial.norm2[from==fromcellid,weight]
        weighting_vec = matrix(weights, nrow=length(weights), ncol = 1)
        nbr_expr_matrix[,fromcell] = cell_expr_matrix[, tocells] %*% weighting_vec
    }
    neighbourExpressionMatrix = nbr_expr_matrix
    rownames(nbr_expr_matrix) <- paste0(rownames(cell_expr_matrix), '.nbr')
    banksyMatrices = list(cellMatrix = cell_expr_matrix, nbrMatrix = nbr_expr_matrix)

    toc()
    return(nbr_expr_matrix)
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

