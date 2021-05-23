# compute.banksyMatrices
# zScaleDatasets
# splitByDataset
# addClippedZCustomExpr

#' @importFrom data.table as.data.table setDT `:=`
#' @importFrom dbscan kNN
#' @importFrom Matrix t
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
                        with = F]
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
    return(banksyMatrices)
}

#' @importFrom Matrix t
zScaleDatasets <- function(dataArray, multiple_datasets = FALSE,
                           zScaleDatasetsIndividually = TRUE,
                           cells_split_by_dataset=NULL){
    # dataArray is a matrix of genes (rows) by cells (columns). The cells are labelled as
    # <datasetID>__cell_1234
    # we use this to subset the matrices for the individual zScaling.
    # Actually, in this initial version, we do not do this, and instead use the
    # cells_split_by_dataset array to subset out the cells.
    zScaledDataArray=NULL

    if (multiple_datasets){
        # print('1')
        # str(zScaledDataArray)
        # str(dataArray)
        if (zScaleDatasetsIndividually){
            # print('2')
            # str(zScaledDataArray)
            # str(dataArray)
            splittedData = splitByDataset(dataArray=dataArray,locsArray = NULL,
                                      cells_split_by_dataset=cells_split_by_dataset)

            datasetNames = names(splittedData$expressionList)
            # print('2.5')
            # str(splittedData)
            # str(datasetNames)
            numDatasets = length(datasetNames)
            for (dataset_name in datasetNames){
                current_dataset = splittedData$expressionList[[dataset_name]] #dataArray[,current_cell_IDs]
                current_dataset = t(scale(t(current_dataset)))
                zScaledDataArray = cbind(zScaledDataArray, current_dataset)
                # print('3')
                # str(zScaledDataArray)
            }
      } else { # zScale all the datasets as one.
            # Not a good idea if the datasets represent different batches, like in the brain organoid data.
            # Possibly a good idea if they simply represent non-adjacent FOVs, like in the MIBI-TOF dataset.
            zScaledDataArray = t(scale(t(dataArray)))
      }
    } else { #if there is only one dataset, then we just zscale it.
          zScaledDataArray = t(scale(t(dataArray)))
    }
    # print('4')
    # str(zScaledDataArray)
    # str(dataArray)
    return(zScaledDataArray)
}

splitByDataset <- function(dataArray=NULL, locsArray=NULL, cells_split_by_dataset=NULL){
    datasetNames = names(cells_split_by_dataset)
    numDatasets = length(datasetNames)
    dataArrayList = vector(mode = "list", length = 0) # very important to not initialize lists with NULL. does not work.
    locsArrayList = vector(mode = "list", length = 0)
    for (dataset_name in datasetNames){
        current_cell_IDs = cells_split_by_dataset[[dataset_name]]
        if (!is.null(dataArray)){
            dataArrayList[[dataset_name]] = dataArray[,current_cell_IDs]
        }
        if (!is.null(locsArray)){
            locsArrayList[[dataset_name]] = locsArray[current_cell_IDs,]
        }
    }
    outputList = list(expressionList = NULL, locsList = NULL)
    if (!is.null(dataArray)){
        outputList$expressionList = dataArrayList
    }
    if (!is.null(locsArray)){
        outputList$locsList = locsArrayList
    }
    return(outputList)
}

addClippedZCustomExpr <- function(gobject,
                                  pairConnectedOutput,
                                  sweepresults,
                                  sweepNum = 1,
                                  upper_clip = 2,
                                  lower_clip = -2){
  i = sweepNum
  gridSize = dim(sweepresults$clusteringNames)
  currName = names(pairConnectedOutput$clusterings)[i]
  current_row_in_grid = floor((i-1)/gridSize[2]+1) # i is counting row wise.
  current_col_in_grid = i - gridSize[2]*(floor((i-1)/gridSize[2]))
  name_from_grid = sweepresults$clusteringNames[[current_row_in_grid, current_col_in_grid]]
  name_from_clusteringDF = currName
  if (name_from_grid != name_from_clusteringDF){
    stop("The clustering being referred to in the grid and the clusteringDF should be the same. Probably a row wise vs col wise indexing issue.")
  }

  own_expr_s = sweepresults$ownExpScaledMatrix
  nbr_expr_s = sweepresults$nbrExpScaledMatrix

  currlam = sweepresults$banksyParamsList[[current_row_in_grid, current_col_in_grid]]$lambda
  gobject@norm_scaled_expr = rbind(sqrt(1-currlam)*own_expr_s, sqrt(currlam)*nbr_expr_s)


  currGexp = gobject@norm_scaled_expr
  currGexp[currGexp > upper_clip] = upper_clip
  currGexp[currGexp < lower_clip] = lower_clip
  gobject@custom_expr = currGexp
  return(gobject)
}
