# RunBanksy
# SweepBanksy
# SubsetBanksy
# ComputeBanksy
# ClusterBanksy
# ConnectClusters

#' Subset a BanksyObject
#' @param x a BanksyObject
#' @param cells cells to filter by
#' @param dimx dimx to filter by - must be valid column
#' @param dimy dimy to filter by - must be valid column
#' @param dimz dimz to filter by - must be valid column
#' @param features genes to filter by
#' @param metadata metadata to filter by - must be valid colulmn
#'
#' @importFrom rlang enquo quo_get_expr
#'
#' @return Subset BanksyObject
#'
#' @export
SubsetBanksy <- function(
  x,
  cells=TRUE,
  dimx=TRUE,
  dimy=TRUE,
  dimz=TRUE,
  features=NULL,
  metadata=TRUE) {

  nfeaturesBef <- nrow(x@own.expr)
  ncellsBef <- ncol(x@own.expr)

  ## Filter cells by dimension and cells
  dimx <- rlang::enquo(dimx)
  x@cell.locs <- subset(data.frame(x@cell.locs), subset = eval(rlang::quo_get_expr(dimx)))

  dimy <- rlang::enquo(dimy)
  x@cell.locs <- subset(data.frame(x@cell.locs), subset = eval(rlang::quo_get_expr(dimy)))

  dimz <- rlang::enquo(dimz)
  x@cell.locs <- subset(data.frame(x@cell.locs), subset = eval(rlang::quo_get_expr(dimz)))

  cells <- rlang::enquo(cells)
  x@own.expr <- subset(data.frame(x@own.expr), select = eval(rlang::quo_get_expr(cells)))

  ## Filter metadata
  metadata <- rlang::enquo(metadata)
  x@meta.data <- subset(x@meta.data, subset = eval(rlang::quo_get_expr(metadata)))

  ## Filter features
  if(!(is.null(features))) {
    x@own.expr <- x@own.expr[rownames(x@own.expr)%in%features,,drop=FALSE]
    x@nbr.expr <- x@nbr.expr[rownames(x@nbr.expr)%in%paste0(features,'.nbr'),,drop=FALSE]
    x@own.norm.scaled.expr <- x@own.norm.scaled.expr[rownames(x@own.norm.scaled.expr)%in%features,,drop=FALSE]
    x@nbr.norm.scaled.expr <- x@nbr.norm.scaled.expr[rownames(x@nbr.norm.scaled.expr)%in%paste0(features,'.nbr'),,drop=FALSE]
    x@custom.expr <- x@custom.expr[rownames(x@custom.expr)%in%features,,drop=FALSE]
  }

  ## Consistent Cell filtering
  surviving_cells <- intersect(x@meta.data$cell_ID,
                               intersect(colnames(x@own.expr), rownames(x@cell.locs)))
  x@own.expr <- x@own.expr[,surviving_cells,drop=FALSE]
  x@nbr.expr <- x@nbr.expr[,surviving_cells,drop=FALSE]
  x@own.norm.scaled.expr <- x@own.norm.scaled.expr[,surviving_cells,drop=FALSE]
  x@nbr.norm.scaled.expr <- x@nbr.norm.scaled.expr[,surviving_cells,drop=FALSE]
  x@custom.expr <- x@custom.expr[,surviving_cells,drop=FALSE]
  x@cell.locs <- x@cell.locs[surviving_cells,,drop=FALSE]
  x@meta.data <- x@meta.data[x@meta.data$cell_ID %in% surviving_cells,,drop=FALSE]
  x@dim.reduction <- lapply(x@dim.reduction, function(dim.red) {
    dim.red <- dim.red[rownames(dim.red)%in%surviving_cells,,drop=FALSE]
    dim.red
  })

  ## Report metrics
  nfeaturesAft <- nrow(x@own.expr)
  ncellsAft <- ncol(x@own.expr)
  message('Before filtering: ', ncellsBef, ' cells and ', nfeaturesBef, ' genes.')
  message('After filtering: ', ncellsAft, ' cells and ', nfeaturesAft, ' genes.')
  message('Filtered ', ncellsBef - ncellsAft, ' cells and ', nfeaturesBef - nfeaturesAft, ' genes.')

  return(x)
}

#' Compute Banksy Matrices
#' @param bank BanksyObject
#' @param normalizeColumns to normalize columns
#' @param normalizeColumnsTo a scale factor
#' @param zScaleRows to zscale rows
#' @param zScaleBeforeAveraging to zscale before avg
#' @param zScaleOwnAfterAveraging to zscale own after avg
#' @param zScaleNbrAfterAveraging to zscale nbr after avg
#' @param sigma sigma
#' @param alpha alpha
#' @param kspatial kspatial
#' @param dimensions dims
#' @param spatialMode spatialMode
#' @param k_geom k_geom
#'
#' @importFrom Matrix t
#' @importFrom data.table data.table setnames
#'
#' @return Banksy Object
#'
#' @export
ComputeBanksy <- function(bank,
                                  ## Normalization
                                  normalizeColumns = TRUE,
                                  normalizeColumnsTo = 100,
                                  ## Scaling
                                  zScaleRows = TRUE,
                                  zScaleBeforeAveraging = FALSE,
                                  zScaleOwnAfterAveraging = TRUE,
                                  zScaleNbrAfterAveraging = TRUE,
                                  ## For computing nbr matrix
                                  sigma = 1.5,
                                  alpha = 0.05,
                                  kspatial = 1000,
                                  dimensions = 'all',
                                  spatialMode = 'kNN_r',
                                  k_geom = 10) {


  gcm <- .sparsify(bank@own.expr)

  message('Performing normalization...')
  if (normalizeColumns) {
    gcm <- t(t(gcm)/colSums(gcm))*normalizeColumnsTo
    bank@own.norm.scaled.expr <- gcm
  }

  if (zScaleRows & zScaleBeforeAveraging) {
    gcm <- t(scale(t(gcm)))
    bank@own.norm.scaled.expr <- gcm
  }

  locs <- data.table(bank@cell.locs, keep.rownames = TRUE)
  setnames(locs, 'rn', 'cell_ID')

  message('Computing Banksy matrices...')
  banksyMatrices = compute.banksyMatrices(gcm, locs,
                                          sigma=sigma,
                                          alpha=alpha,
                                          kspatial=kspatial,
                                          dimensions = dimensions,
                                          spatialMode = spatialMode ,
                                          k_geom = k_geom)

  bank@nbr.expr <- .framify(banksyMatrices$nbrMatrix)     ## called nbr_concat

  if (zScaleRows){
    if (zScaleOwnAfterAveraging) bank@own.norm.scaled.expr <- .framify(zScaleDatasets(banksyMatrices$cellMatrix))     ## called cell_expr_matrix_s
    if (zScaleNbrAfterAveraging) bank@nbr.norm.scaled.expr <- .framify(zScaleDatasets(banksyMatrices$nbrMatrix))       ## called nbr_expr_matrix_s
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
#' @param leiden_iters leiden iters
#' @param pypath pypath
#' @param use_docker use docker
#' @param verbose messages
#'
#' @importFrom Giotto createGiottoInstructions createGiottoObject runPCA runUMAP
#'   createNearestNetwork doLeidenCluster
#'
#' @return BanksyObject with cluster labels in meta.data
#'
#' @export
ClusterBanksy <- function(bank,
                          lambda = 0.5,
                          resolution = 1.2,
                          kneighbours = 30,
                          npcs = 50,
                          leiden_iters = -1,
                          pypath = '/usr/bin/python3',
                          use_docker = TRUE,
                          verbose = FALSE) {

  tic()
  max_iters <- prod(length(lambda), length(resolution), length(kneighbours))
  instrs = createGiottoInstructions(python_path = pypath,
                                    show_plot = FALSE,
                                    return_plot = TRUE,
                                    save_plot = TRUE,
                                    save_dir = '/data',
                                    plot_format = 'png',
                                    dpi = 300, height = 9, width = 9,
                                    is_docker = use_docker)
  iter = 1
  message(paste0('Iteration ', iter, ' out of ', max_iters) )

  for (lam in lambda) {

    joint <- rbind(sqrt(1-lam)*bank@own.norm.scaled.expr,
                   sqrt(lam)*bank@nbr.norm.scaled.expr)
    giotto <- createGiottoObject(raw_exprs = joint,
                                 norm_expr = joint,
                                 norm_scaled_expr = joint,
                                 spatial_locs = bank@cell.locs,
                                 instructions = instrs)
    if (verbose) message(paste0('Dimensionality reduction for Lambda=',lam))
    giotto <- runPCA(giotto,
                     expression_values = 'normalized',
                     scale_unit = FALSE,
                     ncp = npcs)
    giotto <- runUMAP(giotto,
                      expression_values = 'normalized',
                      dimensions_to_use = 1:npcs)
    bank@dim.reduction[[paste0('umap_', lam)]] <- giotto@dimension_reduction$cells$umap$umap$coordinates

    for (res in resolution) {
      for (k in kneighbours) {

        if (verbose) message(paste0('Create Nearest Neighbour Network and perform Leiden clustering
                                              for Lambda=',lam, ', Resolution=', res, ', K Neighbours=', k))
        giotto <- createNearestNetwork(giotto,
                                       dimensions_to_use = 1:npcs,
                                       k = k)
        giotto <- doLeidenCluster(giotto,
                                  resolution = res,
                                  n_iterations = leiden_iters,
                                  name = paste0('res', res, '_lam', lam, '_k', k))
        message(paste0('Finished clustering for Lambda=', lam, ', Resolution=', res, ', K Neighbours=', k))
        iter <- iter + 1
        if (iter <= max_iters) message(paste0('Iteration ', iter, ' out of ', max_iters) )

      }
    }

    bank@meta.data <- cbind(bank@meta.data, giotto@cell_metadata[,-1])

  }
  toc()
  return(bank)

}


#' Harmonise cluster labels among parameter runs
#'
#' @param bank Banksy Object
#' @param verbose verbose or not
#' @param optim optimise or not (KS-test)
#'
#' @importFrom Giotto getDistinctColors
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
  parent <- max(which(numClust == medClust))
  parentClusters <- allClust[[parent]]
  newClust[,parent] <- plyr::mapvalues(clust[,parent],
                                       from = parentClusters,
                                       to = seq_len(medClust),
                                       warn_missing = FALSE)
  parentDist <- as.numeric(table(newClust[,parent]))
  message(paste0('Mapping clusterings to ', clustNames[parent]))
  ## The rest will be children
  children <- setdiff(seq_len(ncol(clust)), parent)
  child<-children[3]
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
          if (verbose) message(paste0('Ambiguous mapping for cluster ', i, ' - using KS-test'))
          hits <- as.numeric(levels(map)[tab %in% topn])
          test1 <- childDist + .init(medClust, hits[1],
                                     tab[which(levels(map)==hits[1])])
          test2 <- childDist + .init(medClust, hits[2],
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
  toMap <- Giotto::getDistinctColors(maxClust)
  ## Define a color mapping
  for (i in seq_len(ncol(newClust))) {
    clustName <- names(newClust)[i]
    colName <- paste0('col_', clustName)
    plotCols <- plyr::mapvalues(newClust[,i],
                                from = fromMap,
                                to = toMap,
                                warn_missing = FALSE)
    newClust[,colName] <- plotCols
  }
  bank@meta.data <- cbind(cell_ID = bank@meta.data$cell_ID, newClust)
  return(bank)
}


#' This function takes the gene-cell matrix and cell locations, and runs banksy.
#'
#' @param gcm gene-cell matrix, rows are features, columns are objects
#' @param locs cell locations with object locations: x and y (and possibly z)
#' coordinates. x and y are necessary, and must have at least two different
#' locations. If data is in a single plane, it must be in a single z plane.
#' The value of this z plane variable may be provided, but it is not necessary.
#' There needs to be a column labeled cell_ID, which contains strings that match
#' the column names of gcm, in the same order as the columns of gcm. locs should
#' be a data.table with rows and columns as follows:
#' sdimx     sdimy  sdimz    cell_ID
#' @param sigma sigma
#' @param lambda lambda
#' @param alpha alpha
#' @param pcs pcs
#' @param k_expr k_expr for shared neighbourhood clustering
#' @param res res
#' @param kspatial kspatial
#' @param normalizeColumns norm
#' @param normalizeColumnsTo norm
#' @param minimumTotalCounts min
#' @param lieden_iters lie
#' @param zScaleRows z
#' @param zScaleBeforeAveraging z
#' @param zScaleOwnAfterAveraging z
#' @param zScaleNbrAfterAveraging z
#' @param zScaleDatasetsIndividually z
#' @param multiple_datasets multiple
#' @param cells_split_by_dataset no
#' @param dimensions dims
#' @param clustering_name name
#' @param k_geom k
#' @param banksyMatrices mat
#' @param spatialMode mode
#' @param compute_own_and_nbr_umap no
#' @param compute_3D_umap no
#' @param instructions instructions
#' @param ... more params
#'
#' @importFrom tictoc tic toc
#' @importFrom Giotto createGiottoObject normalizeGiotto runPCA runUMAP
#'   createNearestNetwork doLeidenCluster addStatistics
#' @importFrom Matrix t colSums
#'
#' @return List of clustering results
#'
#' @export
RunBanksy <- function(gcm, locs,
                      sigma = 3, lambda = 0.25, alpha = 0.05, pcs = 15,
                      k_expr = 30, res = 0.5, kspatial = 1000,
                      normalizeColumns = TRUE,
                      normalizeColumnsTo = 100,
                      minimumTotalCounts = 10,
                      lieden_iters = -1,
                      zScaleRows = TRUE,
                      zScaleBeforeAveraging = FALSE,
                      zScaleOwnAfterAveraging = TRUE,
                      zScaleNbrAfterAveraging = TRUE,
                      zScaleDatasetsIndividually = TRUE,
                      multiple_datasets = FALSE,
                      cells_split_by_dataset = NULL,
                      dimensions = 'all',
                      clustering_name='banksy', k_geom = 10,
                      banksyMatrices = NULL,
                      spatialMode = c('kNN_r', 'kNearestNeighbours',
                                     'kNN_unif','kNN_rank', 'rNN_gauss'),
                      compute_own_and_nbr_umap = FALSE,
                      compute_3D_umap = FALSE,
                      instructions, ...) {

    own_concat = NULL
    nbr_concat = NULL
    locs_concat = NULL
    LAM = lambda

    if (!(multiple_datasets)){
      if (is.null(banksyMatrices)){
        if (normalizeColumns==TRUE){
          gcm <- t(t(gcm)/colSums(gcm))*normalizeColumnsTo
        }
        if (zScaleRows==TRUE){
          if (zScaleBeforeAveraging==TRUE){
            gcm = t(scale(t(gcm)))
          }
        }
        banksyMatrices = compute.banksyMatrices(gcm, locs,
                                                sigma=sigma,
                                                alpha=alpha,
                                                kspatial=kspatial,
                                                dimensions = dimensions,
                                                spatialMode = spatialMode ,
                                                k_geom = k_geom)
        own_concat = banksyMatrices$cellMatrix
        nbr_concat = banksyMatrices$nbrMatrix
      } else {
        own_concat = banksyMatrices$cellMatrix
        nbr_concat = banksyMatrices$nbrMatrix
      }
      locs_concat = locs
    } else {
      # !! april 20 the zscale-before-averaging mode is not great right now.
      # this is because i think the right thing to do is to concat all the matrices,
      # zscale, then separate out and compute the neghbour matrices, then concat again.

      if (is.null(banksyMatrices)){
        gene_names = rownames(gcm[[1]])
        for (dataset_name in names(gcm)[-1]){
          gene_names = intersect(gene_names, rownames(gcm[[dataset_name]]))
        }
        for (dataset_name in names(gcm)){
          current_gcm = gcm[[dataset_name]][gene_names,]
          if (normalizeColumns==TRUE){
            cells_to_keep = colnames(current_gcm)[colSums(current_gcm)>=minimumTotalCounts]
            current_gcm1 = current_gcm[,cells_to_keep] #[,colSums(current_gcm)>=minimumTotalCounts]
            locs_to_keep = locs[[dataset_name]]$cell_ID %in% cells_to_keep
            current_locs = locs[[dataset_name]][locs_to_keep,]
            current_gcm = t(t(current_gcm1)/colSums(current_gcm1))*normalizeColumnsTo


            cells_split_by_dataset[[dataset_name]] = cells_split_by_dataset[[dataset_name]][cells_to_keep]
          } else {
            current_locs = locs[[dataset_name]]
          }
          locs_concat = rbind(locs_concat,current_locs)
          if (zScaleRows==TRUE){
            if (zScaleBeforeAveraging==TRUE){
              # current_gcm = zScaleDatasets(current_gcm,
              #                              multiple_datasets = multiple_datasets,
              #                              zScaleDatasetsIndividually = zScaleDatasetsIndividually,
              #                              cells_split_by_dataset=cells_split_by_dataset)
              current_gcm = t(scale(t(current_gcm))) # need to redo this part to allow for the joint zscaling.
            }
          }
          banksyMatrices = compute.banksyMatrices(current_gcm,
                                                  current_locs,
                                                  sigma=sigma,
                                                  alpha=alpha,
                                                  kspatial=kspatial,
                                                  dimensions = dimensions,
                                                  spatialMode = spatialMode,
                                                  k_geom = k_geom)
          own_concat = cbind(own_concat, banksyMatrices$cellMatrix)
          nbr_concat = cbind(nbr_concat, banksyMatrices$nbrMatrix)
        }
      }  else {
        own_concat = banksyMatrices$cellMatrix
        nbr_concat = banksyMatrices$nbrMatrix
        # print('test-B')
        # str(nbr_concat)
        locs_concat = locs[[1]]
        for (dataset_name in names(locs)[-1]){
          locs_concat = rbind(locs_concat,locs[[dataset_name]])
        }
      }
    }
    # print('test-C')
    # str(nbr_concat)
    if (zScaleRows==TRUE){
      if (zScaleOwnAfterAveraging){
        cell_expr_matrix_s = zScaleDatasets(own_concat,
                                            multiple_datasets = multiple_datasets,
                                            zScaleDatasetsIndividually = zScaleDatasetsIndividually,
                                            cells_split_by_dataset=cells_split_by_dataset)
      } else {
        cell_expr_matrix_s = own_concat
      }
      if (zScaleNbrAfterAveraging){
        # print('test-D')
        # str(nbr_concat)
        # str(cells_split_by_dataset)
        nbr_expr_matrix_s = zScaleDatasets(dataArray = nbr_concat,
                                           multiple_datasets = multiple_datasets,
                                           zScaleDatasetsIndividually = zScaleDatasetsIndividually,
                                           cells_split_by_dataset=cells_split_by_dataset)
        # print('A')
        # str(nbr_expr_matrix_s)
      } else {
        nbr_expr_matrix_s = nbr_concat
      }
    } else {
      cell_expr_matrix_s = own_concat
      nbr_expr_matrix_s = nbr_concat
    }

    ######################################################################
    ##  Create the own and neighbour matrix mixed with parameter lambda ##
    ######################################################################
    joint_matrix_raw = rbind(sqrt(1-LAM)*own_concat, sqrt(LAM)*nbr_concat)
    bo_jointly <- createGiottoObject(raw_exprs = joint_matrix_raw,
                                     spatial_locs =locs_concat,  #[,c('sdimx', 'sdimy')]
                                     instructions = instructions)
    bo_jointly <- normalizeGiotto(gobject = bo_jointly,
                                  library_size_norm = FALSE,log_norm = FALSE,
                                  scale_cells=F, scale_genes = F)
    bo_jointly@norm_scaled_expr <- rbind(sqrt(1-LAM)*cell_expr_matrix_s, sqrt(LAM)*nbr_expr_matrix_s)
    bo_jointly <- addStatistics(gobject = bo_jointly)
    exprs = "scaled"
    nGenes = nrow(joint_matrix_raw)/2
    own_genes = rownames(joint_matrix_raw)[1:nGenes]
    nbr_genes = rownames(joint_matrix_raw)[(nGenes+1):(2*nGenes)]
    ncp = 50 # !! april16 reduced from the default of 200 to 50. check that it does not mess anything up
    if (ncp<(pcs+10)){
      ncp = pcs+12
    }

    #############################
    ## Run PCA on joint matrix ##
    #############################
    message(paste0('Starting PCA for ', clustering_name, '.'))
    tic(paste0('PCA for ', clustering_name, '.'))
    bo_jointly <- runPCA(gobject = bo_jointly,
                         expression_values = exprs,
                         name = 'pca_joint_genes',
                         genes_to_use = NULL,
                         scale_unit = FALSE, # this is required for J to run irlba on local
                         ncp = ncp)
    toc()

    ########################################
    ## Run UMAP on the own and nbr matrix ##
    ########################################
    if (compute_own_and_nbr_umap == TRUE){
      bo_jointly <- runPCA(gobject = bo_jointly,
                           expression_values = exprs,
                           genes_to_use = own_genes,
                           name = 'pca_own_genes',
                           scale_unit = FALSE,
                           ncp = ncp)
      bo_jointly <- runPCA(gobject = bo_jointly,
                           expression_values = exprs,
                           genes_to_use = nbr_genes,
                           name = 'pca_nbr_genes',
                           scale_unit = FALSE,
                           ncp = ncp)
    }

    ##############################
    ## Run UMAP on joint matrix ##
    ##############################
    tic(paste0('UMAP for ', clustering_name, '.'))
    bo_jointly <- runUMAP(bo_jointly, dim_reduction_to_use = "pca",
                          dim_reduction_name = 'pca_joint_genes',
                          expression_values = exprs,
                          dimensions_to_use = 1:pcs,
                          name = paste0('joint_2D_', clustering_name))


    if(compute_own_and_nbr_umap == TRUE){
      bo_jointly <- runUMAP(bo_jointly, dim_reduction_to_use = "pca",
                            dim_reduction_name = 'pca_own_genes',
                            expression_values = exprs,
                            dimensions_to_use = 1:pcs,
                            name = paste0('own_2D_', clustering_name))
      bo_jointly <- runUMAP(bo_jointly, dim_reduction_to_use = "pca",
                            dim_reduction_name = 'pca_nbr_genes',
                            expression_values = exprs,
                            dimensions_to_use = 1:pcs,
                            name = paste0('nbr_2D_', clustering_name))
    }

    #################
    ## Run 3D UMAP ##
    #################
    if (compute_3D_umap==TRUE){
      # the giotto object will get large if 3D own and neighbour UMAPS are required,
      # and if the sweep is large. So if visualizing the own and neigbour umaps, and
      # especially in 3D, keep the sweep small.
      bo_jointly <- runUMAP(bo_jointly,dim_reduction_to_use = "pca",
                            dim_reduction_name = 'pca_joint_genes',
                            expression_values = exprs,
                            dimensions_to_use = 1:pcs,n_components = 3,
                            name = paste0('joint_3D_', clustering_name))
      if( compute_own_and_nbr_umap == TRUE ){
        bo_jointly <- runUMAP(bo_jointly, dim_reduction_to_use = "pca",
                              dim_reduction_name = 'pca_own_genes',
                              expression_values = exprs,
                              dimensions_to_use = 1:pcs,n_components = 3,
                              name = paste0('own_3D_', clustering_name))

        bo_jointly <- runUMAP(bo_jointly, dim_reduction_to_use = "pca",
                              dim_reduction_name = 'pca_nbr_genes',
                              expression_values = exprs,
                              dimensions_to_use = 1:pcs,n_components = 3,
                              name = paste0('nbr_3D_', clustering_name))
      }
    }
    toc()


    ######################################################
    ## Run shared Nearest Neighbourhood on joint matrix ##
    ######################################################
    tic(paste0('sNN for ', clustering_name, '.'))
    bo_jointly <- createNearestNetwork(gobject = bo_jointly,
                                       dim_reduction_to_use = "pca",
                                       dim_reduction_name = "pca_joint_genes",
                                       dimensions_to_use = 1:pcs,
                                       k = k_expr)
    toc()

    #######################
    ## Leiden clustering ##
    #######################
    tic(paste0('Clustering for ', clustering_name, '.'))
    bo_jointly <- doLeidenCluster(gobject = bo_jointly, resolution = res,
                                  n_iterations = lieden_iters, name = clustering_name)
    toc()

    ###################
    ## Banksy params ##
    ###################
    banksy.params = list(sigma = sigma, lambda = lambda,
                         alpha = alpha,
                         pcs = pcs,
                         k_expr = k_expr,
                         res = res,
                         kspatial = kspatial,
                         normalizeColumns = normalizeColumns ,
                         normalizeColumnsTo = normalizeColumnsTo,
                         zScaleRows = zScaleRows ,
                         lieden_iters = lieden_iters,
                         zScaleBeforeAveraging = zScaleBeforeAveraging ,
                         zScaleOwnAfterAveraging = zScaleOwnAfterAveraging,
                         zScaleNbrAfterAveraging = zScaleNbrAfterAveraging,
                         multiple_datasets = FALSE,
                         zScaleDatasetsIndividually = TRUE,
                         dimensions = dimensions,
                         clustering_name=clustering_name,
                         k_geom = k_geom,
                         spatialMode = spatialMode,
                         multiple_datasets = multiple_datasets,
                         compute_own_and_nbr_umap = compute_own_and_nbr_umap,
                         compute_3D_umap = compute_3D_umap)

    ###################################
    ## Conslidate clustering results ##
    ###################################
    clustering_results = list(clustering = bo_jointly@cell_metadata,
                              umap = bo_jointly@dimension_reduction$cells$umap, # this actually stores all the UMAPS.
                              gobject = bo_jointly,
                              parameterSet = bo_jointly@parameters,
                              ownExpScaledMatrix = cell_expr_matrix_s,
                              nbrExpScaledMatrix = nbr_expr_matrix_s,
                              zScaleRows=zScaleRows,
                              zScaleBeforeAveraging=zScaleBeforeAveraging,
                              banksy.params = banksy.params)

    return(clustering_results)
}


#' Sweep gene-cell matrix
#'
#' @param gcm gene-cell matrix, rows are features, columns are objects
#' @param locs cell locations with object locations: x and y (and possibly z)
#' coordinates. x and y are necessary, and must have at least two different
#' locations. If data is in a single plane, it must be in a single z plane.
#' The value of this z plane variable may be provided, but it is not necessary.
#' There needs to be a column labeled cell_ID, which contains strings that match
#' the column names of gcm, in the same order as the columns of gcm. locs should
#' be a data.table with rows and columns as follows:
#' sdimx     sdimy  sdimz    cell_ID
#' @param rowparam rowparam
#' @param rowvalues rowvals
#' @param colparam colparam
#' @param colvalues colvals
#' @param sigma sigma
#' @param lambda lambda
#' @param alpha alpha
#' @param pcs pcs
#' @param k_expr k_expr
#' @param res res
#' @param kspatial kspatial
#' @param k_geom k_geom
#' @param normalizeColumns norm
#' @param normalizeColumnsTo norm
#' @param minimumTotalCounts min
#' @param lieden_iters niters
#' @param zScaleRows z
#' @param zScaleBeforeAveraging z
#' @param zScaleOwnAfterAveraging z
#' @param zScaleNbrAfterAveraging z
#' @param zScaleDatasetsIndividually z
#' @param multiple_datasets m
#' @param cells_split_by_dataset cells
#' @param dimensions dims
#' @param spatialMode spatial
#' @param banksyMatrices bank
#' @param compute_own_and_nbr_umap no
#' @param compute_3D_umap no
#' @param instructions instrs
#'
#' @importFrom Matrix t
#' @importFrom data.table copy
#' @importFrom pracma meshgrid
#'
#' @return Sweep List
#'
#' @export
SweepBanksy <- function(gcm, locs,
                        rowparam = 'lam', rowvalues = c(0, 0.25, 0.5),  #, 0.4, 0.5
                        colparam = 'res' , colvalues = c(0.5, 1), #, 70
                        sigma = 3, lambda = 0.4, alpha = 0.05,
                        pcs = 15, k_expr = 30, res = 0.5,
                        kspatial=1000, k_geom = 10,
                        normalizeColumns = TRUE, normalizeColumnsTo = 100,
                        minimumTotalCounts = 0,
                        lieden_iters = -1,
                        zScaleRows = TRUE,
                        zScaleBeforeAveraging = FALSE,
                        zScaleOwnAfterAveraging = TRUE,
                        zScaleNbrAfterAveraging = TRUE,
                        zScaleDatasetsIndividually = TRUE,
                        multiple_datasets = FALSE,
                        cells_split_by_dataset = NULL,
                        dimensions = 'all',
                        spatialMode = c('kNN_r', 'kNearestNeighbours',
                                        'kNN_unif', 'kNN_rank', 'rNN_gauss'),
                        banksyMatrices = NULL,
                        compute_own_and_nbr_umap = FALSE,
                        compute_3D_umap = FALSE,
                        instructions) {

    nr = length(rowvalues)
    nc = length(colvalues)
    numpoints = nr*nc
    lambdaM = matrix(rep(lambda, numpoints), nr, nc)
    pcsM = matrix(rep(pcs, numpoints), nr, nc)
    k_exprM = matrix(rep(k_expr, numpoints), nr, nc)
    resM = matrix(rep(res, numpoints), nr, nc)
    meshes = meshgrid(colvalues, rowvalues)
    rowparams = meshes$Y
    colparams = meshes$X
    switch(rowparam, lam={lambdaM=rowparams}, pcs={pcsM = rowparams},
           k_expr={k_exprM=rowparams}, res={resM=rowparams}
    )
    switch(colparam,
           lam={lambdaM=colparams}, pcs={pcsM = colparams},
           k_expr={k_exprM=colparams}, res={resM=colparams}
    )

    clusteringList = vector("list", nr*nc)
    umapList = vector("list", nr*nc)
    clusteringNames = vector("list", nr*nc)
    dim(clusteringList)<- c(nr, nc)
    dim(umapList)<- c(nr, nc)
    dim(clusteringNames)<- c(nr, nc)
    parameterSetList = vector("list", nr*nc)
    dim(parameterSetList)<- c(nr, nc)
    banksyParamsList = vector("list", nr*nc)
    dim(banksyParamsList)<- c(nr, nc)
    own_concat = NULL
    nbr_concat = NULL
    locs_concat = NULL
    gcm_raw = vector(mode = "list", length = 0)
    gcm_scaled = vector(mode = "list", length = 0)
    if (is.null(banksyMatrices)){
        #############################################
        if (!(multiple_datasets)){
            # all dem politicians be the same.
            if (normalizeColumns==TRUE){
                gcm = t(t(gcm)/colSums(gcm))*normalizeColumnsTo
            }
            if (zScaleRows==TRUE){
                if (zScaleBeforeAveraging==TRUE){
                    gcm = t(scale(t(gcm)))
                }
            }
            banksyMatrices = compute.banksyMatrices(gcm, locs, sigma=sigma, alpha=alpha, kspatial=kspatial,
                                                  dimensions = dimensions, spatialMode = spatialMode , k_geom = k_geom)
            own_concat = banksyMatrices$cellMatrix
            nbr_concat = banksyMatrices$nbrMatrix
        } else {
            gene_names = rownames(gcm[[1]])
            for (dataset_name in names(gcm)[-1]){
                gene_names = intersect(gene_names, rownames(gcm[[dataset_name]]))
            }
            for (dataset_name in names(gcm)){
                # get the current gene cell matrix.
                current_gcm = gcm[[dataset_name]][gene_names,]
                # print(dataset_name)
                if (normalizeColumns==TRUE){
                    # normalizeColumns = TRUE means that we are looking at counts data. We want to filter that
                    # by some minimum counts threshold.
                    # current_gcm1 = current_gcm[,colSums(current_gcm)>=minimumTotalCounts]
                    # gcm_raw[[dataset_name]] = current_gcm1
                    # current_locs = locs[[dataset_name]][colSums(current_gcm)>=minimumTotalCounts,]
                    # current_gcm = t(t(current_gcm1)/colSums(current_gcm1))*normalizeColumnsTo
                    #  ##!!! !! !! this block was not commented out. Interesting.

                    cells_to_keep = colnames(current_gcm)[colSums(current_gcm)>=minimumTotalCounts]
                    current_gcm1 = current_gcm[,cells_to_keep] #[,colSums(current_gcm)>=minimumTotalCounts]
                    locs_to_keep = locs[[dataset_name]]$cell_ID %in% cells_to_keep
                    current_locs = locs[[dataset_name]][locs_to_keep,]
                    current_gcm = t(t(current_gcm1)/colSums(current_gcm1))*normalizeColumnsTo
                    cells_split_by_dataset[[dataset_name]] = cells_split_by_dataset[[dataset_name]][cells_to_keep]
                    gcm_raw[[dataset_name]] = current_gcm
                } else {
                    current_locs = locs[[dataset_name]]
                    gcm_raw[[dataset_name]] = current_gcm
                }
                locs_concat = rbind(locs_concat,current_locs)
                locs[[dataset_name]] = current_locs
                # we will feed this into run.banksy.
                # we want this to be the "raw" counts, which is why we dont use current_gcm,
                # which is the normalized counts when normalizeColumns = TRUE.
                if (zScaleRows && zScaleBeforeAveraging && zScaleDatasetsIndividually){
                    current_gcm = t(scale(t(current_gcm)))
                }
                gcm_scaled[[dataset_name]] = current_gcm
            }
            if (zScaleRows && zScaleBeforeAveraging && !(zScaleDatasetsIndividually)){
                gcm_concat = NULL
                for (dataset_name in names(gcm)){
                    gcm_concat = cbind(gcm_concat, gcm_scaled[[dataset_name]])
                }
                gcm_concat = t(scale(t(gcm_concat)))
                splittedData = splitByDataset(dataArray=gcm_concat,
                                              cells_split_by_dataset=cells_split_by_dataset)
                gcm_scaled = splittedData$expressionList
            }
            for (dataset_name in names(gcm)){
                banksyMatrices = compute.banksyMatrices(gcm_scaled[[dataset_name]],
                                                    locs[[dataset_name]],
                                                    sigma=sigma,
                                                    alpha=alpha,
                                                    kspatial=kspatial,
                                                    dimensions = dimensions,
                                                    spatialMode = spatialMode,
                                                    k_geom = k_geom)
                # print(dim(banksyMatrices$cellMatrix))
                # print(dim(banksyMatrices$nbrMatrix))
                # print(dim(gcm_scaled[[dataset_name]]))
                # print(dim(locs[[dataset_name]]))
                own_concat = cbind(own_concat, banksyMatrices$cellMatrix)
                nbr_concat = cbind(nbr_concat, banksyMatrices$nbrMatrix)
                # print(dim(own_concat))
                # print(dim(nbr_concat))
            }
        }
        banksyMatrices$cellMatrix = own_concat
        banksyMatrices$nbrMatrix = nbr_concat
        ##############################################

        # banksyMatrices <- compute.banksyMatrices(gcm, locs, sigma=sigma,
        #                                          alpha=alpha, kspatial=kspatial,
        #                                          dimensions = dimensions,
        #                                          spatialMode = spatialMode, k_geom = k_geom )
    }

    for (j in 1:nc){
        for (i in 1:nr){
            lambda = lambdaM[i,j]
            pcs = pcsM[i,j]
            k_expr = k_exprM[i,j]
            res = resM[i,j]
            currentclustname = paste0('c','_', i, '_',j, '__',
                                rowparam,
                                '_', rowvalues[i],
                                '_', colparam, '_',
                                colvalues[j])

            tic(paste0('Full BANKSY Run: ', currentclustname, '.'))
            clustering_results <- RunBanksy(gcm_raw, locs, lambda = lambda,
                                       pcs = pcs, k_expr = k_expr, res = res,
                                       clustering_name=currentclustname,
                                       zScaleRows = zScaleRows, lieden_iters=lieden_iters,
                                       zScaleBeforeAveraging = zScaleBeforeAveraging,
                                       zScaleOwnAfterAveraging = zScaleOwnAfterAveraging,
                                       zScaleNbrAfterAveraging = zScaleNbrAfterAveraging,
                                       banksyMatrices = banksyMatrices,
                                       multiple_datasets = multiple_datasets,
                                       compute_own_and_nbr_umap = compute_own_and_nbr_umap,
                                       compute_3D_umap = compute_3D_umap,
                                       cells_split_by_dataset = cells_split_by_dataset,
                                       instructions = instructions)
            # apr 28-- i just forgot to add cells_split_by_dataset to this. that is why the neighbour matrix was null.
            clusteringList[[i,j]] <- clustering_results$clustering[,c('cell_ID',
                                                                currentclustname),
                                                             with = FALSE]
            umapList[[i,j]] <- clustering_results$umap
            clusteringNames[[i,j]] = currentclustname
            parameterSetList[[i,j]] = clustering_results$parameterSet
            banksyParamsList[[i,j]] = clustering_results$banksy.params
            toc()
        }
    }
    sweepList = list(clusteringList=clusteringList,
                   umapList=umapList,
                   clusteringNames = clusteringNames,
                   parameterSetList = parameterSetList,
                   samplegobject = copy(clustering_results$gobject),
                   ownExpScaledMatrix = clustering_results$ownExpScaledMatrix,
                   nbrExpScaledMatrix = clustering_results$nbrExpScaledMatrix,
                   banksyParamsList = banksyParamsList
    )
    return(sweepList)
}
