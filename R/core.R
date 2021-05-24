# RunBanksy
# SweepBanksy

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
#' @param k_expr k_expr
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
