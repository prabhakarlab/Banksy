# clusteringsAsDF
# updateGobjectUmapClust
# update_indiv_umap

#' @importFrom data.table setDF
clusteringsAsDF<- function(sweepListArray){
  # convert the clustering results from a 2d list array of data tables to a data table.

  nr = dim(sweepListArray$clusteringList)[1]
  nc = dim(sweepListArray$clusteringList)[2]

  clusteringsDT = sweepListArray$clusteringList[[1,1]]
  for (i in 1:nr){
    for (j in 1:nc){
      if (!(i==1 && j==1)){
        incomingClustering = sweepListArray$clusteringList[[i, j]]
        incomingClusteringName = sweepListArray$clusteringNames[[i, j]]
        if (!all(incomingClustering$cell_ID == clusteringsDT$cell_ID)){
          stop(paste0('Not all the cell IDs in the clustering (', i, ',
                    ', j, ') are the same as those the previous clusterings.'))
        } else {
          clusteringsDT = cbind(clusteringsDT, incomingClustering[,incomingClusteringName, with=FALSE])
        }
      }
    }
  }# so add row wise.
  rownames(clusteringsDT) = sweepListArray$clusteringList[[1,1]]$cell_ID
  # do type conversion from numeric to character.
  clusteringsDF = clusteringsDT
  rownames(clusteringsDF) = sweepListArray$clusteringList[[1,1]]$cell_ID
  for (i in 1:nrow(clusteringsDT)){
    for (j in 2:ncol(clusteringsDT)){ #i think the first colunm is cellID, dont convert
      clusteringsDF[[i, j]] = as.character(clusteringsDT[[i,j]])
    }
  }
  setDF(clusteringsDF)
  rownames(clusteringsDF) <- rownames(clusteringsDT)
  clusteringsDF = clusteringsDF[,!(names(clusteringsDF)%in%c('cell_ID'))]
  return(clusteringsDF)
}

#' Updates Giotto object UMAP clusters with sweep results and connected output
#'
#' @param gobject Giotto object
#' @param sweepresults returned by SweepResults
#' @param pairConnectedOutput returned by gridConnect
#' @param compute_own_and_nbr_umap no
#' @param compute_3D_umap no
#'
#' @return updated Giotto object
#'
#' @export
updateGobjectUmapClust <- function(gobject, sweepresults, pairConnectedOutput, compute_own_and_nbr_umap = FALSE,
                                   compute_3D_umap=FALSE){
  # !!april20
  if (!(all(gobject@cell_metadata$cell_ID==rownames(rownames(pairConnectedOutput$clusterings))))){
    stop('The cell IDs in the metadata and the clustering results data frame do not match. Something has gone wrong.')
  } else {
    message('The cell IDs in the metadata and the clustering results data frame match; appending clustering results')
    gobject@cell_metadata = cbind(gobject@cell_metadata, pairConnectedOutput$clusterings) # only do this once!
  }
  for (i in 1:prod(dim(sweepresults$clusteringNames))){
    message(sweepresults$clusteringNames[[i]]) # counts down columns.
    current_cluster_parameters = sweepresults$parameterSetList[[i]][grep('cluster',
                                                                         names(sweepresults$parameterSetList[[i]]))][[1]]
    current_clustering_name <- sweepresults$clusteringNames[[i]]
    ##
    umap_count = 1
    current_umap_type = paste0('joint_2D_', current_clustering_name)
    gobject = update_indiv_umap(gobject,sweepresults,
                                sweepNum=i,
                                umap_count=umap_count,
                                current_umap_type = current_umap_type)
    umap_count = umap_count+1


    if( compute_own_and_nbr_umap == TRUE ){
      current_umap_type = paste0('own_2D_', current_clustering_name)
      gobject = update_indiv_umap(gobject,sweepresults,
                                  sweepNum=i,
                                  umap_count=umap_count,
                                  current_umap_type = current_umap_type)
      umap_count = umap_count+1

      current_umap_type = paste0('nbr_2D_', current_clustering_name)
      gobject = update_indiv_umap(gobject,sweepresults,
                                  sweepNum=i,
                                  umap_count=umap_count,
                                  current_umap_type = current_umap_type)
      umap_count = umap_count+1
    }
    if (compute_3D_umap==TRUE){
      current_umap_type = paste0('joint_3D_', current_clustering_name)
      gobject = update_indiv_umap(gobject,sweepresults,
                                  sweepNum=i,
                                  umap_count=umap_count,
                                  current_umap_type = current_umap_type)
      umap_count = umap_count+1

      if( compute_own_and_nbr_umap == TRUE ){
        current_umap_type = paste0('own_3D_', current_clustering_name)
        gobject = update_indiv_umap(gobject,sweepresults,
                                    sweepNum=i,
                                    umap_count=umap_count,
                                    current_umap_type = current_umap_type)
        umap_count = umap_count+1

        current_umap_type = paste0('nbr_3D_', current_clustering_name)
        gobject = update_indiv_umap(gobject,sweepresults,
                                    sweepNum=i,
                                    umap_count=umap_count,
                                    current_umap_type = current_umap_type)
        umap_count = umap_count+1

      }
    }
    #
    # gobject@dimension_reduction[["cells"]][["umap"]][[current_clustering_name]] <- sweepresults$umapList[[i]][[current_clustering_name]]
    # parameters_list = gobject@parameters
    # number_of_rounds = length(parameters_list)
    # update_name = paste0(number_of_rounds, "_umap")
    # gobject@parameters[[update_name]] = current_umap_parameters
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds, "_cluster")
    gobject@parameters[[update_name]] = current_cluster_parameters
  }
  return(gobject)
}

update_indiv_umap <- function(gobject,sweepresults, sweepNum=1, umap_count=1, current_umap_type = NULL){
  current_umap_parameters = sweepresults$parameterSetList[[sweepNum]][grep('umap',
                                                                           names(sweepresults$parameterSetList[[sweepNum]]))][[umap_count]]
  if (current_umap_parameters["name for umap"]!=current_umap_type){
    stop('The UMAP being used is not correct. Check what went wrong.')
  }
  gobject@dimension_reduction[["cells"]][["umap"]][[current_umap_type]] <- sweepresults$umapList[[sweepNum]][[current_umap_type]]
  parameters_list = gobject@parameters
  number_of_rounds = length(parameters_list)
  update_name = paste0(number_of_rounds, "_umap")
  gobject@parameters[[update_name]] = current_umap_parameters
  return(gobject)
}


