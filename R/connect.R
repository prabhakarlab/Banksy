# gridConnect
# pairConnect

#' Connects the colours in the grid of the parameter sweep
#'
#' @param sweepresults returned by SweepResults
#' @param rootnode root node used
#'
#' @return PairedConnectedOutput
#'
#' @export
gridConnect <- function(sweepresults, rootnode=c(5,2)){
  # !TODO: check what happens if the root node is at either edge / corner. seems like c(2, 1) does not work.

  clusteringsDF = clusteringsAsDF(sweepresults) ## clustering with cell ID as rowname
  gridSize = dim(sweepresults$clusteringNames)
  rootColumn = rootnode[2]
  rootRow = rootnode[1]
  rightLength = gridSize[2] - rootColumn # take this many steps rightwards
  leftLength = rootColumn - 1 # take this many steps leftwards
  upLength = rootRow - 1 # take this many steps upwards
  downLength = gridSize[1] - rootRow #take this many steps downwards

  # set up the initial clustering array. This will be updated every time.
  clusteringsDF_state = copy(clusteringsDF) # initialize
  usedLabels_state = NULL
  #### expand in the left direction
  if (leftLength != 0){
    # do the first one here.
    parentNodeIx = gridSize[2]*(rootRow-1)+rootColumn
    childNodeIx = parentNodeIx-1 # because we are counting row wise in the data frame of clusterings.
    pairConnectedOutput = pairConnect(clusteringsDF_state, childNode = colnames(clusteringsDF_state)[childNodeIx],
                                      parentNode = colnames(clusteringsDF_state)[parentNodeIx],
                                      usedLabels = usedLabels_state)
    clusteringsDF_state = pairConnectedOutput$clusterings
    usedLabels_state = pairConnectedOutput$usedLabels

    if (leftLength != 1){
      for (i in 2:leftLength){
        parentNodeIx = childNodeIx
        childNodeIx = parentNodeIx-1 # because we are counting row wise in the data frame of clusterings.
        pairConnectedOutput = pairConnect(clusteringsDF_state,
                                          childNode = colnames(clusteringsDF_state)[childNodeIx],
                                          parentNode = colnames(clusteringsDF_state)[parentNodeIx],
                                          usedLabels = usedLabels_state)
        clusteringsDF_state = pairConnectedOutput$clusterings
        usedLabels_state = pairConnectedOutput$usedLabels

      }
    }
  }
  # print('b')
  #### expand in the right direction .
  if (rightLength != 0){
    # do the first one here.
    parentNodeIx = gridSize[2]*(rootRow-1)+rootColumn
    childNodeIx = parentNodeIx+1 # because we are counting row wise in the data frame of clusterings.

    pairConnectedOutput = pairConnect(clusteringsDF_state,
                                      childNode = colnames(clusteringsDF_state)[childNodeIx],
                                      parentNode = colnames(clusteringsDF_state)[parentNodeIx],
                                      usedLabels = usedLabels_state)
    clusteringsDF_state = pairConnectedOutput$clusterings
    usedLabels_state = pairConnectedOutput$usedLabels

    if (rightLength != 1){
      for (i in 2:rightLength){
        parentNodeIx = childNodeIx
        childNodeIx = parentNodeIx+1 # because we are counting row wise in the data frame of clusterings.
        pairConnectedOutput = pairConnect(clusteringsDF_state,
                                          childNode = colnames(clusteringsDF_state)[childNodeIx],
                                          parentNode = colnames(clusteringsDF_state)[parentNodeIx],
                                          usedLabels = usedLabels_state)
        clusteringsDF_state = pairConnectedOutput$clusterings
        usedLabels_state = pairConnectedOutput$usedLabels

      }
    }
  }
  # print('d')
  # now expand out in the columns from this row:
  for (j in 1:gridSize[2]){
    # print(j)
    current_column_parent_Ix = gridSize[2]*(rootRow-1)+j
    if (upLength != 0){
      childNodeIx = current_column_parent_Ix - gridSize[2]
      pairConnectedOutput = pairConnect(clusteringsDF_state,
                                        childNode = colnames(clusteringsDF_state)[childNodeIx],
                                        parentNode = colnames(clusteringsDF_state)[current_column_parent_Ix],
                                        usedLabels = usedLabels_state)
      clusteringsDF_state = pairConnectedOutput$clusterings
      usedLabels_state = pairConnectedOutput$usedLabels

      if (upLength != 1){
        for (i in 2:upLength){
          parentNodeIx = childNodeIx
          childNodeIx = parentNodeIx - gridSize[2] # because we are counting row wise in the data frame of clusterings.
          pairConnectedOutput = pairConnect(clusteringsDF_state,
                                            childNode = colnames(clusteringsDF_state)[childNodeIx],
                                            parentNode = colnames(clusteringsDF_state)[parentNodeIx],
                                            usedLabels = usedLabels_state)
          clusteringsDF_state = pairConnectedOutput$clusterings
          usedLabels_state = pairConnectedOutput$usedLabels

        }
      }
    }
    if (downLength != 0){
      childNodeIx = current_column_parent_Ix + gridSize[2]
      pairConnectedOutput = pairConnect(clusteringsDF_state,
                                        childNode = colnames(clusteringsDF_state)[childNodeIx],
                                        parentNode = colnames(clusteringsDF_state)[current_column_parent_Ix],
                                        usedLabels = usedLabels_state)
      clusteringsDF_state = pairConnectedOutput$clusterings
      usedLabels_state = pairConnectedOutput$usedLabels

      if (downLength != 1){
        for (i in 2:downLength){
          parentNodeIx = childNodeIx
          childNodeIx = parentNodeIx + gridSize[2] # because we are counting row wise in the data frame of clusterings.
          pairConnectedOutput = pairConnect(clusteringsDF_state,
                                            childNode = colnames(clusteringsDF_state)[childNodeIx],
                                            parentNode = colnames(clusteringsDF_state)[parentNodeIx],
                                            usedLabels = usedLabels_state)
          clusteringsDF_state = pairConnectedOutput$clusterings
          usedLabels_state = pairConnectedOutput$usedLabels

        }
      }
    }
  }
  return(pairConnectedOutput)
}


#' @importFrom dplyr `%>%` count arrange
#'
pairConnect <- function(clusterings, parentNode, childNode, usedLabels=NULL, mappingMode = 'parentCentric'){
  ################# Combined Method ###############
  # Get both the parent centric and child centric transfers, and merge to get the best outcome.
  ######## child centric ############
  subset_clustering_results = clusterings[,c(parentNode, childNode)]
  rownames(subset_clustering_results) = rownames(clusterings)
  # subset_clustering_results[,mget(childNode)] =
  #   factor(paste0(subset_clustering_results[,mget(childNode)], '_'))
  subset_clustering_results[,childNode] = paste0(subset_clustering_results[,childNode], '_')

  current_linkage_freqs = subset_clustering_results %>% count(!!as.name(parentNode),
                                                              !!as.name(childNode),
                                                              name='freq') # need dplyr
  # ---- not just a simple majority, but actually over some threshold percentage ----
  proportion_threshold = 0.5
  clust1 = sort(as.character((unique(current_linkage_freqs[,1]))[[1]]))
  clust2 = sort(as.character((unique(current_linkage_freqs[,2]))[[1]]))
  two_relto_one_thresh = matrix(0L, nrow = length(clust1), ncol = length(clust2))
  newChildLabels1 = rep(NA, length(clust2))
  names(newChildLabels1) = clust2
  parentClust <- childClust <- NULL # binding for variable
  parentLabelFrequencies1 = data.frame(parentClust = character(),
                                       childCLust = character(),
                                       clustFreq = numeric())
  for (clust2_ix in 1:length(clust2)){
    current_cluster_rows = as.vector(current_linkage_freqs[,2]==clust2[clust2_ix])
    current_frequencies_tibble = current_linkage_freqs[current_cluster_rows,]
    # print(current_frequencies_tibble)
    current_frequencies = current_frequencies_tibble$freq
    maxrow = which.max(current_frequencies)
    sum_of_frequencies = sum(current_frequencies)
    max_of_frequencies = current_frequencies[maxrow]
    proportion_of_max = max_of_frequencies / sum_of_frequencies
    color_replaced_list = c()
    if (!(proportion_of_max < proportion_threshold)){
      parent_label_to_use = as.character(current_frequencies_tibble[maxrow, parentNode])
      child_label_replaced = clust2[clust2_ix]
      parentLabelFrequencies1 = rbind(parentLabelFrequencies1, data.frame(parentClust = as.character(parent_label_to_use),
                                                                          childClust = as.character(child_label_replaced),
                                                                          clustFreq = max_of_frequencies,
                                                                          stringsAsFactors=FALSE))
    }
  }
  uniqueParentsToTransfer = unique(parentLabelFrequencies1$parentClust)
  # print(parentLabelFrequencies1)
  for (uniqueParent in uniqueParentsToTransfer){
    currentParentFrequencies = parentLabelFrequencies1[parentLabelFrequencies1$parentClust==uniqueParent, ]
    maxrowparents = which.max(currentParentFrequencies$clustFreq)
    newChildLabels1[currentParentFrequencies[maxrowparents,]$childClust] = currentParentFrequencies[maxrowparents,]$parentClust
  }

  ############ do parent centric ##############
  subset_clustering_results = clusterings[,c(parentNode, childNode)]

  rownames(subset_clustering_results) = rownames(clusterings)
  # subset_clustering_results[,mget(childNode)] =
  #   factor(paste0(subset_clustering_results[,mget(childNode)], '_'))
  subset_clustering_results[,childNode] = paste0(subset_clustering_results[,childNode], '_')
  current_linkage_freqs = subset_clustering_results %>% count(!!as.name(parentNode),
                                                              !!as.name(childNode),
                                                              name='freq') # need dplyr
  # ---- not just a simple majority, but actually over some threshold percentage ----
  proportion_threshold = 0.4
  clust1 = sort(as.character((unique(current_linkage_freqs[,1]))[[1]]))
  clust2 = sort(as.character((unique(current_linkage_freqs[,2]))[[1]]))
  two_relto_one_thresh = matrix(0L, nrow = length(clust1), ncol = length(clust2))
  newChildLabels2 = rep(NA, length(clust2))
  names(newChildLabels2) = clust2
  parentLabelFrequencies2 = data.frame(parentClust = character(),
                                       childCLust = character(),
                                       clustFreq = numeric())

  # this is the part that is really different from the childCentric mode.
  # we go over each parent cluster one at a time.
  for (clust1_ix in 1:length(clust1)){
    current_cluster_rows = as.vector(current_linkage_freqs[,1]==clust1[clust1_ix])
    current_frequencies_tibble = current_linkage_freqs[current_cluster_rows,]
    # print(current_frequencies_tibble)
    current_frequencies = current_frequencies_tibble$freq
    maxrow = which.max(current_frequencies)
    sum_of_frequencies = sum(current_frequencies)
    max_of_frequencies = current_frequencies[maxrow]
    proportion_of_max = max_of_frequencies / sum_of_frequencies
    color_replaced_list = c()
    if (!(proportion_of_max < proportion_threshold)){
      parent_label_to_use = as.character(current_frequencies_tibble[maxrow, parentNode])
      child_label_replaced = as.character(current_frequencies_tibble[maxrow, childNode])
      parentLabelFrequencies2 = rbind(parentLabelFrequencies2, data.frame(parentClust = as.character(parent_label_to_use),
                                                                          childClust = as.character(child_label_replaced),
                                                                          clustFreq = max_of_frequencies,
                                                                          stringsAsFactors=FALSE))
    }
  }
  uniqueChildrenToReceive = unique(parentLabelFrequencies2$childClust)
  # print(parentLabelFrequencies2)
  for (uniqueChild in uniqueChildrenToReceive){
    currentChildFrequencies = parentLabelFrequencies2[parentLabelFrequencies2$childClust==uniqueChild, ]
    maxrowparents = which.max(currentChildFrequencies$clustFreq)
    newChildLabels2[currentChildFrequencies[maxrowparents,]$childClust] = currentChildFrequencies[maxrowparents,]$parentClust
  }

  # > parentLabelFrequencies1
  # parentClust childClust clustFreq
  # 1           2         1_      2030
  # 2           2        10_         3
  # 3           4        11_         4
  # 4           4         2_      3123
  # 5           1         3_      2328
  # 6           3         4_      1006
  # 7           1         6_       405
  # 8           6         7_       330
  # 9           7         9_       110
  # > parentLabelFrequencies2
  # parentClust childClust clustFreq
  # 1           1         3_      2328
  # 2           2         1_      2030
  # 3           4         2_      3123
  # 4           5         1_       353
  # 5           6         7_       330
  # 6           7         9_       110
  # 7           8         3_        25
  # 8           9         2_        21
  # >
  mergedParentLabelFrequencies = arrange(unique(rbind(parentLabelFrequencies1,
                                                      parentLabelFrequencies2)),
                                         parentClust)
  # parentClust childClust clustFreq
  # 1            1         3_      2328
  # 2            1         6_       405
  # 3            2         1_      2030
  # 4            2        10_         3
  # 5            3         4_      1006
  # 6            4        11_         4
  # 7            4         2_      3123
  # 8            5         1_       353
  # 9            6         7_       330
  # 10           7         9_       110
  # 11           8         3_        25
  # 12           9         2_        21

  # now:
  # for each parent, find the row with the max freq:
  parentLabelFrequencies3 = data.frame(parentClust = character(),
                                       childCLust = character(),
                                       clustFreq = numeric())
  for (parentID in unique(mergedParentLabelFrequencies$parentClust)){
    subsetByParent = mergedParentLabelFrequencies[mergedParentLabelFrequencies$parentClust==parentID,]
    current_frequencies = subsetByParent$clustFreq
    maxrow = which.max(current_frequencies)
    parentLabelFrequencies3 = rbind(parentLabelFrequencies3, subsetByParent[maxrow,])
  }

  # parentClust childClust clustFreq
  # 1            1         3_      2328
  # 3            2         1_      2030
  # 5            3         4_      1006
  # 7            4         2_      3123
  # 8            5         1_       353
  # 9            6         7_       330
  # 10           7         9_       110
  # 11           8         3_        25
  # 12           9         2_        21
  parentLabelFrequencies3 = arrange(parentLabelFrequencies3, childClust)
  parentLabelFrequencies4 = data.frame(parentClust = character(),
                                       childCLust = character(),
                                       clustFreq = numeric())
  for (childID in unique(parentLabelFrequencies3$childClust)){
    subsetByChild= parentLabelFrequencies3[parentLabelFrequencies3$childClust==childID,]
    current_frequencies = subsetByChild$clustFreq
    maxrow = which.max(current_frequencies)
    parentLabelFrequencies4 = rbind(parentLabelFrequencies4, subsetByChild[maxrow,])
  }
  parentLabelFrequencies4 = arrange(parentLabelFrequencies4, parentClust)
  # print(parentLabelFrequencies4)
  newChildLabels4 = rep(NA, length(clust2))
  names(newChildLabels4) = clust2

  # print(parentLabelFrequencies4)
  for (row_id in 1:nrow(parentLabelFrequencies4)){
    newChildLabels4[parentLabelFrequencies4$childClust[row_id]] = parentLabelFrequencies4$parentClust[row_id]
  }
  # newChildLabels4
  #################### final merging ##############
  numberOfNewLabelsNeeded = sum(is.na(newChildLabels4))
  usedLabels2 = union(as.character(clust1), as.character(usedLabels))
  usedLabelsNumeric = sort(as.numeric(usedLabels2))
  candidateLabels = 1:(length(usedLabelsNumeric)+numberOfNewLabelsNeeded) #try max (usedLabelsNumeric)?
  newLabels = sort(setdiff(candidateLabels, usedLabelsNumeric))[1:numberOfNewLabelsNeeded]
  newChildLabels4[is.na(newChildLabels4)] = as.character(newLabels)
  newUsedLabels = union(usedLabels2, newChildLabels4)
  # finally, update the child clustering with the new labels
  # also replace all instances of the labels in the cell labels in the child clustering.
  subset_clustering_results2 = subset_clustering_results
  for (i in 1:length(newChildLabels4)){
    ix_to_replace = subset_clustering_results2[,childNode]==names(newChildLabels4)[i]
    subset_clustering_results2[ix_to_replace, childNode] = newChildLabels4[i]
  }
  # Now update the original all_clustering_results
  updatedClusterings = clusterings
  updatedClusterings[, childNode] = subset_clustering_results2[,childNode]
  # output results
  pairConnectedOutput = list(clusterings=updatedClusterings,
                             usedLabels=newUsedLabels,
                             parentToChildMapping = newChildLabels4
  )
  return(pairConnectedOutput)




}

