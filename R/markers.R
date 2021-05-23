# getMarkerGenes
# orderMarkerGenes
# filterMarkers

#' Get marker genes
#'
#' @param gobject Giotto object
#' @param sweepresults returned by SweepBanksy
#' @param pairConnectedOutput returned by gridConnect
#' @param numMarkers number of markers
#' @param upper_clip upp
#' @param lower_clip low
#' @param markerMethod method used
#'
#' @importFrom Giotto findMarkers_one_vs_all
#' @importFrom utils head
#'
#' @return List of marker genes
#'
#' @export
getMarkerGenes <- function(gobject=NULL, sweepresults=NULL, pairConnectedOutput = NULL,
                           numMarkers = 2, upper_clip = 2, lower_clip = -2,
                           markerMethod = c('scran', 'scran_and_gini','gini')){
  # !! april20
  gridSize = dim(sweepresults$clusteringNames)
  topgene2a = NULL
  for (i in 1:prod(gridSize)){
    # print(i)
    # gobj_hm = copy(gobject)
    currName = names(pairConnectedOutput$clusterings)[i]
    # gene_names = rownames(gobject@norm_scaled_expr)[1:(length(rownames(gobject@norm_scaled_expr))/2)]
    # d_gene <- dist(as.matrix(gobject@norm_scaled_expr[gene_names,]))
    # hc_gene <- hclust(d_gene)
    # gene_order = hc_gene$labels[hc_gene$order]
    # gene_order_banksy = c(gene_order, paste0(gene_order, '.nbr'))

    gobject <- addClippedZCustomExpr(gobject,pairConnectedOutput,sweepresults,
                                     sweepNum = i, upper_clip = upper_clip, lower_clip = lower_clip)
    if (markerMethod == 'scran_and_gini'){

      markers_banksy_gini = findMarkers_one_vs_all(gobject,method = 'gini',
                                                   expression_values = 'custom',
                                                   cluster_column = currName,min_genes = 5, rank_score = 2)
      markers_banksy_scran = findMarkers_one_vs_all(gobject,
                                                    method = 'scran',
                                                    expression_values = 'custom',
                                                    cluster_column = currName,
                                                    min_genes = 5, rank_score = 2)

      .SD <- NULL # global binding
      topgini_banksy_genes_2 = unique(markers_banksy_gini[, head(.SD, numMarkers), by = 'cluster']$genes)
      topscran_banksy_genes_2 = unique(markers_banksy_scran[, head(.SD, numMarkers), by = 'cluster']$genes)
      topgene2a = union(topgene2a, union(topgini_banksy_genes_2, topscran_banksy_genes_2))
    } else if (markerMethod == 'scran'){
      markers_banksy_scran = findMarkers_one_vs_all(gobject,
                                                    method = 'scran',
                                                    expression_values = 'custom',
                                                    cluster_column = currName,
                                                    min_genes = 5, rank_score = 2)
      topscran_banksy_genes_2 = unique(markers_banksy_scran[, head(.SD, numMarkers), by = 'cluster']$genes)
      topgene2a = union(topgene2a, topscran_banksy_genes_2)
    } else if (markerMethod == 'gini'){
      markers_banksy_gini = findMarkers_one_vs_all(gobject,method = 'gini',
                                                   expression_values = 'custom',
                                                   cluster_column = currName,min_genes = 5, rank_score = 2)
      topgini_banksy_genes_2 = unique(markers_banksy_gini[, head(.SD, numMarkers), by = 'cluster']$genes)
      topgene2a = union(topgene2a, topgini_banksy_genes_2)
    }
  }

  topgene2b = unique(gsub(".nbr", replacement = '', topgene2a))
  topgene2c = c(topgene2b, paste0(topgene2b, '.nbr'))
  d_gene_marker <- dist(as.matrix(gobject@norm_scaled_expr[topgene2b,]))
  hc_gene_marker <- hclust(d_gene_marker)
  gene_order_marker = hc_gene_marker$labels[hc_gene_marker$order]
  gene_order_banksy_marker2 = c(gene_order_marker, paste0(gene_order_marker, '.nbr'))
  # !! call orderMarkerGenes to do this.
  return(list(marker2 = topgene2c, markerOrder2 = gene_order_banksy_marker2))
}

#' @importFrom stats dist hclust
orderMarkerGenes <- function(gene_cell_matrix, markers){
  d_gene_marker <- dist(as.matrix(gene_cell_matrix[markers,]))
  hc_gene_marker <- hclust(d_gene_marker)
  gene_order_marker = hc_gene_marker$labels[hc_gene_marker$order]
  gene_order_banksy_marker2 = c(gene_order_marker, paste0(gene_order_marker, '.nbr'))
  return(list(marker2 = c(markers, paste0(markers, '.nbr')), markerOrder2 = gene_order_banksy_marker2))
}

filterMarkers <- function(markers, clusteringLabels, expression_matrix){
  # remove markers that are not at least expressed in 10 percent of cells in at least one cluster.

}
