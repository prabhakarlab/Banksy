#' Mouse Hippocampus VeraFISH data
#'
#' This dataset comprises VeraFISH profiling of cells in the mouse hippocampus.
#' Gene expression and cell centroids for 10,944 cells and 129 genes in 2 
#' spatial dimensions are provided.
#'
#' @format A list with 2 entries:
#' \describe{
#' \item{expression}{(matrix) gene expression matrix}
#' \item{locations}{(data.frame) cell centroids in 2D}
#' }
#'
#' @usage data(hippocampus)
#'
"hippocampus"

#' An unrealistic simulation of spatially-resolved omics data.
#'
#' This dataset comprises gene expression and spatial coordinates for 50 genes
#' and 308 cells from 4 clusters (\code{rings$clusters}). Generate this dataset
#' with \cr
#'
#' \code{set.seed(2023)} \cr
#' \code{rings <- simulateDataset(n_cells=300, n_rings=4, n_genes=50)} \cr
#'
#' @format A SpatialExperiment object.
#'
#' @usage data(rings)
#'
"rings"
