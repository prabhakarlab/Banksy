#' Mouse Hippocampus VeraFISH data
#'
#' This dataset comprises VeraFISH profiling of cells in the mouse hippocampus.
#' Gene expression and cell centroids for 10,944 cells and 129 genes in 2
#' spatial dimensions are provided. For details on how this dataset was 
#' generated, refer to Supplementary Information section 2.2 of our 
#' \href{https://www.biorxiv.org/content/10.1101/2022.04.14.488259v1.supplementary-material}{preprint}. 
#'
#' @format A list with 2 entries:
#' \describe{
#' \item{expression}{(matrix) gene expression matrix}
#' \item{locations}{(data.frame) cell centroids in 2D}
#' }
#'
#' @usage data(hippocampus)
#' 
#' @return List with expression and locations
#'
"hippocampus"

#' An unrealistic simulation of spatially-resolved omics data.
#'
#' This dataset comprises gene expression and spatial coordinates for 50 genes
#' and 308 cells from 4 clusters (\code{rings$clusters}). See \code{system.file('scripts/rings.R', package='Banksy')} on how this dataset was generated. 
#'
#' @format A SpatialExperiment object.
#'
#' @usage data(rings)
#' 
#' @return A SpatialExperiment object
#'
"rings"
