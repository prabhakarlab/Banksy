
#' Mouse Hippocampus VeraFISH data
#'
#' This dataset comprises VeraFISH profiling of cells in the mouse hippocampus.
#' Gene expression and cell centroids for 7,944 cells and 124 genes in 2 spatial
#' dimensions are provided.
#'
#' @format A list with 2 entries:
#' \describe{
#' \item{expression}{(matrix) gene-cell matrix}
#' \item{locations}{data.frame} cell centroids in 2D)
#' }
#'
"hippocampus"

#' Mouse brain organoid MERFISH data
#'
#' This data consists of 3 separate data sets comprising MERFISH profiling of
#' mouse brain organoids. Gene expression and cell centroids for 3,600 cells,
#' at least 557 genes in 2 spatial dimensions are provided.
#'
#' @format A list with 2 entries:
#' \describe{
#' \item{expression}{(list) a list with 3 entries, each a gene-cell matrix}
#' \item{locations}{list} a list with 3 entries with corresponding cell
#'   centroids in 2D)
#' }
#'
"organoid"

#' Mouse cerebellum Slide-seq V2 data
#'
#' This data set comprise Slide-seq V2 profiling of cells of the mouse
#' cerebellum. Gene expression and cell centroids for 9,953 cells and 8,365
#' genes in 2 spatial dimensions are provided.
#'
#' @format A list with 2 entries:
#' \describe{
#' \item{expression}{(matrix) gene-cell matrix}
#' \item{locations}{(data.frame) cell centroids in 2D}
#' }
#'
#' @source \url{https://singlecell.broadinstitute.org/single_cell/study/SCP948}
#'
"cerebellum"

#' Mouse Hypothalamus MERFISH data
#'
#' This data set comprises MERFISH profiling of cells in the pre-optic region
#' of the mouse hypothalamus (Moffitt et al. (2018), Molecular, spatial, and
#' functional single-cell profiling of the hypothalamic preoptic region.
#' doi: 10.1126/science.aau5324). The data is subset to 2 z-planes from a single
#' animal (Animal 1, bregmas 0.21 and 0.26). Gene expression, cell centroids
#' and metadata for 11,162 cells and 161 genes in 3 spatial dimensions are
#' provided.
#'
#' @format A list with 3 entries:
#' \describe{
#' \item{expression}{(matrix) gene-cell matrix}
#' \item{locations}{(data.frame) cell centroids in 3D}
#' \item{metadata}{(data.frame) metadata for cells, including cell class,
#'   neuronal cluster, animal ID, behavioral class}
#' }
#'
#' @source \url{https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248}
#'
"hypothalamus"


#' Mouse visual cortex SeqFISH data
#'
#' This dataset comprises SeqFISH profiling of cells in the mouse visual cortex
#' (Qian Zhu et al. (2018) Identification of spatially associated subpopulations
#' by combining scRNA-seq and sequential florescence in situ hybridization data.
#' doi: 10.1038/nbt.4260). Gene expression, cell centroids and metadata are
#' provided for 1,597 cells, 125 genes in 2 spatial dimensions are provided.
#' In addition, the spatial genes used for HMRF analysis in the
#' reference are included.
#'
#' @format A list with 3 entries:
#' \describe{
#' \item{expression}{(matrix) gene-cell matrix}
#' \item{locations}{(data.frame) cell centroids in 2D}
#' \item{metadata}{(data.frame) metadata for cells, including HMRF domain}
#' \item{spatialgenes}{(character) spatial genes for HMRF analysis}
#' }
#'
#' @source \url{https://bitbucket.org/qzhudfci/smfishhmrf-r/src}
#'
"visualcortex"




