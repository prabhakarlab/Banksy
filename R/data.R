
#' Mouse Hippocampus VeraFISH data
#'
#' This dataset comprises VeraFISH profiling of cells in the mouse hippocampus.
#' Gene expression and cell centroids for 7,944 cells and 124 genes in 2 spatial
#' dimensions are provided.
#'
#' @format A list with 2 entries:
#' \describe{
#' \item{expression}{(matrix) gene expression matrix}
#' \item{locations}{data.frame} cell centroids in 2D)
#' }
#'
"hippocampus"

#' Mouse cerebellum Slide-seq V2 data
#'
#' This dataset comprise Slide-seq V2 profiling of cells of the mouse
#' cerebellum. Gene expression and cell centroids for 9,953 cells and 8,365
#' genes in 2 spatial dimensions are provided.
#'
#' @format A list with 2 entries:
#' \describe{
#' \item{expression}{(matrix) gene expression matrix}
#' \item{locations}{(data.frame) cell centroids in 2D}
#' }
#'
#' @source \url{https://singlecell.broadinstitute.org/single_cell/study/SCP948}
#'
"cerebellum"

#' Mouse Hypothalamus MERFISH data
#'
#' This dataset comprises MERFISH profiling of cells in the pre-optic region
#' of the mouse hypothalamus (Moffitt et al. (2018), Molecular, spatial, and
#' functional single-cell profiling of the hypothalamic preoptic region.
#' doi: 10.1126/science.aau5324). The data is subset to 2 z-planes from a single
#' animal (Animal 1, bregmas 0.21 and 0.26). Gene expression, cell centroids
#' and metadata for 11,162 cells and 161 genes in 3 spatial dimensions are
#' provided.
#'
#' @format A list with 3 entries:
#' \describe{
#' \item{expression}{(matrix) gene expression matrix}
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
#' \item{expression}{(matrix) gene expression matrix}
#' \item{locations}{(data.frame) cell centroids in 2D}
#' \item{metadata}{(data.frame) metadata for cells, including HMRF domain}
#' \item{spatialgenes}{(character) spatial genes for HMRF analysis}
#' }
#'
#' @source \url{https://bitbucket.org/qzhudfci/smfishhmrf-r/src}
#'
"visualcortex"

#' Human dorsolateral prefrontal cortex 10x Visium data (sample 151673)
#'
#' This dataset comprises 10x Visium profiling of the human DLPFC (Maynard et
#' al. (2021) Transcriptome-scale spatial gene expression in the human
#' dorsolateral prefrontal cortex. doi.org: 10.1038/s41593-020-00787-0) for
#' subject 3 sample 151673. Gene expression for 3,639 spots in 2 dimensions was
#' count normalized and subset to 7,319 genes obtained by taking the union of
#' the top 3,000 genes in each sample for subject 3.
#'
#' @format A list with 3 entries:
#' \describe{
#' \item{expression}{(matrix) gene expression matrix}
#' \item{locations}{(data.frame) spot locations in 2D}
#' }
#'
#' @source \url{http://research.libd.org/spatialLIBD/}
#'
"dlpfc151673"

#' Human dorsolateral prefrontal cortex 10x Visium data (sample 151674)
#'
#' This dataset comprises 10x Visium profiling of the human DLPFC (Maynard et
#' al. (2021) Transcriptome-scale spatial gene expression in the human
#' dorsolateral prefrontal cortex. doi.org: 10.1038/s41593-020-00787-0) for
#' subject 3 sample 151674. Gene expression for 3,673 spots in 2 dimensions was
#' count normalized and subset to 7,319 genes obtained by taking the union of
#' the top 3,000 genes in each sample for subject 3.
#'
#' @format A list with 3 entries:
#' \describe{
#' \item{expression}{(matrix) gene expression matrix}
#' \item{locations}{(data.frame) spot locations in 2D}
#' }
#'
#' @source \url{http://research.libd.org/spatialLIBD/}
#'
"dlpfc151674"

#' Human dorsolateral prefrontal cortex 10x Visium data (sample 151675)
#'
#' This dataset comprises 10x Visium profiling of the human DLPFC (Maynard et
#' al. (2021) Transcriptome-scale spatial gene expression in the human
#' dorsolateral prefrontal cortex. doi.org: 10.1038/s41593-020-00787-0) for
#' subject 3 sample 151675. Gene expression for 3,592 spots in 2 dimensions was
#' count normalized and subset to 7,319 genes obtained by taking the union of
#' the top 3,000 genes in each sample for subject 3.
#'
#' @format A list with 3 entries:
#' \describe{
#' \item{expression}{(matrix) gene expression matrix}
#' \item{locations}{(data.frame) spot locations in 2D}
#' }
#'
#' @source \url{http://research.libd.org/spatialLIBD/}
#'
"dlpfc151675"

#' Human dorsolateral prefrontal cortex 10x Visium data (sample 151676)
#'
#' This dataset comprises 10x Visium profiling of the human DLPFC (Maynard et
#' al. (2021) Transcriptome-scale spatial gene expression in the human
#' dorsolateral prefrontal cortex. doi.org: 10.1038/s41593-020-00787-0) for
#' subject 3 sample 151676. Gene expression for 3,460 spots in 2 dimensions was
#' count normalized and subset to 7,319 genes obtained by taking the union of
#' the top 3,000 genes in each sample for subject 3.
#'
#' @format A list with 3 entries:
#' \describe{
#' \item{expression}{(matrix) gene expression matrix}
#' \item{locations}{(data.frame) spot locations in 2D}
#' }
#'
#' @source \url{http://research.libd.org/spatialLIBD/}
#'
"dlpfc151676"
