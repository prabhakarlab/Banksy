#'
#' The Banksy package
#'
#' Banksy is a library and R package for network analysis.
#'
#' @rdname Banksy-package
#' @name Banksy-package
#' @keywords internal
#' @aliases Banksy-package Banksy
#' @docType package
#'
#' @section Description:
#' Banksy is an R package that incorporates spatial information to cluster 
#' cells in a feature space (e.g. gene expression). To incorporate spatial
#' information, BANKSY computes the mean neighborhood expression and azimuthal
#' Gabor filters that capture gene expression gradients. These features are
#' combined with the cell's own expression to embed cells in a 
#' neighbor-augmented product space which can then be clustered, allowing for
#' accurate and spatially-aware cell typing and tissue domain segmentation. 
#' 
#' @section Details:
#' For a quick start to the package, please refer to the GitHub page at 
#' \url{https://github.com/prabhakarlab/Banksy}. For in-depth guides to package 
#' functionality and use cases, refer to the package webpage at 
#' \url{https://prabhakarlab.github.io/Banksy}. 
#' 
"_PACKAGE"