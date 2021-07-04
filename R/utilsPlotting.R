# Index of helper functions and where they are called
# init ---------------------- ConnectClusters
# getPalette ---------------- Plotting
# getAssay ------------------ Plotting Heatmap
# getHeatmapPalette --------- Plotting Heatmap

init <- function(n, idx, val) {
  x <- rep(0, n)
  x[idx] <- val
  x
}

#' @importFrom pals kelly glasbey polychrome
getPalette <- function(n) {
  all.cols <- as.character(pals::kelly()[-c(1,3)],
                           pals::glasbey(),
                           pals::polychrome())
  return(all.cols[seq_len(n)])
}

getAssay <- function(bank, assay, dataset) {

  if (!is.na(pmatch(assay, c('own.expr', 'nbr.expr', 'custom.expr')))) {
    slot <- get(assay, mode = 'function')
    mat <- slot(bank)
    if (is.list(mat)) {
      if (is.null(dataset)) {
        message('Dataset not specified. Choosing first dataset.')
        mat <- mat[[1]]
      } else if (!(dataset %in% names(mat))) {
        stop(paste0('Dataset ', dataset, ' not found.'))
      } else {
        mat <- mat[[dataset]]
      }
    }
    if (is.null(mat)) stop('Assay chosen is empty.')
  } else {
    stop('Specify a valid assay.')
  }
  return(mat)
}

#' @importFrom circlize colorRamp2
#' @importFrom stats quantile
getHeatmapPalette <- function(mat, col, col.breaks) {

  if (!is.null(col) &
      !is.null(col.breaks) &
      (length(col) != length(col.breaks))) {
    stop('Unequal number of colours and breaks')
  }
  if (is.null(col)) {
    col <- c('blue', 'white', 'red')
  }
  if (is.null(col.breaks)) {
    lim.high <- quantile(mat, 0.9)
    lim.low <- quantile(mat, 0.1)
    col.breaks <- seq(lim.low, lim.high, length.out = length(col))
  }

  col.fun <- colorRamp2(col.breaks, col)

  return(col.fun)
}














