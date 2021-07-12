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
  all.cols <- c(pals::kelly()[-c(1,3)],
                pals::glasbey(),
                pals::polychrome())
  all.cols <- as.character(all.cols)
  return(all.cols[seq_len(n)])
}

getAssay <- function(bank, assay, dataset, lambda) {

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

  } else if (!is.na(pmatch(assay, 'banksy'))) {

    if (is.null(lambda)) stop('Lambda not specified.')

    if (is.list(bank@own.expr)) {
      if (is.null(dataset)) {
        message('Dataset not specified. Choosing first dataset.')
        index <- 1
      } else if (!(dataset %in% names(bank@own.expr))) {
        stop(paste0('Dataset', dataset, ' not found.'))
      }
      own <- bank@own.expr[[dataset]]
      nbr <- bank@nbr.expr[[dataset]]
      mat <- rbind(sqrt(1-lambda)*own,
                     sqrt(lambda)*nbr)

    } else {
      mat <- getBanksyMatrix(bank, lambda)$expr
    }

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
    col <- c('magenta', 'black', 'yellow')
  }
  if (is.null(col.breaks)) {
    lim.high <- quantile(mat, 0.9)
    lim.low <- quantile(mat, 0.1)
    col.breaks <- c(lim.low, 0, lim.high)
  }

  col.fun <- colorRamp2(col.breaks, col)

  return(col.fun)
}














