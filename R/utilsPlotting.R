# Index of helper functions and where they are called
# init ---------------------- ConnectClusters
# getPalette ---------------- Plotting
# getAssay ------------------ Plotting Heatmap
# getHeatmapPalette --------- Plotting Heatmap
# getCellAnnotation --------- Plotting Heatmap
# appendBarplots ------------ Plotting Heatmap
# sampleMatrix -------------- Plotting Heatmap

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

getAssay <- function(bank, assay, dataset, lambda, cells, features) {

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
  if (is.null(features)) features <- rownames(mat)
  if (is.null(cells)) cells <- colnames(mat)
  features <- c(features, paste0(features, '.nbr'))
  mat <- mat[rownames(mat) %in% features, colnames(mat) %in% cells]
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
    lim.low <- -1
    lim.high <- 1
    col.breaks <- c(lim.low, 0, lim.high)
  }

  col.fun <- colorRamp2(col.breaks, col)

  return(col.fun)
}


getClusterColors <- function(x) {

  clusters <- unique(x)
  if (is.numeric(x)) {
    n <- max(x)
    cluster.cols <- getPalette(n)[sort(unique(clusters))]
  } else if (is.character(x)) {
    n <- length(clusters)
    cluster.cols <- getPalette(n)
  }
  names(cluster.cols) <- clusters

  return(cluster.cols)
}

#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @importFrom grid gpar
getCellAnnotation <- function(bank, mat,
                              annotate.by,
                              annotation.name,
                              annotation.size,
                              annotation.pos,
                              order.by) {

  if (is.null(annotate.by) | !all(annotate.by %in% names(bank@meta.data))) {
    stop('Specify a valid metadata column to annotate by')
  }

  if (is.null(order.by)) {
    order.by <- annotate.by[1]
    if (length(annotate.by) > 1) {
      message(paste0('Metadata to order cells by not specified. Using ', order.by))
    }
  } else if (!(order.by %in% names(bank@meta.data))) {
    stop('Specify a valid metadata column to order by.')
  }

  ## Order by clusters
  mdata <- bank@meta.data
  mdata <- mdata[order(mdata[[order.by]]),]
  mdata <- mdata[mdata$cell_ID %in% colnames(mat), ]
  cell.order <- mdata$cell_ID
  mdata <- mdata[,annotate.by,drop=FALSE]
  cell.split <- mdata[[order.by]]

  ## Generate simple annotations
  hc <- lapply(mdata, getClusterColors)
  mdata[] <- lapply(mdata, factor)
  ha <- HeatmapAnnotation(df = mdata,
                          show_annotation_name = annotation.name,
                          annotation_name_gp = gpar(fontsize = annotation.size),
                          annotation_name_side = annotation.pos,
                          col = hc)

  return(list(anno = ha,
              cell.order = cell.order,
              cell.split = cell.split))
}

#' @importFrom ComplexHeatmap anno_barplot HeatmapAnnotation `%v%`
#' @importFrom grid gpar
appendBarplots <- function(bank, mat, cell.order,
                           barplot.by,
                           barplot.border,
                           barplot.width,
                           annotation.name,
                           annotation.size,
                           annotation.pos,
                           heatmap) {

  if (!all(barplot.by %in% names(bank@meta.data))) {
    stop('Specify valid metadata columns to plot barplots.')
  }

  mdata <- bank@meta.data
  mdata <- mdata[cell.order, ]

  nbar <- length(barplot.by)
  for (i in seq_len(nbar)) {
    bar.anno <- anno_barplot(mdata[[barplot.by[i]]],
                             border = barplot.border,
                             bar_width = barplot.width)
    bar <- HeatmapAnnotation(bar = bar.anno,
                             show_annotation_name = annotation.name,
                             annotation_name_gp = gpar(fontsize = annotation.size),
                             annotation_name_side = annotation.pos)
    bar@anno_list$bar@label <- barplot.by[i]
    heatmap <- bar %v% heatmap
  }

  return(heatmap)
}

sampleMatrix <- function(mat, max.cols, seed) {
  set.seed(seed)
  message(paste0('Sampling ', max.cols, ' columns with seed ', seed))
  message(paste0('Keeping ', round(max.cols/ncol(mat) * 100,2), '% of cells'))
  keep <- sample(colnames(mat), min(length(colnames(mat)), max.cols))
  mat <- mat[, keep]
  return(mat)
}
