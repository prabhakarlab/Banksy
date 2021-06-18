# plotUMAP
# plotSpatialDims

#' Plot UMAP
#'
#' @param bank BanksyObject
#' @param params clustering run to plot by
#' @param col cluster colours
#' @param legend cluster legend
#' @param main title
#' @param pt.size size of points
#' @param main.size size of title
#'
#' @importFrom ggplot2 ggplot geom_point aes xlab ylab theme_minimal ggtitle
#'   scale_color_manual theme element_text
#' @importFrom plyr mapvalues
#'
#' @return NULL
#'
#' @export
plotUMAP <- function(bank, params,
                     col = NULL, legend = TRUE,
                     main = NULL, pt.size = 1, main.size = 5) {

  clustnames <- names(bank@meta.data)
  clustnames <- clustnames[grep('^res', clustnames)]
  if (!(params %in% clustnames)) {
    stop(paste0('Invalid parameters. One of ',
                paste(clustnames, collapse = ' ')))
  }

  lam <- gsub('.*lam|_k.*', '', params)
  umap_name <- paste0('umap_', lam)

  clusters <- bank@meta.data[[params]]
  n <- length(unique(clusters))

  if (is.null(col)) {
    plotCols <- bank@meta.data[[paste0('col_', params)]]
    if (is.null(plotCols)) plotCols <- getPalette(n)
    else {
      colMap <- bank@meta.data[!duplicated(bank@meta.data[[params]]),]
      plotCols <- colMap[[paste0('col_', params)]][order(colMap[[params]])]
    }
  } else {
    if (length(col) < n) stop('Not enough colours')
    plotCols <- col[seq_len(n)]
  }

  umap_1 <- umap_2 <- Clust <- NULL
  plotData <- as.data.frame(bank@dim.reduction[[umap_name]])
  plotData <- cbind(plotData, Clust = as.factor(clusters))
  names(plotData)[seq_len(2)] <- c('umap_1', 'umap_2')

  p <- ggplot(plotData, aes(x=umap_1, y=umap_2, col=Clust)) +
    xlab('umap 1') + ylab('umap 2') + theme_minimal() + geom_point(size=pt.size) +
    scale_color_manual(values=plotCols)

  if (!legend) p <- p + theme(legend.position = 'none')
  if (is.null(main)) p <- p + ggtitle(gsub('_', ' ', params)) +
    theme(plot.title = element_text(size = main.size))
  return(p)
}

#' Plot Spatial dims
#'
#' @param bank BanksyObject
#' @param params clustering run to plot by
#' @param dataset if multiple dataset are run
#' @param col colours
#' @param main title
#' @param pt.size size of points
#' @param main.size size of title
#'
#' @importFrom ggplot2 ggplot geom_point aes xlab ylab theme_minimal ggtitle
#' @importFrom plyr mapvalues
#'
#' @return NULL
#'
#' @export
plotSpatialDims <- function(bank, params, dataset = NULL,
                            col = NULL, main = NULL,
                            pt.size = 1, main.size = 5) {

  clustnames <- names(bank@meta.data)
  clustnames <- clustnames[grep('^res', clustnames)]
  if (!(params %in% clustnames)) {
    stop(paste0('Invalid parameters. One of ',
                paste(clustnames, collapse = ' ')))
  }

  if (is.list(bank@own.expr)) {
    if (is.null(dataset)) {
      warning('No dataset specified. Using first dataset.')
      dataset <- 1
    }
    plotData <- bank@cell.locs[[dataset]]
  } else {
    plotData <- bank@cell.locs
  }

  clusters <- bank@meta.data[[params]][bank@meta.data$cell_ID %in%
                                         rownames(plotData)]
  n <- max(clusters)
  if (is.null(col)) {
    plotCols <- bank@meta.data[[paste0('col_', params)]]
    if (is.null(plotCols)) {
      plotCols <- plyr::mapvalues(clusters,
                                  from = seq_len(n),
                                  to = getPalette(n),
                                  warn_missing = FALSE)
    }
  } else {
    if (length(col) < n) stop('Not enough colours')
    plotCols <- plyr::mapvalues(clusters,
                                from = seq_len(n),
                                to = col[seq_len(n)],
                                warn_missing = FALSE)
  }
  dimx <- dimy <- NULL
  names(plotData)[seq_len(2)] <- c('dimx','dimy')
  p <- ggplot(plotData, aes(x=dimx, y=dimy)) +
    geom_point(color=plotCols, size=pt.size) +
    xlab('x coordinates') + ylab('y coordinates') + theme_minimal()
  if (is.null(main)) p <- p + ggtitle(gsub('_', ' ', params)) +
    theme(plot.title = element_text(size = main.size))

  return(p)
}

