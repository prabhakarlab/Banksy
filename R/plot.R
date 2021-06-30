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
  n <- max(clusters)
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
    xlab('umap 1') + ylab('umap 2') + theme_minimal() +
    geom_point(size=pt.size) +
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


#' Plot Heatmap (wrapper for ComplexHeatmap)
#'
#' @param bank BanksyObject
#' @param assay assay to plot heatmap
#' @param dataset dataset to plot heatmap
#' @param col colours to use in heatmap
#' @param col.breaks color breaks to use in heatmap (same number as col is
#'   specified)
#' @param cluster.row cluster rows
#' @param cluster.column cluster columns
#' @param row.dend draw row dendrograms
#' @param column.dend draw column dendrograms
#' @param cex.row row label size
#' @param cex.column column label size
#' @param annotate add annotation for cells
#' @param annotate.by metadata to annotate cells by
#' @param annotate.col colours to annotate cells by
#' @param name name of heatmap legend
#' @param ... parameters to pass to ComplexHeatmap::Heatmap
#'
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap
#' @importFrom grid gpar
#'
#' @return Heatmap of class ComplexHeatmap
#'
#' @export
plotHeatmap <- function(bank, assay = 'own.expr', dataset = NULL,
                        col = NULL, col.breaks = NULL,
                        cluster.row = TRUE, cluster.column = FALSE,
                        row.dend = TRUE, column.dend = FALSE,
                        cex.row = 4, cex.column = 4,
                        annotate = FALSE,
                        annotate.by = NULL,
                        annotate.col = NULL,
                        name = 'Expression',
                        ...) {

  mat <- getAssay(bank, assay, dataset)
  col.fun <- getHeatmapPalette(mat, col, col.breaks)

  if (annotate) {

    if(is.null(annotate.by) | !(annotate.by %in% names(bank@meta.data))) {
      stop('Specify a valid metadata column to annotate by.')
    }

    ## Order by clusters
    mdata <- bank@meta.data
    mdata <- mdata[order(mdata[[annotate.by]]),]
    mdata <- mdata[mdata$cell_ID %in% colnames(mat), ]
    cell.order <- mdata$cell_ID
    cell.split <- mdata[[annotate.by]]

    ## Col side colours
    clusters <- unique(cell.split)
    n <- length(clusters)
    cluster.cols <- getPalette(n)
    if (!is.null(annotate.col)) {
      if (length(annotate.col) < n) stop('Not enough colors for annotations.')
      else cluster.cols <- annotate.col[seq_len(n)]
    }
    names(cluster.cols) <- clusters

    ha <- HeatmapAnnotation(Clusters = factor(cell.split),
                            col = list(Clusters = cluster.cols),
                            show_annotation_name = FALSE)

    ht <- Heatmap(mat[, cell.order],
                  cluster_rows = cluster.row,
                  cluster_columns = cluster.column,
                  show_row_dend = row.dend,
                  show_column_dend = column.dend,
                  row_names_gp = gpar(fontsize = cex.row),
                  column_names_gp = gpar(fontsize = cex.column),
                  column_split = cell.split,
                  top_annotation = ha,
                  name = name,
                  col = col.fun,
                  ...)

  } else {

    ht <- Heatmap(mat,
                  cluster_rows = cluster.row,
                  cluster_columns = cluster.column,
                  show_row_dend = row.dend,
                  show_column_dend = column.dend,
                  row_names_gp = gpar(fontsize = cex.row),
                  column_names_gp = gpar(fontsize = cex.column),
                  name = name,
                  col = col.fun,
                  ...)

  }
  return(ht)
}

#' Plot Alluvial plot for different cluster assignments
#'
#' @param bank BanksyObject
#'
#' @importFrom data.table data.table melt
#' @importFrom ggplot2 ggplot aes geom_text scale_fill_manual coord_flip guides
#'   theme_minimal xlab
#' @importFrom ggalluvial geom_flow geom_stratum
#'
#' @return Alluvial plot
#'
#' @export
plotAlluvia <- function(bank) {

  mdata <- bank@meta.data
  adata <- data.table(mdata[grep('^res|cell_ID', colnames(mdata))])
  adata <- melt(adata, id.vars = 'cell_ID')
  n <- max(adata$value)
  adata$value <- factor(adata$value)

  variable <- cell_ID <- value <- NULL
  alluvia <- ggplot(adata,
                    aes(x = variable,
                        alluvium = cell_ID,
                        stratum = value,
                        fill = value,
                        label = value)) +
    geom_flow(reverse=FALSE) +
    geom_stratum(alpha = .5, reverse = FALSE) +
    geom_text(stat = "stratum", reverse = FALSE) +
    scale_fill_manual(values = getPalette(n)) +
    coord_flip() +
    guides(fill = FALSE) +
    theme_minimal() +
    xlab('Assignment')

  return(alluvia)
}

