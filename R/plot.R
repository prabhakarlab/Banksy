# plotUMAP
# plotSpatialDims

#' Plot UMAP
#'
#' @param bank BanksyObject
#' @param by groupings for cells
#' @param reduction reduction to use
#' @param col cluster colours
#' @param legend cluster legend
#' @param main title
#' @param pt.size size of points
#' @param main.size size of title
#' @param legend.text.size size of legend text
#' @param legend.pt.size size of legent point
#'
#' @importFrom ggplot2 ggplot geom_point aes xlab ylab theme_minimal ggtitle
#'   scale_color_manual theme element_text element_blank guides guide_legend
#' @importFrom plyr mapvalues
#'
#' @return NULL
#'
#' @export
plotUMAP <- function(bank, by, reduction,
                     col = NULL, legend = TRUE,
                     main = NULL,
                     pt.size = 1,
                     main.size = 5,
                     legend.text.size = 6,
                     legend.pt.size = 3) {

  mnames <- names(bank@meta.data)
  mnames <- mnames[!grepl('cell_ID|n_features', mnames)]
  if (length(by) > 1) {
    warning('More than one plot parameter supplied. Using first.')
    by <- by[1]
  }
  if (!(by %in% names(bank@meta.data))) {
    stop(paste0('Invalid parameter to plot UMAP by. One of ',
                paste(mnames, collapse = ' ')))
  }

  umaps <- names(bank@dim.reduction)[grep('umap', names(bank@dim.reduction))]
  if (length(reduction) > 1) {
    warning('More than one reduction supplied. Using first.')
    reduction <- reduction[1]
  }
  if (!(reduction %in% names(bank@dim.reduction))) {
    stop(paste0('Invalid umap name. One of ',
                paste(umaps, collapse = ' ')))
  }

  ## Get plotting colours
  clusters <- bank@meta.data[[by]]
  if (is.numeric(clusters)) n <- max(clusters)
  if (is.character(clusters)) n <- length(unique(clusters))
  if (is.null(col)) {
    plotCols <- getPalette(n)
    if (is.numeric(clusters)) plotCols <- plotCols[sort(unique(clusters))]
  } else {
    if (length(col) < n) stop('Not enough colours')
    plotCols <- col[seq_len(n)]
  }

  ## Initialize plot
  umap_1 <- umap_2 <- Clust <- NULL
  plotData <- as.data.frame(bank@dim.reduction[[reduction]])
  plotData <- cbind(plotData, Clust = as.factor(clusters))
  names(plotData)[seq_len(2)] <- c('umap_1', 'umap_2')

  p <- ggplot(plotData, aes(x=umap_1, y=umap_2, col=Clust)) +
    xlab('umap 1') + ylab('umap 2') + theme_minimal() +
    geom_point(size=pt.size) +
    scale_color_manual(values=plotCols) +
    theme(legend.text = element_text(size = legend.text.size),
          legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = legend.pt.size)))


  if (!legend) p <- p + theme(legend.position = 'none')
  if (is.null(main)) p <- p + ggtitle(by) +
    theme(plot.title = element_text(size = main.size))

  return(p)
}

#' Plot Spatial dims
#'
#' @param bank BanksyObject
#' @param by clustering run to plot by
#' @param dataset if multiple dataset are run
#' @param col colours
#' @param legend show legend
#' @param main title
#' @param pt.size size of points
#' @param main.size size of title
#' @param legend.text.size size of legend text
#' @param legend.pt.size size of legend point
#'
#' @importFrom ggplot2 ggplot geom_point aes xlab ylab theme_minimal ggtitle
#'   facet_wrap guides guide_legend element_text element_blank scale_color_manual
#' @importFrom plyr mapvalues
#'
#' @return NULL
#'
#' @export
plotSpatialDims <- function(bank, by, dataset = NULL,
                            col = NULL, legend = TRUE,
                            main = NULL,
                            pt.size = 1,
                            main.size = 5,
                            legend.text.size = 6,
                            legend.pt.size = 3) {

  mnames <- names(bank@meta.data)
  mnames <- mnames[!grepl('cell_ID|n_features', mnames)]
  if (length(by) > 1) {
    warning('More than one plot parameter supplied. Using first.')
    by <- by[1]
  }
  if (!(by %in% names(bank@meta.data))) {
    stop(paste0('Invalid parameter to plot UMAP by. One of ',
                paste(mnames, collapse = ' ')))
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

  meta <- bank@meta.data[bank@meta.data$cell_ID %in% rownames(plotData), ]
  clusters <- bank@meta.data[[by]][bank@meta.data$cell_ID %in%
                                     rownames(plotData)]
  if (is.numeric(clusters)) n <- max(clusters)
  if (is.character(clusters)) n <- length(unique(clusters))
  if (is.null(col)) {
    plotCols <- getPalette(n)
    if (is.numeric(clusters)) plotCols <- plotCols[sort(unique(clusters))]
  } else {
    if (length(col) < n) stop('Not enough colours')
    plotCols <- col[seq_len(n)]
  }

  dimx <- dimy <- Clust <- NULL
  names(plotData)[seq_len(2)] <- c('dimx','dimy')
  plotData <- cbind(plotData, Clust = as.factor(clusters), meta)
  p <- ggplot(plotData, aes(x=dimx, y=dimy, col=Clust)) +
    xlab('x coordinates') + ylab('y coordinates') + theme_minimal() +
    geom_point(size=pt.size) +
    scale_color_manual(values=plotCols) +
    theme(legend.text = element_text(size = legend.text.size),
          legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = legend.pt.size)))


  if ('sdimz' %in% names(plotData)) {
    plotData$sdimz <- as.factor(plotData$sdimz)
    p <- p + facet_wrap( ~ sdimz)
  }

  if (!legend) p <- p + theme(legend.position = 'none')
  if (is.null(main)) p <- p + ggtitle(by) +
    theme(plot.title = element_text(size = main.size))

  return(p)
}

#' Plot Heatmap (wrapper for ComplexHeatmap)
#'
#' @param bank BanksyObject
#' @param assay assay to plot heatmap
#' @param dataset dataset to plot heatmap
#' @param lambda lambda if assay == banksy
#' @param col colours to use in heatmap
#' @param col.breaks color breaks to use in heatmap (same number as col is
#'   specified)
#' @param cluster.row cluster rows
#' @param cluster.column cluster columns
#' @param row.dend draw row dendrograms
#' @param column.dend draw column dendrograms
#' @param cex.row row label size
#' @param annotate add annotation for cells
#' @param annotate.by metadata to annotate cells by
#' @param order.by metadata to order cells by (one of annotate.by)
#' @param annotation.name show annotation name
#' @param annotation.size size of annotation labels
#' @param annotation.pos position of annotation labels
#' @param name name of heatmap legend
#' @param barplot.by metadata to plot barplots (numeric)
#' @param barplot.border show borders of barplot
#' @param barplot.width barplot width
#' @param ... parameters to pass to ComplexHeatmap::Heatmap
#'
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap
#' @importFrom grid gpar
#'
#' @return Heatmap of class ComplexHeatmap
#'
#' @export
plotHeatmap <- function(bank, assay = 'own.expr',
                        dataset = NULL, lambda = NULL,
                        col = NULL, col.breaks = NULL,
                        cluster.row = TRUE, cluster.column = FALSE,
                        row.dend = FALSE, column.dend = FALSE,
                        cex.row = 4,
                        annotate = FALSE,
                        annotate.by = NULL,
                        order.by = NULL,
                        annotation.name = TRUE,
                        annotation.size = 6,
                        annotation.pos = 'right',
                        name = 'Expression',
                        barplot.by = NULL,
                        barplot.border = FALSE,
                        barplot.width = 0.6,
                        ...) {

  mat <- getAssay(bank, assay, dataset, lambda)
  col.fun <- getHeatmapPalette(mat, col, col.breaks)

  if (annotate) {

    anno <- getCellAnnotation(bank = bank, mat = mat,
                              annotate.by = annotate.by,
                              annotation.name = annotation.name,
                              annotation.size = annotation.size,
                              annotation.pos = annotation.pos,
                              order.by = order.by)
    ha <- anno$anno
    group <- NULL
    ra <- NULL

    if (assay == 'banksy') {
      group <- rep('own', nrow(mat))
      group[grepl('.nbr$', rownames(mat))] <- 'nbr'
      ra <- HeatmapAnnotation(Clusters = factor(group),
                              which = 'row',
                              col = list(Clusters=c('own'='transparent',
                                                    'nbr'='transparent')),
                              show_annotation_name = FALSE,
                              show_legend = FALSE)

    }

    ht <- Heatmap(mat[, anno$cell.order],
                  cluster_rows = cluster.row,
                  cluster_columns = cluster.column,
                  show_row_dend = row.dend,
                  show_column_dend = column.dend,
                  row_names_gp = gpar(fontsize = cex.row),
                  show_column_names = FALSE,
                  column_title = NULL,
                  column_split = anno$cell.split,
                  row_split = group,
                  top_annotation = ha,
                  left_annotation = ra,
                  name = name,
                  col = col.fun,
                  ...)

    if (!is.null(barplot.by)) {

      ht <- appendBarplots(bank = bank, mat = mat,
                           cell.order = anno$cell.order,
                           annotation.name = annotation.name,
                           annotation.size = annotation.size,
                           annotation.pos = annotation.pos,
                           barplot.by = barplot.by,
                           barplot.border = barplot.border,
                           barplot.width = barplot.width,
                           heatmap = ht)
    }

  } else {
    ht <- Heatmap(mat,
                  cluster_rows = cluster.row,
                  cluster_columns = cluster.column,
                  show_row_dend = row.dend,
                  show_column_dend = column.dend,
                  row_names_gp = gpar(fontsize = cex.row),
                  show_column_names = FALSE,
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

