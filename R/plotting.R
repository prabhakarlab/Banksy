# plotUMAP
# plotSpatialDims

#' Plot dimensionality reduction assays
#'
#' @param bank BanksyObject
#' @param reduction reduction to visualize
#' @param components which components to visualize
#' @param by groupings for cells
#' @param type type of groupings (one of 'discrete' or 'continuous')
#' @param pt.size size of points
#' @param pt.alpha transparency of points
#' @param col.midpoint for continuous labels - midpoint of color gradient
#' @param col.low for continuous labels - color gradient low
#' @param col.mid for continuous labels - color gradient mid
#' @param col.high for continuous labels - color gradient high
#' @param col.discrete for discrete labels - colors for groupings
#' @param main title
#' @param main.size size of title
#' @param legend show legend
#' @param legend.text.size size of legend text
#' @param legend.pt.size size of legend point
#'
#' @importFrom ggplot2 ggplot geom_point aes xlab ylab theme_minimal ggtitle
#'   facet_wrap guides guide_legend element_text element_blank scale_color_manual
#'   element_line ensym scale_color_gradient2
#'
#' @return NULL
#'
#' @export
#' 
#' @examples 
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' bank <- NormalizeBanksy(bank)
#' bank <- ScaleBanksy(bank)
#' bank <- ComputeBanksy(bank)
#' bank <- RunPCA(bank, lambda = 0.2)
#' names(reduction(bank))
#' plotReduction(bank, reduction = 'pca_0.2', by = 'Label', type = 'discrete')
#' 
plotReduction <- function(bank, reduction, components = c(1,2),
                          by = NULL, type = c('discrete', 'continuous'),
                          pt.size = 0.5, pt.alpha = 0.7,
                          col.midpoint = NULL, col.low = 'blue',
                          col.mid = 'gray95', col.high = 'red',
                          col.discrete = NULL,
                          main = NULL, main.size = 5,
                          legend = TRUE, legend.text.size = 6,
                          legend.pt.size = 3) {

  data <- getReduction(bank, reduction, components)
  x <- colnames(data)[1]
  y <- colnames(data)[2]

  if (is.null(by)) plot <- ggplot(data, aes(x = !!ensym(x), y = !!ensym(y)))
  else {

    feature <- getFeature(bank, by = by, dataset = NULL)
    checkType(type)
    
    if (type == 'continuous') {
      if (is.null(col.midpoint)) col.midpoint <- median(feature)
      plot <- ggplot(data, aes(x = !!ensym(x), y = !!ensym(y), col = feature)) +
        scale_color_gradient2(midpoint = col.midpoint, low = col.low,
                              mid = col.mid, high = col.high)
    }

    if (type == 'discrete') {
      plot <- ggplot(data, aes(x = !!ensym(x), y = !!ensym(y), 
                               col = as.factor(feature))) +
        scale_color_manual(values = getDiscretePalette(feature, col.discrete)) +
        guides(color = guide_legend(override.aes = list(size = legend.pt.size)))
    }

  }

  legend.pos <- ifelse(legend, 'right', 'none')
  plot <- plot + geom_point(size = pt.size, alpha = pt.alpha) + ggtitle(main) +
    theme_blank(legend.text.size = legend.text.size,
                main.size = main.size,
                legend.pos = legend.pos)

  return(plot)
}

#' Plot Spatial dims
#'
#' @param bank BanksyObject
#' @param dataset if multiple dataset are run
#' @param by groupings for cells
#' @param type type of groupings (one of 'discrete' or 'continuous')
#' @param pt.size size of points
#' @param pt.alpha transparency of points
#' @param col.midpoint for continuous labels - midpoint of color gradient
#' @param col.low for continuous labels - color gradient low
#' @param col.mid for continuous labels - color gradient mid
#' @param col.high for continuous labels - color gradient high
#' @param col.highpoint for continuous labels - high limit of scale
#' @param col.lowpoint for continuous labels - low limit of scale
#' @param na.value for continuous labels - colors to use for missing values
#' @param col.discrete for discrete labels - colors for groupings
#' @param main title
#' @param main.size size of title
#' @param legend show legend
#' @param legend.text.size size of legend text
#' @param legend.pt.size size of legend point
#' @param wrap wrap by feature
#'
#' @importFrom ggplot2 ggplot geom_point aes xlab ylab theme_minimal ggtitle
#'   facet_wrap facet_grid guides guide_legend element_text element_blank 
#'   scale_color_manual element_line scale_color_gradient2
#'
#' @return Spatial plot
#'
#' @export
#' 
#' @examples 
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' plotSpatial(bank, by = 'Label', type = 'discrete')
#' 
plotSpatial <- function(bank, dataset = NULL,
                        by = NA, type = c('discrete', 'continuous'),
                        pt.size = 0.5, pt.alpha = 0.7, col.midpoint = NULL,
                        col.lowpoint = NULL, col.highpoint = NULL,
                        col.low = 'blue', col.mid = 'gray95', col.high = 'red',
                        na.value = 'gray', col.discrete = NULL, main = NULL,
                        main.size = 5, legend = TRUE,
                        legend.text.size = 6, legend.pt.size = 3,
                        wrap = FALSE) {

  data <- getLocations(bank, dataset = dataset)
  sdimx <- sdimy <- NULL


  if (is.na(by)) plot <- ggplot(data, aes(x = sdimx, y = sdimy))
  else {

    feature <- getFeature(bank, by = by, dataset = dataset)
    checkType(type)
    
    if (type == 'continuous') {
      if (is.null(col.midpoint)) col.midpoint <- median(feature)
      if ((is.null(col.highpoint)) | (is.null(col.lowpoint))){
        plot <- ggplot(data, aes(x = sdimx, y = sdimy, col = feature)) +
          scale_color_gradient2(midpoint = col.midpoint, low = col.low,
                                mid = col.mid, high = col.high)
      } else {
        feature[feature > col.highpoint] <- col.highpoint
        feature[feature < col.lowpoint] <- col.lowpoint
        plot <- ggplot(data, aes(x = sdimx, y = sdimy, col = feature)) +
          scale_color_gradient2(midpoint = col.midpoint, 
                                limits = c(col.lowpoint, col.highpoint),
                                low = col.low, mid = col.mid, high = col.high, 
                                na.value = na.value)
      }
    }

    if (type == 'discrete') {
      data <- cbind(data, feature = as.factor(feature))
      plot <- ggplot(data, aes(x = sdimx, y = sdimy, col = as.factor(feature))) +
        scale_color_manual(values = getDiscretePalette(feature, col.discrete)) +
        guides(color = guide_legend(override.aes = list(size = legend.pt.size)))
    }

  }

  legend.pos <- ifelse(legend, 'right', 'none')
  plot <- plot + geom_point(size = pt.size, alpha = pt.alpha) + ggtitle(main) +
    xlab('x coordinates') + ylab('y coordinates') +
    theme_blank(legend.text.size = legend.text.size,
                main.size = main.size,
                legend.pos = legend.pos)
  
  wrap.z <- 'sdimz' %in% names(data)
  if (wrap.z & wrap) plot <- plot + facet_grid(sdimz ~ feature) 
  else if (wrap) plot <- plot + facet_wrap(~ feature)
  else if (wrap.z) plot <- plot + facet_wrap(~ sdimz)
  
  return(plot)
}

#' Wrapper for plotSpatial
#'
#' @param bank BanksyObject
#' @param by features to plot
#' @param type type corresponding to features
#' @param dataset if multiple datasets are used
#' @param nrow nrow for grob
#' @param ncol ncol for grob
#' @param return.plot return plots
#' @param ... parameters to pass to plotSpatial
#'
#' @importFrom gridExtra grid.arrange
#'
#' @return Spatial plots
#' 
#' @export
#' 
#' @examples 
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' bank <- NormalizeBanksy(bank)
#' bank <- ScaleBanksy(bank)
#' bank <- ComputeBanksy(bank)
#' bank <- RunPCA(bank, lambda = 0.2)
#' bank <- ClusterBanksy(bank, lambda = 0.2, npcs = 20, k.neighbors = 50, resolution = 0.8)
#' plotSpatialFeatures(bank, by = c('Label', clust.names(bank)), 
#'   type = rep('discrete', 2), nrow = 1, ncol = 2, pt.size = 2)
#' 
plotSpatialFeatures <- function(bank, by, type, dataset = NULL,
                                nrow = NULL, ncol= NULL, return.plot = FALSE, 
                                ...) {

  valid <- length(by) == length(type)
  if (!valid) stop('Ensure a 1-1 correspondence between features to plot (by)
                   and type of features (type)')

  plots <- Map(f = function(feature, feature.type, ...) {
    plotSpatial(bank, dataset = dataset, by = feature, type = feature.type, ...)
  }, by, type, ...)

  do.call(grid.arrange, c(plots, nrow = nrow, ncol = ncol))
  
  if (return.plot) return(plots)
}

#' Plot Heatmap (wrapper for ComplexHeatmap)
#'
#' @param bank BanksyObject
#' @param assay assay to plot heatmap (one of own.expr, nbr.expr or banksy)
#' @param dataset dataset to plot heatmap
#' @param lambda lambda if assay == banksy
#' @param cells specific cells to plot
#' @param features specific features to plot
#' @param col colours to use in heatmap
#' @param col.breaks color breaks to use in heatmap (same number as col is
#'   specified)
#' @param col.discrete cluster colors, named color array, names are clusters,
#'   values are colors in hexadecimal
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
#' @param max.cols max columns to display - will subsample if not NULL
#' @param seed used for sampling
#' @param rasterize rasterize if TRUE
#' @param ... parameters to pass to ComplexHeatmap::Heatmap
#'
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap
#' @importFrom grid gpar
#'
#' @return Heatmap of class ComplexHeatmap
#'
#' @export
#' 
#' @examples 
#' # Generate a simulated dataset
#' d <- simulateDataset(n_cells = c(20,20,20), n_genes = 30)
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' bank <- NormalizeBanksy(bank)
#' bank <- ScaleBanksy(bank)
#' plotHeatmap(bank, annotate = TRUE, annotate.by = 'Label')
#' 
plotHeatmap <- function(bank, assay = 'own.expr',
                        dataset = NULL, lambda = NULL,
                        cells = NULL, features = NULL,
                        col = NULL, col.breaks = NULL,
                        col.discrete = NULL,
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
                        max.cols = NULL,
                        seed = 42,
                        rasterize = FALSE,
                        ...) {

  mat <- getAssay(bank, assay, dataset, lambda, cells, features)
  if (!is.null(max.cols)) mat <- sampleMatrix(mat, max.cols, seed)
  col.fun <- getHeatmapPalette(mat, col, col.breaks)

  if (annotate) {

    anno <- getCellAnnotation(bank = bank, mat = mat,
                              annotate.by = annotate.by,
                              annotation.name = annotation.name,
                              annotation.size = annotation.size,
                              annotation.pos = annotation.pos,
                              order.by = order.by, col.discrete = col.discrete)
    ha <- anno$anno
    group <- NULL
    ra <- NULL

    if (assay == 'banksy') {
      group <- rep('own', nrow(mat))
      group[grepl('.nbr$', rownames(mat))] <- 'nbr'
      group <- factor(group, levels = c('own', 'nbr'))
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
                  use_raster = rasterize,
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
                  use_raster = rasterize,
                  ...)

  }

  return(ht)
}

#' Plot ARI matrix
#'
#' @param bank BanksyObject with harmonised clusterings
#' @param col.low color gradient low
#' @param col.high color gradient high
#' @param label TRUE if to label plot with ARI
#' @param digits number of digits to round ARI to
#'
#' @importFrom ggplot2 ggplot geom_tile scale_fill_gradient theme_minimal theme
#'   element_blank labs geom_text element_text
#'
#' @export
#' 
#' @examples 
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' bank <- NormalizeBanksy(bank)
#' bank <- ScaleBanksy(bank)
#' bank <- ComputeBanksy(bank)
#' bank <- RunPCA(bank, lambda = 0.2)
#' bank <- ClusterBanksy(bank, lambda = 0.2, npcs = 20, k.neighbors = 50, resolution = c(0.5,1.5))
#' plotARI(bank)
#' 
plotARI <- function(bank, col.low = 'white', col.high = 'red', label = TRUE, digits = 3) {
  ari.mat <- getARI(bank, digits = digits)
  data <- reshape2::melt(ari.mat)
  data$Var1 <- factor(data$Var1, levels = colnames(ari.mat))
  data$Var2 <- factor(data$Var2, levels = rev(colnames(ari.mat)))
  Var1 <- Var2 <- value <- NULL
  tile <- ggplot(data, aes(x=Var1, y=Var2, fill = value)) + geom_tile()+
    scale_fill_gradient(low = 'white', high = 'red') + theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 90)) +
    labs(fill = 'ARI')
  if (label) tile <- tile + geom_text(aes(label = round(value, digits)))
  return(tile)
}


#' Plot Alluvial plot for different cluster assignments
#'
#' @param bank BanksyObject
#' @param max.cells Number of cells to display
#' @param seed Seed for sampling cells if max.cells is not NULL
#' @param flow.colors Colors for alluvial flow
#'
#' @importFrom ggplot2 ggplot aes theme theme_minimal element_blank
#'   element_line element_text scale_x_discrete scale_fill_manual
#' @importFrom ggalluvial geom_flow geom_stratum StatStratum
#'
#' @return Alluvial plot
#'
#' @export
#' 
#' @examples 
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' bank <- NormalizeBanksy(bank)
#' bank <- ScaleBanksy(bank)
#' bank <- ComputeBanksy(bank)
#' bank <- RunPCA(bank, lambda = 0.2)
#' bank <- ClusterBanksy(bank, lambda = 0.2, npcs = 20, k.neighbors = 50, resolution = c(0.5,1.5))
#' plotAlluvia(bank, max.cells = 20)
#' 
plotAlluvia <- function(bank, max.cells = 500, seed = 42, flow.colors = NULL) {

  df <- meta.data(bank)[,clust.names(bank),drop=FALSE]

  if (ncol(df) < 2) {
    stop('Alluvia will only be plotted for at least 2 clustering runs.')
  }

  if (!is.null(max.cells)) {
    set.seed(seed)
    df <- df[sample(1:nrow(df), max.cells), ]
  }

  cell <- rep(1:nrow(df), times = ncol(df))
  run <- rep(names(df), each = nrow(df))
  cluster <- factor(unlist(df))
  plotdf <- data.frame(cell = cell, run = run, cluster = cluster)
  nclust <- length(levels(cluster))

  if (!is.null(flow.colors)) {
    if (length(flow.colors) < nclust) {
      stop('Not enough colors. Need ', nclust,
           ' but supplied ', length(flow.colors))
    } else {
      flow.colors <- flow.colors[1:nclust]
    }
  } else {
    flow.colors <- getPalette(nclust)
  }

  StatStratum <- ggalluvial::StatStratum

  alluvia <- ggplot(plotdf,
                    aes(x = run, stratum = cluster, alluvium = cell,
                        fill = cluster, label = cluster)) +
    geom_flow(stat = 'alluvium', lode.guidance = 'forward',
              color = 'transparent', alpha = 0.4) +
    geom_stratum(alpha = 0.75) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = 'grey40'),
          axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank()) +
    scale_x_discrete(expand = c(0.05, 0.05)) +
    scale_fill_manual(values = flow.colors)

  return(alluvia)
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
        dataset <- names(bank@own.expr)[1]
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
  mat <- mat[rownames(mat) %in% features, colnames(mat) %in% cells, drop = FALSE]
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
    names(cluster.cols) <- sort(clusters)
  } else if (is.character(x)) {
    n <- length(clusters)
    cluster.cols <- getPalette(n)
    names(cluster.cols) <- clusters
  }


  return(cluster.cols)
}

#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @importFrom grid gpar
getCellAnnotation <- function(bank, mat,
                              annotate.by,
                              annotation.name,
                              annotation.size,
                              annotation.pos,
                              col.discrete = NULL,
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
  #hc <- lapply(mdata, getDiscretePalette)
  hc <- lapply(mdata, FUN = getDiscretePalette, col.discrete = col.discrete)
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

#' @importFrom ggplot2 theme element_text element_blank element_line
theme_blank <- function(legend.text.size,
                        main.size,
                        legend.pos) {
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size = main.size),
        legend.text = element_text(size = legend.text.size),
        legend.title = element_blank(),
        legend.position = legend.pos,
        legend.key = element_blank())
}

getFeature <- function(bank, by, dataset) {

  if (is.null(by)) return(NULL)
  if (length(by) > 1) {
      warning('More than one variable to plot by supplied. Using first.')
      by <- by[1]
  }
  if (is.null(dataset)) {
    found.meta <- by %in% names(bank@meta.data)
    found.own <- by %in% rownames(bank@own.expr)
    found.nbr <- by %in% rownames(bank@nbr.expr)
    found <- found.meta | found.own | found.nbr
    if (!found) stop('Invalid parameter to plot by.')
    if (found.meta) return(bank@meta.data[[by]])
    if (found.own) return(bank@own.expr[by,])
    if (found.nbr) return(bank@nbr.expr[by,])

  } else {
    found.meta <- by %in% names(bank@meta.data)
    found.own <- by %in% rownames(bank@own.expr[[1]])
    found.nbr <- by%in% rownames(bank@nbr.expr[[1]])
    found <- found.meta | found.own | found.nbr
    if (!found) stop('Invalid parameter to plot by.')
    if (found.meta) return(bank@meta.data[[by]][bank@meta.data[['dataset']] ==
                                                dataset])
    if (found.own) return(bank@own.expr[[dataset]][by,])
    if (found.nbr) return(bank@nbr.expr[[dataset]][by,])
  }
}

getLocations <- function(bank, dataset) {

  if (is.list(bank@own.expr)) {  # Multiple dataset case
    if (is.null(dataset)) stop('No dataset specified.')
    data <- bank@cell.locs[[dataset]]
  } else {  # Single dataset case
    data <- bank@cell.locs
  }
  return(data)

}

getReduction <- function(bank, reduction, components) {

  found <- reduction %in% names(bank@reduction)
  if (!found) stop('No reduction named ', reduction)
  valid <- length(components) == 2
  if (!valid) stop('Only 2 components may be visualized, not ', length(components))
  data <- NULL
  if (grepl('pca', reduction)) data <- bank@reduction[[reduction]]$x[,components]
  if (grepl('umap', reduction)) data <- bank@reduction[[reduction]][,components]
  if (is.null(data)) stop('Specify a valid reduction')
  return(data.frame(data))

}

checkType <- function(type) {

  message <- 'Specify type - must be either discrete or continuous'
  if (length(type) > 1) stop(message)
  if (!(type %in% c('discrete', 'continuous'))) stop(message)

}

getDiscretePalette <- function(feature, col.discrete = NULL) {

  num <- is.numeric(feature)
  char <- is.character(feature)

  if (num) n <- max(feature)
  if (char) n <- length(unique(feature))

  if (is.null(col.discrete)) {
    pal <- getPalette(n)
    if (num) {
      pal <- pal[sort(unique(feature))]
      names(pal) <- sort(unique(feature))
    } else if (char) {
      names(pal) <- unique(feature)
    }
  } else {
    if (length(col.discrete) < n) stop('Not enough colors')
    pal <- col.discrete
    if (is.null(names(col.discrete))) stop('
       Color palatte must be supplied with valid cluster names.')
    if (!all(unique(feature) %in% names(col.discrete))) stop(
      'Not all clusters have an assigned color.
       Ensure that names(col.discrete) has all the cluster labels.')

    if (num) {
      clust.as.char <- as.character(unique(feature))
      pal <- col.discrete[clust.as.char]
    } else if (char) {
      clust.as.char <- unique(feature)
      pal <- col.discrete[clust.as.char]
    }
  }

  return(pal)

}

sampleMatrix <- function(mat, max.cols, seed) {
  set.seed(seed)
  message(paste0('Sampling ', max.cols, ' columns with seed ', seed))
  message(paste0('Keeping ', round(max.cols/ncol(mat) * 100,2), '% of cells'))
  keep <- sample(colnames(mat), min(length(colnames(mat)), max.cols))
  mat <- mat[, keep, drop=FALSE]
  return(mat)
}
