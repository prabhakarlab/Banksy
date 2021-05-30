# plotMetaDataHeatmap_vs
# plotHeatmap_banksy
# plotSpatDimPlots


#' Plots expression heatmaps
#'
#' @param gobject Giotto object
#' @param expression_values one of normalized, scaled, custom
#' @param metadata_cols mdata
#' @param selected_genes genes to use
#' @param first_meta_col 1st meta col
#' @param second_meta_col 2nd meta col
#' @param show_values one of zscores, rescaled, original, zscores_rescaled
#' @param custom_cluster_order cell clustering order
#' @param clus_cor_method cell correlation method
#' @param clus_cluster_method cell clustering method
#' @param custom_gene_order gene order
#' @param gene_cor_method gene correlation method
#' @param gene_cluster_method gene clustering method
#' @param gradient_color grad colour
#' @param gradient_midpoint grad midpoint
#' @param gradient_limits gradlimits
#' @param x_text_size x text size
#' @param x_text_angle x text angle
#' @param y_text_size y text size
#' @param strip_text_size strip text size
#' @param show_plot show plot
#' @param return_plot return plot
#' @param save_plot save plot
#' @param save_param save parameters
#' @param default_save_name default saving name
#'
#' @importFrom data.table `:=`
#' @importFrom ggplot2 ggplot geom_tile aes_string scale_fill_gradient2
#'   theme_classic theme element_text element_blank labs facet_grid
#' @importFrom Giotto calculateMetaTable readGiottoInstructions
#' @importFrom scales rescale
#' @importFrom stats as.dist cor hclust reformulate
#'
#' @return Heatmap
#'
#' @export
plotMetaDataHeatmap_vs <- function (gobject, expression_values = c("normalized", "scaled", "custom"),
                                    metadata_cols = NULL, selected_genes = NULL,
                                    first_meta_col = NULL, second_meta_col = NULL,
                                    show_values = c("zscores", "rescaled", "original", "zscores_rescaled"),
                                    custom_cluster_order = NULL, clus_cor_method = "pearson", clus_cluster_method = "complete",
                                    custom_gene_order = NULL, gene_cor_method = "pearson", gene_cluster_method = "complete",
                                    gradient_color = c("blue", "white", "red"), gradient_midpoint = 0.5, gradient_limits = NULL,
                                    x_text_size = 10, x_text_angle = 45,
                                    y_text_size = 10, strip_text_size = 8,
                                    show_plot = NA, return_plot = NA,
                                    save_plot = NA, save_param = list(), default_save_name = "plotMetaDataHeatmap")
{
  metaDT = calculateMetaTable(gobject = gobject, expression_values = expression_values,
                              metadata_cols = metadata_cols, selected_genes = selected_genes)
  # Global binding for variables
  zscores <- value <- zscores_rescaled_per_gene <- value_rescaled_per_gene <-
    factor_column <- variable <- factor_1_column <- factor_2_column <- NULL
  metaDT[, `:=`(zscores, scale(value)), by = c("variable")]
  metaDT[, `:=`(zscores_rescaled_per_gene, rescale(zscores, to = c(-1, 1))), by = c("variable")]
  metaDT[, `:=`(value_rescaled_per_gene, rescale(value,to = c(0, 1))), by = c("variable")]
  show_values = match.arg(show_values, choices = c("zscores", "rescaled", "original", "zscores_rescaled"))
  if (show_values == "zscores") {
    show_values = "zscores"
  }
  else if (show_values == "original") {
    show_values = "value"
  }
  else if (show_values == "rescaled") {
    show_values = "value_rescaled_per_gene"
    gradient_midpoint = 0.5
  }
  else {
    show_values = "zscores_rescaled_per_gene"
  }
  if (length(metadata_cols) > 2) {
    cat("\n visualization is only possible for 1 or 2 metadata annotations, data.table is returned \n")
    return(metaDT)
  }
  main_factor = ifelse(length(metadata_cols) == 1, metadata_cols, first_meta_col)
  testmain = metaDT[, mean(value), by = c("variable", main_factor)]
  testmain_matrix = dfunction(d = testmain, col_name1 = "variable",
                              col_name2 = main_factor, value.var = "V1")
  testmain_mat = as.matrix(testmain_matrix[, -1])
  rownames(testmain_mat) = testmain_matrix$variable
  if (is.null(custom_cluster_order)) {
    cormatrix = cor(x = testmain_mat, method = clus_cor_method)
    cordist = as.dist(1 - cormatrix, diag = T, upper = T)
    corclus = hclust(d = cordist, method = clus_cluster_method)
    clus_names = rownames(cormatrix)
    names(clus_names) = 1:length(clus_names)
    clus_sort_names = clus_names[corclus$order]
  }
  else {
    clus_sort_names = unique(as.character(custom_cluster_order))
    if (all(colnames(testmain_mat) %in% clus_sort_names) ==
        FALSE) {
      stop("\n custom cluster order is given, but not all clusters are represented \n")
    }
  }
  if (is.null(custom_gene_order)) {
    gene_cormatrix = cor(x = t(testmain_mat), method = gene_cor_method)
    gene_cordist = as.dist(1 - gene_cormatrix, diag = T,
                                  upper = T)
    gene_corclus = hclust(d = gene_cordist, method = gene_cluster_method)
    gene_names = rownames(gene_cormatrix)
    names(gene_names) = 1:length(gene_names)
    gene_sort_names = gene_names[gene_corclus$order]
  }
  else {
    gene_sort_names = unique(as.character(custom_gene_order))
    if (all(rownames(testmain_mat) %in% gene_sort_names) ==
        FALSE) {
      stop("\n custom gene order is given, but not all genes are represented \n")
    }
  }
  if (length(metadata_cols) == 1) {
    metaDT[, `:=`(factor_column, factor(get(metadata_cols),
                                        levels = clus_sort_names))]
    metaDT[, `:=`(variable, factor(get("variable"), levels = gene_sort_names))]
    if (!is.null(gradient_limits) & is.vector(gradient_limits) &
        length(gradient_limits) == 2) {
      lower_lim = gradient_limits[[1]]
      upper_lim = gradient_limits[[2]]
      numeric_data = metaDT[[show_values]]
      limit_numeric_data = ifelse(numeric_data > upper_lim,
                                  upper_lim, ifelse(numeric_data < lower_lim,
                                                    lower_lim, numeric_data))
      metaDT[[show_values]] = limit_numeric_data
    }
    pl <- ggplot()
    pl <- pl + geom_tile(data = metaDT, aes_string(x = "factor_column",
                                                   y = "variable", fill = show_values), color = "black")
    pl <- pl + scale_fill_gradient2(low = gradient_color[[1]],
                                    mid = gradient_color[[2]], high = gradient_color[[3]],
                                    midpoint = gradient_midpoint)
    pl <- pl + theme_classic()
    pl <- pl + theme(axis.text.x = element_text(size = x_text_size,
                                                angle = x_text_angle, hjust = 1, vjust = 1),
                     axis.text.y = element_text(size = y_text_size),
                     legend.title = element_blank())
    pl <- pl + labs(x = metadata_cols, y = "genes")
    pl
    show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject,
                                                                param = "show_plot"), show_plot)
    save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject,
                                                                param = "save_plot"), save_plot)
    return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject,
                                                                    param = "return_plot"), return_plot)
    if (show_plot == TRUE) {
      print(pl)
    }
    if (save_plot == TRUE) {
      do.call("all_plots_save_function", c(list(gobject = gobject,
                                                plot_object = pl, default_save_name = default_save_name),
                                           save_param))
    }
    if (return_plot == TRUE) {
      return(pl)
    }
  }
  else {
    if (is.null(first_meta_col) | is.null(second_meta_col)) {
      cat("\n both first_meta_col and second_meta_col need to be defined, return data.table \n")
      return(metaDT)
    }
    else {
      metaDT[, `:=`(factor_1_column, factor(get(first_meta_col),
                                            clus_sort_names))]
      metaDT[, `:=`(factor_2_column, as.factor(get(second_meta_col)))]
      metaDT[, `:=`(variable, factor(get("variable"),
                                     levels = gene_sort_names))]
      if (!is.null(gradient_limits) & is.vector(gradient_limits) &
          length(gradient_limits) == 2) {
        lower_lim = gradient_limits[[1]]
        upper_lim = gradient_limits[[2]]
        numeric_data = metaDT[[show_values]]
        limit_numeric_data = ifelse(numeric_data > upper_lim,
                                    upper_lim, ifelse(numeric_data < lower_lim,
                                                      lower_lim, numeric_data))
        metaDT[[show_values]] = limit_numeric_data
      }
      pl <- ggplot()
      pl <- pl + geom_tile(data = metaDT, aes_string(x = "factor_1_column",
                                                     y = "variable", fill = show_values), color = "black")
      pl <- pl + scale_fill_gradient2(low = gradient_color[[1]],
                                      mid = gradient_color[[2]], high = gradient_color[[3]],
                                      midpoint = gradient_midpoint)
      pl <- pl + facet_grid(reformulate("factor_2_column"))
      pl <- pl + theme_classic()
      pl <- pl + theme(axis.text.x = element_text(size = x_text_size,
                                                  angle = x_text_angle, hjust = 1, vjust = 1),
                       axis.text.y = element_text(size = y_text_size),
                       strip.text = element_text(size = strip_text_size),
                       legend.title = element_blank())
      pl <- pl + labs(x = first_meta_col, y = "genes", title = second_meta_col)
      show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject,
                                                                  param = "show_plot"), show_plot)
      save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject,
                                                                  param = "save_plot"), save_plot)
      return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject,
                                                                      param = "return_plot"), return_plot)
      if (show_plot == TRUE) {
        print(pl)
      }
      if (save_plot == TRUE) {
        do.call("all_plots_save_function", c(list(gobject = gobject,
                                                  plot_object = pl, default_save_name = default_save_name),
                                             save_param))
      }
      if (return_plot == TRUE) {
        return(pl)
      }
    }
  }
}

#' Plot Banksy Heatmap
#'
#' @param gobject Giotto object
#' @param expression_values one of normalized, scaled, custom
#' @param plotValuesRange plot value range
#' @param genes genes to use
#' @param cluster_column column to cluster by
#' @param cluster_order custom order for clusters
#' @param cluster_cor_method method for cluster correlation
#' @param cluster_hclust_method method for hierarchical clustering of clusters
#' @param cluster_custom_order custom cluster order
#' @param cluster_color_code custom colors for clusters
#' @param gene_order method to determine gene order
#' @param gene_custom_order custom order for genes
#' @param gene_cor_method method for gene correlation
#' @param gene_hclust_method method for hierarchical clustering of genes
#' @param show_values show values or not
#' @param size_vertical_lines size of vertical lines
#' @param gradient_colors grad colours
#' @param gene_label_selection select genes
#' @param axis_text_y_size size of y axis text
#' @param legend_nrows nrows
#' @param show_plot show plot or not
#' @param return_plot return plot or not
#' @param save_plot save plot or not
#' @param save_param save parameters
#' @param default_save_name save name
#'
#' @importFrom cowplot align_plots get_legend plot_grid
#' @importFrom data.table setorder
#' @importFrom ggplot2 aes aes_string ggplot geom_raster geom_tile geom_vline
#'   guide_legend guide_colorbar theme element_blank margin labs .pt
#'   scale_fill_manual scale_fill_gradient2 scale_x_continuous scale_y_continuous
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid unit
#' @importFrom Giotto getDistinctColors readGiottoInstructions
#'
#' @return Combined heatmap plot
#'
#' @export
plotHeatmap_banksy <- function(gobject, expression_values = c("normalized", "scaled", "custom"),
                               plotValuesRange = NULL,
                               genes, cluster_column = NULL,
                               cluster_order = c("size","correlation", "custom"),
                               cluster_custom_order = NULL, cluster_color_code = NULL,
                               cluster_cor_method = "pearson", cluster_hclust_method = "ward.D",
                               gene_order = c("correlation", "custom"), gene_custom_order = NULL,
                               gene_cor_method = "pearson", gene_hclust_method = "complete",
                               show_values = c("rescaled", "z-scaled", "original"), size_vertical_lines = 1.1,
                               gradient_colors = c("blue", "yellow", "red"), gene_label_selection = NULL,
                               axis_text_y_size = NULL, legend_nrows = 1, show_plot = NA,
                               return_plot = NA, save_plot = NA, save_param = list(), default_save_name = "plotHeatmap") {
    show_values = match.arg(show_values, choices = c("rescaled",
                                                   "z-scaled", "original"))
    heatmap_data = createHeatmap_DT(gobject = gobject, expression_values = expression_values,
                                           genes = genes, cluster_column = cluster_column, cluster_order = cluster_order,
                                           cluster_custom_order = cluster_custom_order, cluster_cor_method = cluster_cor_method,
                                           cluster_hclust_method = cluster_hclust_method, gene_order = gene_order,
                                           gene_custom_order = gene_custom_order, gene_cor_method = gene_cor_method,
                                           gene_hclust_method = gene_hclust_method)
    cell_order_DT = heatmap_data[["cell_DT"]]
    subset_values_DT = heatmap_data[["DT"]]
    x_lines = heatmap_data[["x_lines"]]

    if (is.null(cluster_color_code)) {
        clus_values = unique(cell_order_DT[[cluster_column]])
        clus_colors = getDistinctColors(n = length(clus_values))
        names(clus_colors) = clus_values
    } else {
        clus_colors = cluster_color_code
    }

    clus_pl <- ggplot() + geom_raster(data = cell_order_DT, aes_string(x = "cells", y = "1", fill = cluster_column)) +
        geom_vline(xintercept = x_lines, color = "white", size = size_vertical_lines) +
        scale_fill_manual(values = clus_colors, guide = guide_legend(title = "", nrow = legend_nrows)) +
        theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
            legend.position = "top", plot.margin = margin(0, 0, 0, 5.5, "pt")) +
        labs(x = "", y = "clusters")

    if (!(is.null(plotValuesRange))){
        if (length(plotValuesRange)==2){
            midpoint = mean(plotValuesRange)
        } else {
            stop("plotValueRange should only have two numeric elements.")
        }

        if (show_values == "original") {
          value_column = "expression"
        } else if (show_values == "z-scaled") {
          value_column = "z_scores"
        } else {
          value_column = "scale_scores"
        }
    } else {
        if (show_values == "original") {
            value_column = "expression"
            midpoint = max(subset_values_DT$expression)/2
        } else if (show_values == "z-scaled") {
            value_column = "z_scores"
            midpoint = max(subset_values_DT$z_scores)/2
        } else {
            value_column = "scale_scores"
            midpoint = 0.5
        }
    }

    empty = ggplot() + theme_classic() +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), plot.margin = margin(0, 0, 0, 0, "cm"))

    if (length(gradient_colors) > 3)
        cat("\n only first 3 colors will be used for gradient \n")
        low_color = gradient_colors[[1]]
        mid_color = gradient_colors[[2]]
        high_color = gradient_colors[[3]]

    hmap <- ggplot() + geom_tile(data = subset_values_DT, aes_string(x = "cells", y = "genes", fill = value_column)) +
        geom_vline(xintercept = x_lines, color = "white", size = size_vertical_lines) +
        scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color, midpoint = midpoint, guide = guide_colorbar(title = ""))

    if (is.null(gene_label_selection)) {

        hmap <- hmap + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                             axis.text.y = element_text(size = axis_text_y_size), plot.margin = margin(0, 0, 5.5, 5.5, "pt"))

        aligned <- align_plots(clus_pl, empty, hmap + theme(legend.position = "none"), align = "v", axis = "l")
        aligned <- append(aligned, list(get_legend(hmap)))
        combplot = plot_grid(plotlist = aligned, ncol = 2, rel_widths = c(1, 0.2), nrow = 2, rel_heights = c(0.2,1))

        show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = "show_plot"), show_plot)
        save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = "save_plot"), save_plot)
        return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = "return_plot"), return_plot)

        if (show_plot == TRUE) {
            print(combplot)
        }
        if (save_plot == TRUE) {
            do.call("all_plots_save_function",
                c(list(gobject = gobject, plot_object = combplot, default_save_name = default_save_name), save_param))
        }
        if (return_plot == TRUE) {
            return(combplot)
        }
    } else {
        if (is.null(axis_text_y_size)) axis_text_y_size = 0.8 * 11/.pt
        hmap <- hmap + theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = margin(0, 0, 5.5, 5.5, "pt"))

        # Global binding for variables
        geneOrder <- subset_genes <- .N <- NULL
        geneDT = subset_values_DT[, c("genes"), with = F]
        geneDT = unique(setorder(geneDT, genes))
        geneDT[, `:=`(geneOrder, 1:.N)]
        geneDT[, `:=`(subset_genes, ifelse(genes %in% gene_label_selection,
                                           as.character(genes), ""))]

        axis <- ggplot(data = geneDT, aes(x = 0, y = geneOrder, label = subset_genes)) +
          geom_text_repel(min.segment.length = unit(0, "pt"), color = "grey30", size = axis_text_y_size) +
          scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
          scale_y_continuous(limits = c(0, nrow(geneDT)), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
          theme(panel.background = element_blank(), plot.margin = margin(0, 0, 0, 0, "cm"))

        aligned <- align_plots(clus_pl, empty, empty, hmap + theme(legend.position = "none"), axis, align = "h", axis = "b")
        aligned <- append(aligned, list(get_legend(hmap)))
        combplot = plot_grid(plotlist = aligned, ncol = 3, rel_widths = c(1, 0.2, 0.1), nrow = 2, rel_heights = c(0.2,1))

        show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = "show_plot"), show_plot)
        save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = "save_plot"), save_plot)
        return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = "return_plot"), return_plot)

        if (show_plot == TRUE) {
            print(combplot)
        }
        if (save_plot == TRUE) {
            do.call("all_plots_save_function", c(list(gobject = gobject, plot_object = combplot, default_save_name = default_save_name), save_param))
        }
        if (return_plot == TRUE) {
            return(combplot)
        }
    }
}

#' Plots spatial dimension plots
#'
#' @param gobject Giotto object
#' @param sweepresults sweep results
#' @param pairConnectedOutput paired connected output
#' @param upper_clip upper clip
#' @param lower_clip lower clip
#' @param upper_clip2 upper clip
#' @param lower_clip2 lower clip
#' @param markerList marker list
#' @param point_size point
#' @param full_color_palate color palete
#' @param custom_rowwise_ind custom row wise indices
#' @param additional_string additional string for name
#' @param plot_spatial_dim_plots whether to plot or not
#'
#' @importFrom Giotto getDistinctColors spatDimPlot
#'
#' @return Spatial dim plot
#'
#' @export
plotSpatDimPlots <- function(gobject, sweepresults, pairConnectedOutput,
                             upper_clip = 3, lower_clip = -3,
                             upper_clip2 = 2, lower_clip2 = -2,
                             markerList = NULL,
                             point_size = 0.5,
                             full_color_palate = NULL,
                             custom_rowwise_ind = NULL,
                             additional_string = '',
                             show_plot = TRUE,
                             plot_spatial_dim_plots = TRUE){

    if (is.null(full_color_palate)){
        number_colors = length(pairConnectedOutput$usedLabels)
        full_color_palate = getDistinctColors(n = number_colors)
        names(full_color_palate) = pairConnectedOutput$usedLabels
    }

    minx = min(gobject@spatial_locs$sdimx)
    maxx = max(gobject@spatial_locs$sdimx)
    miny = min(gobject@spatial_locs$sdimy)
    maxy = max(gobject@spatial_locs$sdimy)

    gridSize = dim(sweepresults$clusteringNames)

    # heatmap plotting parameters
    if (is.null(custom_rowwise_ind)){
        custom_rowwise_ind = 1:prod(gridSize)
    }

    for (i in custom_rowwise_ind){ #c(2, 5, 8)
        gobj_hm = copy(gobject)
        gobj_hm2 = copy(gobject)
        currName = names(pairConnectedOutput$clusterings)[i]

        # compute cluster and gene order
        # actually we compute the gene order based only on the correlation of
        # the first half of the rows, ie, on the genes without the neighbour genes.

        # check that the current clustering is indeed the one whose matrix is being used
        current_row_in_grid = floor((i-1)/gridSize[2]+1) # i is counting row wise.
        # print(current_row_in_grid)
        current_col_in_grid = i - gridSize[2]*(floor((i-1)/gridSize[2]))
        # print(current_col_in_grid)
        name_from_grid = sweepresults$clusteringNames[[current_row_in_grid, current_col_in_grid]]
        # print(name_from_grid)
        name_from_clusteringDF = currName
        # print(name_from_clusteringDF)

        if (name_from_grid != name_from_clusteringDF){
            stop("The clustering being referred to in the grid and the clusteringDF should be the same. Probably a row wise vs col wise indexing issue.")
        }
        #!!
        gobj_hm@norm_scaled_expr = sweepresults$scaledNormedMatrixList[[current_row_in_grid, current_col_in_grid]]
        gobj_hm2@norm_scaled_expr = sweepresults$scaledNormedMatrixList[[current_row_in_grid, current_col_in_grid]]


        gobj_hm <- addClippedZCustomExpr(gobj_hm,pairConnectedOutput,sweepresults,
                                         sweepNum = i,
                                         upper_clip = upper_clip,
                                         lower_clip = lower_clip)


        ###### spatial and UMAP plots #########
        dim_point_size = (maxx-minx)/50000
        spat_point_size = point_size
        spatDimPlot(gobject = gobj_hm,
                    cell_color = currName,
                    cell_color_code = full_color_palate,
                    dim_reduction_to_use = 'umap',
                    dim_reduction_name = paste0("joint_2D_", currName),
                    dim_point_size = dim_point_size,
                    spat_point_size = spat_point_size,
                    spat_point_shape = "no_border",
                    dim_point_shape = "no_border",
                    plot_alignment = 'horizontal',
                    show_plot = show_plot,
                    save_param = list(save_name = paste('c', currName, additional_string),
                                      save_folder = 'spatDimPlots',
                                      base_width = max(15,2*(maxx-minx)/3000),base_height = max((maxy-miny)/3000, 7.5)))

    }
}
