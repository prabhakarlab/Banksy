% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotHeatmap}
\alias{plotHeatmap}
\title{Plot Heatmap (wrapper for ComplexHeatmap)}
\usage{
plotHeatmap(
  bank,
  assay = "own.expr",
  dataset = NULL,
  lambda = NULL,
  cells = NULL,
  features = NULL,
  col = NULL,
  col.breaks = NULL,
  col.discrete = NULL,
  cluster.row = TRUE,
  cluster.column = FALSE,
  row.dend = FALSE,
  column.dend = FALSE,
  cex.row = 4,
  annotate = FALSE,
  annotate.by = NULL,
  order.by = NULL,
  annotation.name = TRUE,
  annotation.size = 6,
  annotation.pos = "right",
  name = "Expression",
  barplot.by = NULL,
  barplot.border = FALSE,
  barplot.width = 0.6,
  max.cols = NULL,
  seed = 42,
  rasterize = FALSE,
  ...
)
}
\arguments{
\item{bank}{BanksyObject}

\item{assay}{assay to plot heatmap (one of own.expr, nbr.expr or banksy)}

\item{dataset}{dataset to plot heatmap}

\item{lambda}{lambda if assay == banksy}

\item{cells}{specific cells to plot}

\item{features}{specific features to plot}

\item{col}{colours to use in heatmap}

\item{col.breaks}{color breaks to use in heatmap (same number as col is
specified)}

\item{col.discrete}{cluster colors, named color array, names are clusters,
values are colors in hexadecimal}

\item{cluster.row}{cluster rows}

\item{cluster.column}{cluster columns}

\item{row.dend}{draw row dendrograms}

\item{column.dend}{draw column dendrograms}

\item{cex.row}{row label size}

\item{annotate}{add annotation for cells}

\item{annotate.by}{metadata to annotate cells by}

\item{order.by}{metadata to order cells by (one of annotate.by)}

\item{annotation.name}{show annotation name}

\item{annotation.size}{size of annotation labels}

\item{annotation.pos}{position of annotation labels}

\item{name}{name of heatmap legend}

\item{barplot.by}{metadata to plot barplots (numeric)}

\item{barplot.border}{show borders of barplot}

\item{barplot.width}{barplot width}

\item{max.cols}{max columns to display - will subsample if not NULL}

\item{seed}{used for sampling}

\item{rasterize}{rasterize if TRUE}

\item{...}{parameters to pass to ComplexHeatmap::Heatmap}
}
\value{
Heatmap of class ComplexHeatmap
}
\description{
Plot Heatmap (wrapper for ComplexHeatmap)
}
\examples{
# Generate a simulated dataset
d <- simulateDataset(n_cells = c(20,20,20), n_genes = 30)
bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
bank <- NormalizeBanksy(bank)
bank <- ScaleBanksy(bank)
plotHeatmap(bank, annotate = TRUE, annotate.by = 'Label')

}
