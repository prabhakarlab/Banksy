% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotReduction}
\alias{plotReduction}
\title{Plot dimensionality reduction assays}
\usage{
plotReduction(
  bank,
  reduction,
  components = c(1, 2),
  by = NULL,
  type = c("discrete", "continuous"),
  pt.size = 0.5,
  pt.alpha = 0.7,
  col.midpoint = NULL,
  col.low = "blue",
  col.mid = "gray95",
  col.high = "red",
  col.discrete = NULL,
  main = NULL,
  main.size = 5,
  legend = TRUE,
  legend.text.size = 6,
  legend.pt.size = 3
)
}
\arguments{
\item{bank}{BanksyObject}

\item{reduction}{reduction to visualize}

\item{components}{which components to visualize}

\item{by}{groupings for cells}

\item{type}{type of groupings (one of 'discrete' or 'continuous')}

\item{pt.size}{size of points}

\item{pt.alpha}{transparency of points}

\item{col.midpoint}{for continuous labels - midpoint of color gradient}

\item{col.low}{for continuous labels - color gradient low}

\item{col.mid}{for continuous labels - color gradient mid}

\item{col.high}{for continuous labels - color gradient high}

\item{col.discrete}{for discrete labels - colors for groupings}

\item{main}{title}

\item{main.size}{size of title}

\item{legend}{show legend}

\item{legend.text.size}{size of legend text}

\item{legend.pt.size}{size of legend point}
}
\value{
Dimensionality reduction plot
}
\description{
Plot dimensionality reduction assays
}
\examples{
# Generate a simulated dataset
d <- simulateDataset()
bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
bank <- NormalizeBanksy(bank)
bank <- ScaleBanksy(bank)
bank <- ComputeBanksy(bank)
bank <- RunPCA(bank, lambda = 0.2)
names(reduction(bank))
plotReduction(bank, reduction = 'pca_0.2', by = 'Label', type = 'discrete')

}
