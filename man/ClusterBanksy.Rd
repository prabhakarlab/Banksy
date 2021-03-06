% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{ClusterBanksy}
\alias{ClusterBanksy}
\title{Cluster based on joint expression matrix}
\usage{
ClusterBanksy(
  bank,
  lambda = 0.25,
  pca = TRUE,
  npcs = 30,
  method = c("leiden", "louvain", "mclust", "kmeans"),
  k.neighbors = NULL,
  resolution = NULL,
  leiden.iter = -1,
  mclust.G = NULL,
  kmeans.centers = NULL,
  kmeans.iter.max = 10,
  seed = 42,
  ...
)
}
\arguments{
\item{bank}{BanksyObject}

\item{lambda}{weighting parameter}

\item{pca}{if TRUE, runs clustering on PCA, else runs on Banksy matrix}

\item{npcs}{number of pcs to use for clustering}

\item{method}{one of leiden, louvain, mclust, kmeans}

\item{k.neighbors}{leiden / louvain - parameter for constructing shared nearest neighbor network}

\item{resolution}{leiden - parameter used for clustering}

\item{leiden.iter}{leiden - number of leiden iterations}

\item{mclust.G}{mclust - number of mixture components G}

\item{kmeans.centers}{kmeans - number of clusters}

\item{kmeans.iter.max}{kmeans - max number of iterations}

\item{seed}{seed}

\item{...}{to pass to methods}
}
\value{
BanksyObject with cluster labels in meta.data
}
\description{
Cluster based on joint expression matrix
}
\examples{
# Generate a simulated dataset
d <- simulateDataset()
bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
bank <- NormalizeBanksy(bank)
bank <- ComputeBanksy(bank)
bank <- ScaleBanksy(bank)
bank <- RunPCA(bank, lambda = 0.3)
bank <- ClusterBanksy(bank, lambda = 0.3, npcs = 20, k.neighbors = 50, resolution = 0.5)

}
