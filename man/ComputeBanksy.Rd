% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computation.R
\name{ComputeBanksy}
\alias{ComputeBanksy}
\title{Compute Banksy Matrices}
\usage{
ComputeBanksy(
  bank,
  spatialMode = "kNN_r",
  k_geom = 10,
  n = 2,
  sigma = 1.5,
  alpha = 0.05,
  kspatial = 10,
  dimensions = "all",
  verbose = FALSE
)
}
\arguments{
\item{bank}{BanksyObject}

\item{spatialMode}{spatial mode to use - one of kNN_r, kNN_rn, kNN_rank,
kNN_unif, rNN_gauss}

\item{k_geom}{kNN parameter - number of neighbors to use}

\item{n}{kNN_rn parameter - exponent of radius}

\item{sigma}{rNN parameter - standard deviation of Gaussian kernel}

\item{alpha}{rNN parameter - determines radius used (larger alphas give
smaller radii)}

\item{kspatial}{rNN parameter - number of neighbors to use}

\item{dimensions}{dimensions to use when computing neighborhood - one of
column name in cell.locs, or 'all'}

\item{verbose}{messages}
}
\value{
BanksyObject
}
\description{
Compute Banksy Matrices
}
\examples{
d <- simulateDataset()
bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
bank <- NormalizeBanksy(bank)
# Compute neighbors 
bank <- ComputeBanksy(bank)

}
