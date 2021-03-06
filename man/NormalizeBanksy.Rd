% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computation.R
\name{NormalizeBanksy}
\alias{NormalizeBanksy}
\title{Normalize columns to a normFactor}
\usage{
NormalizeBanksy(
  bank,
  assay = "both",
  normFactor = 100,
  logNorm = FALSE,
  pseudocount = 0.1
)
}
\arguments{
\item{bank}{BanksyObject}

\item{assay}{assay to scale - one of both, own, nbr (default: both)}

\item{normFactor}{normalization factor}

\item{logNorm}{if TRUE, log transforms data}

\item{pseudocount}{for log normalization}
}
\value{
normalized Banksy object
}
\description{
Normalize columns to a normFactor
}
\examples{
# Generate a simulated dataset
d <- simulateDataset()
bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
# Normalize the own.expr matrix
bank <- NormalizeBanksy(bank)

}
