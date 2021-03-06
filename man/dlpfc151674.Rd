% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{dlpfc151674}
\alias{dlpfc151674}
\title{Human dorsolateral prefrontal cortex 10x Visium data (sample 151674)}
\format{
A list with 2 entries:
\describe{
\item{expression}{(matrix) gene expression matrix}
\item{locations}{(data.frame) spot locations in 2D}
}
}
\source{
\url{http://research.libd.org/spatialLIBD/}
}
\usage{
data(dlpfc151674)
}
\description{
This dataset comprises 10x Visium profiling of the human DLPFC (Maynard et
al. (2021) Transcriptome-scale spatial gene expression in the human
dorsolateral prefrontal cortex. doi.org: 10.1038/s41593-020-00787-0) for
subject 3 sample 151674. Gene expression for 3,673 spots in 2 dimensions was
count normalized and subset to 7,319 genes obtained by taking the union of
the top 3,000 genes in each sample for subject 3.
}
\keyword{datasets}
