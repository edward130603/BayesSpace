% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{getRDS}
\alias{getRDS}
\title{Download a processed sample from our S3 bucket}
\usage{
getRDS(dataset, sample, cache = TRUE)
}
\arguments{
\item{dataset}{Dataset identifier}

\item{sample}{Sample identifier}

\item{cache}{If true, cache the dataset locally with \code{BiocFileCache}.
Otherwise, download directly from our S3 bucket. Caching saves time on
subsequent loads, but consumes disk space.}
}
\value{
sce A SingleCellExperiment with positional information in colData and
  PCs based on the top 2000 HVGs
}
\description{
Datasets are cached locally using \code{BiocFileCache}. The first time using
this function, you may need to consent to creating a BiocFileCache directory
if one does not already exist.
}
\examples{
sce <- getRDS("2018_thrane_melanoma", "ST_mel1_rep2", cache=FALSE)

}
