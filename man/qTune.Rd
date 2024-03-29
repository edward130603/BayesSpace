% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qTune.R
\name{qTune}
\alias{qTune}
\alias{qPlot}
\title{Tuning the choice of q (number of clusters) before running spatialCluster}
\usage{
qPlot(sce, qs = seq(3, 7), force.retune = FALSE, ...)

qTune(sce, qs = seq(3, 7), burn.in = 100, nrep = 1000, cores = 1L, ...)
}
\arguments{
\item{sce}{A SingleCellExperiment object containing the spatial data.}

\item{qs}{The values of q to evaluate.}

\item{force.retune}{If specified, existing tuning values in \code{sce} will
be overwritten.}

\item{...}{Other parameters are passed to \code{spatialCluster()}.}

\item{burn.in, nrep}{Integers specifying the range of repetitions to
compute.}

\item{cores}{The number of threads to use. The results are invariate to the
value of \code{cores}.}
}
\value{
\code{qTune()} returns a modified \code{sce} with tuning log
  likelihoods stored as an attribute named \code{"q.logliks"}.

  \code{qPlot()} returns a ggplot object.
}
\description{
Before running \code{spatialCluster()}, we recommend tuning the choice of
\code{q} by choosing the \code{q} that minimizes the model's negative log
likelihood over early iterations. \code{qTune()} computes the average
negative log likelihood for a range of q values over iterations 100:1000, and
\code{qPlot()} displays the results.
}
\details{
\code{qTune()} takes the same parameters as \code{spatialCluster()} and will
  run the MCMC clustering algorithm up to \code{nrep} iterations for each
  value of \code{q}. The first \code{burn.in} iterations are discarded as
  burn-in and the log likelihood is averaged over the remaining iterations.

\code{qPlot()} plots the computed negative log likelihoods as a function of
  q. If \code{qTune()} was run previously, i.e. there exists an attribute of
  \code{sce} named \code{"q.logliks"}, the pre-computed results are
  displayed. Otherwise, or if \code{force.retune} is specified,
  \code{qplot()} will automatically run \code{qTune()} before plotting (and
  can take the same parameters as \code{spatialCluster()}.
}
\examples{
set.seed(149)
sce <- exampleSCE()
sce <- qTune(sce, seq(3, 7), burn.in = 10, nrep = 100)
qPlot(sce)

}
