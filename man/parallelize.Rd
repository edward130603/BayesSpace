% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parallelize.R
\name{parallelize}
\alias{parallelize}
\alias{paraLapply}
\title{Parallelization}
\usage{
paraLapply(
  X,
  FUN,
  BPPARAM = NULL,
  cores = 1L,
  type = c("serial", "fork", "sock"),
  verbose = FALSE,
  ...
)
}
\arguments{
\item{X}{Any object for which methods \code{length}, \code{[}, and \code{[[}
are implemented (passed to \code{bplapply}).}

\item{FUN}{The \code{function} to be applied to each element of X
(passed to \code{bplapply}).}

\item{BPPARAM}{An optional \code{BiocParallelParam} instance determining
the parallel back-end to be used during evaluation, or a list of
\code{BiocParallelParam} instances, to be applied in sequence for
nested calls to \code{BiocParallel} functions.}

\item{cores}{The number of threads to use. The results are invariate to the
value of \code{cores}.}

\item{type}{One of "serial", "fork", or "sock". When \code{cores} is one,
\code{type} is always "serial". Both "fork" and "sock" are for
multi-threading. "fork" is faster, but only supports linux and macos.
"sock" supports linux, macos, and windows.}

\item{verbose}{Whether to print debug information or not.}

\item{...}{Additional parameters passed to \code{bplapply}.}
}
\value{
See \code{lapply}.
}
\description{
A convenient wrapper function of \code{BiocParallel} providing easy
  parallelization.
}
