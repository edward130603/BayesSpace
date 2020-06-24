#' BayesSpace: A package for processing spatial transcriptomes
#' 
#' The BayesSpace package provides two primary functions: `cluster()` and `deconvolve()`.
#' 
#' @section Clustering:
#' The clustering algorithm groups neighboring spots together by expression profile.
#' 
#' @section Deconvolution:
#' The deconvolution algorithm deconvolves spot-level expression signals into imputed cell-level signals.
#' 
#' @docType package
#' @name BayesSpace
#' 
#' @export cluster
#' @export deconvolve
#' @export predictExpression
#' @export spatialCluster
#' @import SingleCellExperiment
#' @importFrom Rcpp evalCpp
#' @importFrom stats cov rWishart
#' @useDynLib BayesSpace
#' @useDynLib BayesSpace, .registration=TRUE
NULL