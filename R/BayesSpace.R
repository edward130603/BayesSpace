#' BayesSpace: A package for processing spatial transcriptomes
#' 
#' TODO: flesh out package description (or link to vignette a la scran)
#' 
#' The BayesSpace package provides two primary functions:
#' \code{spatialCluster()} and \code{spatialEnhance()}.
#' 
#' @section Clustering:
#' The clustering algorithm groups neighboring spots together by expression
#' profile.
#' 
#' @section Enhancement:
#' The enhancement algorithm deconvolves spot-level expression into imputed
#' cell-level signals.
#' 
#' @docType package
#' @name BayesSpace
#' 
#' @importFrom Rcpp evalCpp
#' @useDynLib BayesSpace
#' @useDynLib BayesSpace, .registration=TRUE
NULL
