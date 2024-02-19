#' BayesSpace: A package for processing spatial transcriptomes
#' 
#' Tools for clustering and enhancing the resolution of spatial gene expression
#' experiments. BayesSpace clusters a low-dimensional representation of the gene
#' expression matrix, incorporating a spatial prior to encourage neighboring
#' spots to cluster together. The method can enhance the resolution of the
#' low-dimensional representation into "sub-spots", for which features such as
#' gene expression or cell type composition can be imputed.
#' 
#' @details
#' For an overview of the functionality provided by the package, please see the
#' vignette:
#' \code{vignette("BayesSpace", package="BayesSpace")}
#'
#' 
#' @name BayesSpace
#' 
#' @keywords internal
#'
#' @importFrom Rcpp evalCpp
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SingleCellExperiment altExpNames
#' @useDynLib BayesSpace
#' @useDynLib BayesSpace, .registration=TRUE
"_PACKAGE"
