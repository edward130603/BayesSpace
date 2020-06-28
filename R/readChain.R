#' Read MCMC chain associated with a BayesSpace clustering or enhancement
#' 
#' BayesSpace stores the MCMC chain associated with a clustering or enhancement 
#' on disk in an HDF5 file. The `readChain` function reads any parameters 
#' specified by the user into a `coda::mcmc` object compatible with TidyBayes.
#' 
#' @details 
#' To interact with the HDF5 file directly, obtain the filename from the 
#' SingleCellExperiment's metadata: `metadata(sce)$chain.h5`. Each parameter is 
#' stored as a separate dataset in the file, and is represented as a matrix of 
#' size (n_iterations x n_parameter_indices). 
#' 
#' @param sce SingleCellExperiment with a file path stored in its metadata.
#' @param params List of model parameters to read
#' 
#' @name readChain
NULL

#' @importFrom rhdf5 h5createFile
#' @importFrom rhdf5 h5createDataset
#' @importFrom rhdf5 h5write
.write_chain <- function(chain, h5.fname=NULL, chunk.length=1000) {
  if (is.null(h5.fname)) {
    h5.fname <- tempfile(fileext=".h5")
  }
  
  h5createFile(h5.fname)
  
  # TODO: clean formats ahead of time (n_reps x n_param_idxs)
  # TODO: clean matrix colnames (parname[i,j])
  # for (param in names(chain)) {
  for (param in c("z", "mu")) {
    dims <- dim(chain[[param]])
    h5createDataset(h5.fname, param, dims, 
                           chunk=c(chunk.length, dims[2]))
    
    h5write(chain[[param]], h5.fname, param)
  }
  
  h5.fname
}

#' @importFrom coda mcmc
#' @importFrom purrr map
.read_chain <- function(h5.fname, params=NULL) {
  # TODO: add option to subset last n rows/iterations
  
  if (is.null(params)) {
    params <- h5ls(h5.fname)$name
  }
  
  xs <- map(params, function(x) {as.matrix(h5read(h5.fname, x))})
  x <- do.call(cbind, xs)
  
  # TODO: specify thinning interval; start/end if we drop burn-in earlier
  mcmc(x)
}

.make_index_names <- function(name, m, n=NULL, dim=1) 
# Make colnames for parameter indices. 
{
  if (is.null(n)) {
    paste0(name, "[", rep(1:m), "]")
  } else if (dim == 1) {
    paste0(name, "[", rep(1:m, each=n), ",", rep(1:n, m), "]")
  } else {
    paste0(name, "[", rep(1:m, n), ",", rep(1:n, each=m), "]")
  }
}

#' @importFrom purrr map
.clean_chain <- function(out, method=c("cluster", "enhance")) 
# Tidy C++ outputs before writing to disk.
# 1) Convert each parameter to matrix (n_iterations x n_indices)
# 2) Add appropriate colnames 
# 3) Thin evenly (for enhance)
{
  d <- ncol(out$lambda[[1]])
  q <- ncol(out$mu) / d
    
  colnames(out$z) <- .make_index_names("z", ncol(out$z))
  colnames(out$mu) <- .make_index_names("mu", q, d)
  
  lambdas <- map(out$lambda, function(x) as.vector(t(x)))
  out$lambda <- do.call(rbind, lambdas)
  colnames(out$lambda) <- .make_index_names("lambda", d, d)
  
  out$plogLik <- as.matrix(out$plogLik)
  colnames(out$plogLik) <- c("pLogLikelihood")
  
  out
} 

#' @export
#' @rdname readChain
readChain <- function(sce, params=NULL) {
  if (!("chain.h5" %in% names(metadata(sce)))) {
    stop("Path to chain file not available in object metadata.")
  }
  
  .read_chain(metadata(sce)$chain.h5, params)
}