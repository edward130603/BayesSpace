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

#' @export
#' @rdname readChain
readChain <- function(sce, params=NULL) {
  if (!("chain.h5" %in% names(metadata(sce)))) {
    stop("Path to chain file not available in object metadata.")
  }
  
  .read_chain(metadata(sce)$chain.h5, params)
}