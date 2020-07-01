#' Read MCMC chain associated with a BayesSpace clustering or enhancement
#' 
#' BayesSpace stores the MCMC chain associated with a clustering or enhancement
#' on disk in an HDF5 file. The \code{mcmcChain()} function reads any parameters
#' specified by the user into a \code{coda::mcmc} object compatible with
#' TidyBayes.
#' 
#' @details 
#' To interact with the HDF5 file directly, obtain the filename from the
#' SingleCellExperiment's metadata: \code{metadata(sce)$chain.h5}. Each
#' parameter is stored as a separate dataset in the file, and is represented as
#' a matrix of size (n_iterations x n_parameter_indices).
#' 
#' @param sce SingleCellExperiment with a file path stored in its metadata.
#' @param params List of model parameters to read
#' 
#' @return Returns an \code{mcmc} object containing the values of the requested
#'   parameters over the constructed chain.
#' 
#' @name mcmcChain
NULL

#' @importFrom rhdf5 h5createFile h5createDataset h5write
.write_chain <- function(chain, h5.fname = NULL, params = NULL, 
    chunk.length = 1000, chunk.width=4000) {
    
    if (is.null(h5.fname)) {
        h5.fname <- tempfile(fileext=".h5")
    }
    
    if (is.null(params)) {
        params <- names(chain)
    }
    
    h5createFile(h5.fname)
    
    for (param in params) {
        dims <- dim(chain[[param]])
        chunk <- c(min(chunk.length, dims[1]), dims[2])
        h5createDataset(h5.fname, param, dims, chunk=chunk)
        
        attr(chain[[param]], "colnames") <- colnames(chain[[param]])
        ## TODO: write colnames manually to avoid warnings about dimnames
        suppressWarnings(h5write(chain[[param]], h5.fname, param,
            write.attributes=TRUE))
    }
    
    h5.fname
}

#' @importFrom rhdf5 h5ls h5read
#' @importFrom coda mcmc
#' @importFrom purrr map
.read_chain <- function(h5.fname, params = NULL) {
    ## TODO: add option to subset last n rows/iterations
    
    if (is.null(params)) {
        params <- h5ls(h5.fname)$name
    }
    
    .read_param <- function(x) {
        x <- as.matrix(h5read(h5.fname, x, read.attributes=TRUE))
        colnames(x) <- attr(x, "colnames")
        attr(x, "colnames") <- NULL
        x
    }
    
    xs <- map(params, .read_param)
    x <- do.call(cbind, xs)
    
    ## TODO: specify thinning interval; start/end if we drop burn-in earlier
    mcmc(x)
}

## Make colnames for parameter indices.
.make_index_names <- function(name, m, n = NULL, dim = 1) {
    if (is.null(n)) {
        paste0(name, "[", rep(seq_len(m)), "]")
    } else if (dim == 1) {
        paste0(name, "[", rep(seq_len(m), each=n), ",", rep(seq_len(n), m), "]")
    } else {
        paste0(name, "[", rep(seq_len(m), n), ",", rep(seq_len(n), each=m), "]")
    }
}

## Tidy C++ outputs before writing to disk.
##  1) Convert each parameter to matrix (n_iterations x n_indices) 
##  2) Add appropriate colnames 
##  3) Thin evenly (for enhance)
#' @importFrom purrr map
.clean_chain <- function(out, method = c("cluster", "enhance"), thin=100) 
{
    n_iter <- nrow(out$z)
    n <- ncol(out$z)
    d <- ncol(out$lambda[[1]])
    q <- ncol(out$mu)/d
    
    colnames(out$z) <- .make_index_names("z", n)
    colnames(out$mu) <- .make_index_names("mu", q, d)
    out$lambda <- .flatten_matrix_list(out$lambda, "lambda", d, d)
    
    ## Include function-specific chain parameters
    if (method == "cluster") {
        out$plogLik <- as.matrix(out$plogLik)
        colnames(out$plogLik) <- c("pLogLikelihood")
    } else if (method == "enhance") {
        colnames(out$weights) <- .make_index_names("weights", n)
        out$Y <- .flatten_matrix_list(out$Y, "Y", n, d)
        out$Ychange <- as.matrix(out$Ychange)
        colnames(out$Ychange) <- c("Ychange")
    }
    
    ## TODO: optionally thin cluster output too
    if (method == "enhance") {
        for (param in c("z", "mu", "lambda", "weights")) {
            out[[param]] <- out[[param]][seq(thin, n_iter, thin), ]
        }
        
        ## Subset of a one-column matrix is a vector, not a matrix
        out$Ychange <- as.matrix(out$Ychange[seq(thin, n_iter, thin), ])
        colnames(out$Ychange) <- c("Ychange")
        
        ## Y is thinned inside `deconvolve` but includes starting values
        out$Y <- out$Y[seq(2, nrow(out$Y)), ]
    }
    
    out
}

.flatten_matrix_list <- function(xs, ...) {
    xs <- map(xs, function(x) as.vector(t(x)))
    x <- do.call(rbind, xs)
    colnames(x) <- .make_index_names(...)
    
    x
}

#' @export
#' @rdname mcmcChain
mcmcChain <- function(sce, params = NULL) {
    if (!("chain.h5" %in% names(metadata(sce)))) {
        stop("Path to chain file not available in object metadata.")
    }
    
    .read_chain(metadata(sce)$chain.h5, params)
}
