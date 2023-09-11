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
#' a matrix of size (n_iterations x n_parameter_indices). Parameter choices
#' for the spot-level clustering include:
#' * \code{z} (cluster assignments)
#' * \code{weights} (\eqn{w_i})
#' * \code{mu} (mean vectors) 
#' * \code{lambda} (precision matrix) 
#' * \code{plogLik} (pseudo-log-likelihood) 
#' 
#' Parameter choices for the subspot-level enhanced clustering include: 
#' * \code{z} (cluster assignments)
#' * \code{weights} (\eqn{w_i})
#' * \code{Y} (enhanced PCs)
#' * \code{mu} (mean vectors)
#' * \code{lambda} (precision matrix) 
#' * \code{Ychange} (acceptance rate for the jittering of PCs)
#' 
#' For best results, \code{Ychange} should average between 0.25 and 0.40.
#' 
#' @param sce SingleCellExperiment with a file path stored in its metadata.
#' @param params List of model parameters to read
#' 
#' @return Returns an \code{mcmc} object containing the values of the requested
#'   parameters over the constructed chain.
#'   
#' @examples
#' set.seed(149)
#' sce <- exampleSCE()
#' sce <- spatialCluster(sce, 7, nrep=100, burn.in=10, save.chain=TRUE)
#' chain <- mcmcChain(sce)
#' removeChain(sce)
#' 
#' @name mcmcChain
#' @md
NULL

#' @importFrom rhdf5 h5createFile h5createDataset h5write
.write_chain <- function(chain, h5.fname = NULL, params = NULL, 
    chunk.length = 1000) {
    
    if (is.null(h5.fname)) {
        h5.fname <- tempfile(fileext=".h5")
    }
    
    if (is.null(params)) {
        params <- names(chain)
    }
    
    h5createFile(h5.fname)
    
    for (par.name in params) {
        param <- chain[[par.name]]
        dims <- dim(param)
        chunk <- c(min(chunk.length, dims[1]), dims[2])
        h5createDataset(h5.fname, par.name, dims, chunk=chunk)
        
        attr(param, "dims") <- .infer_param_dims(colnames(param))

        suppressWarnings(h5write(param, h5.fname, par.name, write.attributes=TRUE))
    }
    
    h5.fname
}

#' Infer original dimensions of parameter (per iteration) from colnames
#'
#' Used to avoid writing colnames directly to HDF5 as attribute, which fails
#' for large parameters (e.g. Y)
#' 
#' @param cnames List of column names
#' @return Numeric vector (nrow, ncol)
#' 
#' @keywords internal
.infer_param_dims <- function(cnames) {
    n_idxs <- length(cnames)
    dims <- list()
    
    if (n_idxs == 1) {
        dims <- c(1, 1)
    } else if (!grepl(",", cnames[n_idxs])) {
        dims <- c(n_idxs, 1)
    } else {
        dims <- gsub("[a-zA-Z_\\.]|\\[|\\]", "", cnames[n_idxs])
        dims <- as.numeric(strsplit(dims, ",")[[1]])
    }
    
    dims
}

#' Load saved chain from disk.
#' 
#' @param h5.fname Path to hdf5 file containing chain
#' @param params List of parameters to read from file (will read all by default)
#' 
#' @return MCMC chain, represented as a \code{coda::mcmc} object
#' 
#' @keywords internal
#' 
#' @importFrom rhdf5 h5ls h5read
#' @importFrom coda mcmc
#' @importFrom purrr map
.read_chain <- function(h5.fname, params = NULL, is.enhanced = FALSE) {
    if (is.null(params)) {
        params <- h5ls(h5.fname)$name
    }
    
    .read_param <- function(par.name) {
        param <- as.matrix(h5read(h5.fname, par.name, read.attributes=TRUE))
        dims <- attr(param, "dims")
        colnames(param) <- .make_index_names(par.name, dims[1], dims[2])
        
        param
    }
    
    xs <- map(params, .read_param)
    x <- do.call(cbind, xs)
    
    ## Enhanced chain includes initialization and is thinned to every 100 iters
    ## Cluster chain does not include init and is not thinned
    if (is.enhanced)
        mcmc(x, start=0, end=(nrow(x) - 1) * 100, thin=100)
    else
        mcmc(x)
}

#' Make colnames for parameter indices.
#' 
#' Scalar parameters are named \code{"name"}.
#' Vector parameters are named \code{"name[i]"}.
#' Matrix parameters are named \code{"name[i,j]"}.
#' 
#' @param name Parameter name
#' @param m,n Dimensions of parameter (m=nrow, n=ncol)
#' @param dim Dimensionality of parameter (0=scalar, 1=vector, 2=matrix)
#' 
#' @return List of names for parameter values
#' 
#' @keywords internal
.make_index_names <- function(name, m = NULL, n = NULL, dim = 1) {
    if (is.null(m) || m == 1) {
        name
    } else if (is.null(n) || n == 1) {
        paste0(name, "[", rep(seq_len(m)), "]")
    } else if (dim == 1) {
        paste0(name, "[", rep(seq_len(m), each=n), ",", rep(seq_len(n), m), "]")
    } else {
        paste0(name, "[", rep(seq_len(m), n), ",", rep(seq_len(n), each=m), "]")
    }
}

#' Tidy C++ outputs before writing to disk.
#' 
#' 1) Convert each parameter to matrix (n_iterations x n_indices) 
#' 2) Add appropriate colnames 
#' 3) Thin evenly (for enhance)
#'
#' @param out List returned by \code{cluster()} or \code{deconvolve()}.
#' @param method Whether the output came from clustering or enhancement.
#'   (Different params are included in each.)
#' @param thin Thinning rate. Some enhanced parameters are thinned within C++
#'   loop, others (\code{mu} and \code{Ychange}) need to be thinned afterwards.
#'   
#' @return List with standardized parameters
#'   
#' @keywords internal
#' 
#' @importFrom purrr map
.clean_chain <- function(out, method = c("cluster", "enhance"), thin=100) 
{
    method <- match.arg(method)
    n_iter <- nrow(out$z)  # this is technically n_iters / 100 for enhance

    ## Only one iteration included; need to cast vectors back to matrices
    if (is.null(n_iter)) {
        out$z <- matrix(out$z, nrow=1)
        out$mu <- matrix(out$mu, nrow=1)
        if ("weights" %in% names(out)) {
            out$weights <- matrix(out$mu, nrow=1)
        }
    }

    n <- ncol(out$z)
    d <- ncol(out$lambda[[1]])
    q <- ncol(out$mu)/d

    colnames(out$z) <- .make_index_names("z", n)
    colnames(out$mu) <- .make_index_names("mu", q, d)
    out$lambda <- .flatten_matrix_list(out$lambda, "lambda", d, d)
    
    if ("weights" %in% names(out)) {
        colnames(out$weights) <- .make_index_names("weights", n)
    }
    
    ## Include function-specific chain parameters
    out$plogLik <- as.matrix(out$plogLik)
    colnames(out$plogLik) <- c("pLogLikelihood")
    if (method == "enhance") {
        out$Y <- .flatten_matrix_list(out$Y, "Y", n, d)
        
        out$Ychange <- as.matrix(out$Ychange)
        colnames(out$Ychange) <- c("Ychange")
        
        out$jitterScale <- as.matrix(out$jitterScale)
        colnames(out$jitterScale) <- c("jitterScale")
    }
    
    if (method == "enhance") {
        ## manually thin mu until updated in c++; 
        ## keep init values for consistency with others
        thinned_idx <- c(1, seq(thin, (n_iter - 1) * thin, thin))
        out$mu <- out$mu[thinned_idx, ]
        
        ## Subset of a one-column matrix is a vector, not a matrix
        out$Ychange <- as.matrix(out$Ychange[thinned_idx, ])
        colnames(out$Ychange) <- c("Ychange")
        
        out$plogLik <- as.matrix(out$plogLik[thinned_idx, ])
        colnames(out$plogLik) <- c("pLogLikelihood")
    }
    
    out
}

#' Convert a list of matrices to a single matrix, where each row is a flattened
#' matrix from the original list
#' 
#' @param xs List of matrices
#' @return Matrix
#' 
#' @keywords internal
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

    is.enhanced <- metadata(sce)$BayesSpace.data$is.enhanced
    .read_chain(metadata(sce)$chain.h5, params, is.enhanced)
}

#' @export
#' @rdname mcmcChain
removeChain <- function(sce) {
    if ("chain.h5" %in% names(metadata(sce))) {
        if (file.exists(metadata(sce)$chain.h5)) {
            file.remove(metadata(sce)$chain.h5)
        }
        metadata(sce)$chain.h5 <- NULL
    }
}
