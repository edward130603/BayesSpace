#' Spatial clustering
#'
#' Cluster a spatial expression dataset.
#' 
#' @param sce A SingleCellExperiment object containing the spatial data.
#' @param q The number of clusters.
#' @param platform Spatial transcriptomic platform. Specify 'Visium' for hex 
#'   lattice geometry or 'ST' for square lattice geometry.
#' @param use.dimred Name of a reduced dimensionality result in 
#'   \code{reducedDims(sce)}. If provided, cluster on these features directly. 
#' @param d Number of top principal components to use when clustering.
#' @param init Initial cluster assignments for spots.
#' @param init.method If \code{init} is not provided, cluster the top \code{d} 
#'   PCs with this method to obtain initial cluster assignments.
#' @param model Error model. ('normal' or 't')
#' @param precision Covariance structure. ('equal' or 'variable' for EEE and 
#'   VVV covariance models, respectively.)
#' @param nrep The number of MCMC iterations.
#' @param gamma Smoothing parameter. (Values in range of 1-3 seem to work well.)
#' @param mu0 Prior mean hyperparameter for mu. If not provided, mu0 is set to
#'   the mean of PCs over all spots.
#' @param lambda0 Prior precision hyperparam for mu. If not provided, lambda0
#'   is set to a diagonal matrix \eqn{0.01 I}.
#' @param alpha Hyperparameter for Wishart distributed precision lambda.
#' @param beta Hyperparameter for Wishart distributed precision lambda.
#' @param save.chain If true, save the MCMC chain to an HDF5 file. 
#' @param chain.fname File path for saved chain. Tempfile used if not provided.
#' 
#' @return Returns a modified \code{sce} with cluster assignments stored in
#'   \code{colData} under the name \code{spatial.cluster}.
#'         
#' @details 
#' TODO describe method in detail
#' TODO add details or link to mcmcChain
#' 
#' @examples
#' set.seed(149)
#' sce <- exampleSCE()
#' sce <- spatialCluster(sce, 7, nrep=200)
#'
#' @name spatialCluster
NULL

#' Wrapper around C++ \code{iterate_*()} functions
#' 
#' @return List of clustering parameter values at each iteration
#' 
#' @keywords internal
#' @importFrom purrr map
cluster <- function(Y, q, df_j, init = rep(1, nrow(Y)),
    model = c("t", "normal"), precision = c("equal", "variable"),
    mu0 = colMeans(Y), lambda0 = diag(0.01, nrow = ncol(Y)), 
    gamma = 3, alpha = 1, beta = 0.01, nrep = 1000) {
    
    Y <- as.matrix(Y)
    d <- ncol(Y)
    n <- nrow(Y)
    
    if ((length(mu0) != d) | (ncol(lambda0) != d)) {
        stop("Dimensions of mu0 or lambda0 do not match input data Y")
    }
    
    model <- match.arg(model)
    precision <- match.arg(precision)
    
    if (q == 1) {
        return(list(z=matrix(rep(1, n), nrow=1)))
    }
    
    message("Fitting model...")
    if (model == "normal") {
        if (precision == "equal") {
            cluster.FUN <- iterate
        } else if (precision == "variable") {
            cluster.FUN <- iterate_vvv
        }
    } else if (model == "t") {
        if (precision == "equal") {
            cluster.FUN <- iterate_t
        } else if (precision == "variable") {
            cluster.FUN <- iterate_t_vvv
        }
    }
    
    cluster.FUN(Y=as.matrix(Y), df_j=df_j, nrep=nrep, n=n, d=d, gamma=gamma,
        q=q, init=init, mu0=mu0, lambda0=lambda0, alpha=alpha, beta=beta)
}

## TODO make generic for SCE/matrix instead of wrapping cluster() ?
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors metadata metadata<-
#' 
#' @export
#' @rdname spatialCluster
spatialCluster <- function(sce, q, use.dimred = "PCA", d = 15,
    platform=c("Visium", "ST"),
    init = NULL, init.method = c("mclust", "kmeans"),
    model = c("t", "normal"), precision = c("equal", "variable"), 
    nrep = 50000, gamma = 3, mu0 = NULL, lambda0 = NULL, alpha = 1, 
    beta = 0.01, save.chain = FALSE, chain.fname = NULL) {
    
    if (!(use.dimred %in% reducedDimNames(sce))) 
        stop(sprintf("reducedDim %s not found in input SCE", use.dimred))

    ## Get PCs
    Y <- reducedDim(sce, use.dimred)
    d <- min(ncol(Y), d)
    Y <- Y[, seq_len(d)]
    
    ## Get indices of neighboring spots, and initialize cluster assignments
    df_j <- .find_neighbors(sce, match.arg(platform))
    init <- .init_cluster(Y, q, init, init.method)
    
    ## TODO: pass these through with ...
    model <- match.arg(model)
    precision <- match.arg(precision)
    mu0 <- if (is.null(mu0)) colMeans(Y) else mu0
    lambda0 <- if (is.null(lambda0)) diag(0.01, ncol(Y)) else lambda0
    
    ## Run clustering
    results <- cluster(Y, q, df_j, init=init, 
        model=model, precision=precision, mu0=mu0, 
        lambda0=lambda0, gamma=gamma, alpha=alpha, beta=beta, nrep=nrep)
    
    ## Save MCMC chain
    if (save.chain) {
        results <- .clean_chain(results)
        metadata(sce)$chain.h5 <- .write_chain(results, chain.fname)
    }
    
    ## Save metadata (TODO: add neighbors?)
    sce$cluster.init <- init
    metadata(sce)$BayesSpace.platform <- platform
    metadata(sce)$BayesSpace.is_enhanced <- FALSE
    
    ## Save modal cluster assignments
    ## NOTE: swap below code for this to test against refactoring
    ## colData(sce)$spatial.cluster <- apply(results$z[900:1000, ], 2, Mode)
    iter_from <- ifelse(nrep < 2000, max(2, nrep - 1000), 1000)
    msg <- "Calculating labels using iterations %d through %d"
    message(sprintf(msg, iter_from, nrep))
    labels <- apply(results$z[iter_from:nrep, ], 2, Mode)
    colData(sce)$spatial.cluster <- unname(labels)
    
    sce
}

#' Find neighboring spots based on array coordinates
#' 
#' @param sce SingleCellExperiment
#' @param platform If "Visium", select six neighboring spots around center; if
#'   "ST", select four adjacent spots.
#' @return \code{df_j} a list of neighbor indices (zero-indexed) for each spot
#' 
#' @keywords internal
.find_neighbors <- function(sce, platform) {
    if (platform == "Visium") {
        ## Spots to left and right, two above, two below
        offsets <- data.frame(x.offset=c(-2, 2, -1,  1, -1, 1),
                              y.offset=c( 0, 0, -1, -1,  1, 1))
    } else if (platform == "ST") {
        ## L1 radius of 1 (spots above, right, below, and left)
        offsets <- data.frame(x.offset=c( 0, 1, 0, -1),
                              y.offset=c(-1, 0, 1,  0))
    } else {
        stop(sprintf(".find_neighbors: Unsupported platform %s", platform))
    }
    
    ## Get array coordinates (and label by index of spot in SCE)
    spot.positions <- colData(sce)[, c("col", "row")]
    spot.positions$spot.idx <- seq_len(nrow(spot.positions))
    
    ## Compute coordinates of each possible spot neighbor
    neighbor.positions <- merge(spot.positions, offsets)
    neighbor.positions$x.pos <- neighbor.positions$col + neighbor.positions$x.offset
    neighbor.positions$y.pos <- neighbor.positions$row + neighbor.positions$y.offset
    
    ## Select spots that exist at neighbor coordinates
    neighbors <- merge(as.data.frame(neighbor.positions), 
                       as.data.frame(spot.positions), 
                       by.x=c("x.pos", "y.pos"), by.y=c("col", "row"),
                       suffixes=c(".primary", ".neighbor"))
    
    ## Shift to zero-indexing for C++
    neighbors$spot.idx.neighbor <- neighbors$spot.idx.neighbor - 1
    
    ## Group neighbor indices by spot 
    ## (sort first for consistency with older implementation)
    neighbors <- neighbors[order(neighbors$spot.idx.primary, 
                                 neighbors$spot.idx.neighbor), ]
    df_j <- split(neighbors$spot.idx.neighbor, neighbors$spot.idx.primary)
    
    unname(df_j)
}

#' Initialize cluster assignments
#' 
#' @param sce SingleCellExperiment
#' @param q Number of clusters
#' @param inputs Results from .prepare_inputs() (TODO: store this in sce)
#' @param init Vector of initial cluster assignments
#' @param init.method Initialization clustering algorithm
#' 
#' @return Modified sce with initial cluster assignments stored in
#'   colData$cluster.init (TODO: return vector instead)
#' 
#' @keywords internal
#' 
#' @importFrom stats kmeans
#' @importFrom mclust Mclust mclustBIC
.init_cluster <- function(Y, q, init = NULL, init.method = c("mclust", "kmeans")) {
    if (is.null(init)) {
        init.method <- match.arg(init.method)
        if (init.method == "kmeans") {
            init <- kmeans(Y, centers=q)$cluster
        } else if (init.method == "mclust") {
            init <- Mclust(Y, q, "EEE", verbose=FALSE)$classification
        }
    }
    
    init
}