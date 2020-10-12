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
#' @param burn.in The number of MCMC iterations to exclude as burn-in period.
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
#' The input SCE must have \code{row} and \code{col} columns in its
#' \code{colData}, corresponding to the array row and column coordinates of each
#' spot. These are automatically parsed by \code{\link{readVisium}} or can be
#' added manually when creating the SCE.
#' 
#' Cluster labels are stored in the \code{spatial.cluster} column of the SCE,
#' and the cluster initialization is stored in \code{cluster.init}.
#' 
#' @examples
#' set.seed(149)
#' sce <- exampleSCE()
#' sce <- spatialCluster(sce, 7, nrep=200, burn.in=20)
#' 
#' @seealso \code{\link{spatialPreprocess}} for preparing the SCE for
#'   clustering, \code{\link{spatialEnhance}} for enhancing the clustering
#'   resolution, \code{\link{clusterPlot}} for visualizing the cluster
#'   assignments, \code{\link{featurePlot}} for visualizing expression levels
#'   in spatial context, and \code{\link{mcmcChain}} for examining the full
#'   MCMC chain associated with the clustering.
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
    
    if (length(mu0) != d)
        stop("Dimensions of mu0 do not match input data Y.")
    if (ncol(lambda0) != d)
        stop("Dimensions of lambda0 do not match input data Y.")
    
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

#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom assertthat assert_that
#' 
#' @export
#' @rdname spatialCluster
spatialCluster <- function(sce, q, use.dimred = "PCA", d = 15,
    platform=c("Visium", "ST"),
    init = NULL, init.method = c("mclust", "kmeans"),
    model = c("t", "normal"), precision = c("equal", "variable"), 
    nrep = 50000, burn.in=1000, gamma = 3, mu0 = NULL, lambda0 = NULL,
    alpha = 1, beta = 0.01, save.chain = FALSE, chain.fname = NULL) {
    
    if (!(use.dimred %in% reducedDimNames(sce))) 
        stop("reducedDim \"", use.dimred, "\" not found in input SCE.")

    ## Require at least one iteration and non-negative burn-in
    assert_that(nrep >= 1)
    assert_that(burn.in >= 0)
    if (burn.in >= nrep)
        stop("Please specify a burn-in period shorter than the total number of iterations.")

    ## Get PCs
    Y <- reducedDim(sce, use.dimred)
    d <- min(ncol(Y), d)
    Y <- Y[, seq_len(d)]
    
    ## Get indices of neighboring spots, and initialize cluster assignments
    ## TODO: parse platform from metadata
    platform <- match.arg(platform)
    df_j <- .find_neighbors(sce, platform)
    init <- .init_cluster(Y, q, init, init.method)
    
    ## Set model parameters
    model <- match.arg(model)
    precision <- match.arg(precision)
    if (is.null(mu0))
        mu0 <- colMeans(Y)
    if (is.null(lambda0))
        lambda0 <- diag(0.01, ncol(Y))
    
    ## Run clustering
    ## TODO: set default gamma to 2 if platform=ST and not specified
    results <- cluster(Y, q, df_j, init=init, 
        model=model, precision=precision, mu0=mu0, 
        lambda0=lambda0, gamma=gamma, alpha=alpha, beta=beta, nrep=nrep)
    
    ## Save MCMC chain
    if (save.chain) {
        results <- .clean_chain(results)
        metadata(sce)$chain.h5 <- .write_chain(results, chain.fname)
    }
    
    ## Save metadata
    sce$cluster.init <- init
    if (!exists("BayesSpace.data", metadata(sce)))
        metadata(sce)$BayesSpace.data <- list()
    metadata(sce)$BayesSpace.data$platform <- platform
    metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
    
    ## Save modal cluster assignments, excluding burn-in
    message("Calculating labels using iterations ", burn.in, " through ", nrep, ".")
    zs <- results$z[seq(burn.in + 1, nrep), ]
    if (burn.in + 1 == nrep)
        labels <- matrix(zs, nrow=1)  # if only one iteration kept, return it
    else
        labels <- apply(zs, 2, Mode)  # else take modal assignment
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
#' @importFrom purrr keep discard map
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
        stop(".find_neighbors: Unsupported platform \"", platform, "\".")
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
                       suffixes=c(".primary", ".neighbor"),
                       all.x=TRUE)
    
    ## Shift to zero-indexing for C++
    neighbors$spot.idx.neighbor <- neighbors$spot.idx.neighbor - 1
    
    ## Group neighbor indices by spot 
    ## (sort first for consistency with older implementation)
    neighbors <- neighbors[order(neighbors$spot.idx.primary, 
                                 neighbors$spot.idx.neighbor), ]
    df_j <- split(neighbors$spot.idx.neighbor, neighbors$spot.idx.primary)
    df_j <- unname(df_j)
    
    ## Discard neighboring spots without spot data
    ## This can be implemented by eliminating `all.x=TRUE` above, but
    ## this makes it easier to keep empty lists for spots with no neighbors
    ## (as expected by C++ code)
    df_j <- map(df_j, function(nbrs) discard(nbrs, function(x) is.na(x)))
    
    ## Log number of spots with neighbors
    n_with_neighbors <- length(keep(df_j, function(nbrs) length(nbrs) > 0))
    message("Neighbors were identified for ", n_with_neighbors, " out of ",
            ncol(sce), " spots.")
    
    df_j
}

#' Initialize cluster assignments
#' 
#' @param sce SingleCellExperiment
#' @param q Number of clusters
#' @param inputs Results from \code{.prepare_inputs()}
#' @param init Vector of initial cluster assignments
#' @param init.method Initialization clustering algorithm
#' 
#' @return Vector of cluster assignments.
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
