#' Enhance spot resolution
#' 
#' Backend calls iterate_deconv(), written in Rcpp.
#' Inputs are the same as \code{spatialCluster()} except you have to specify
#' xdist and ydist instead of total dist...(maybe would be better to change
#' \code{spatialCluster} to match this)
#' 
#' @param sce A SingleCellExperiment object containing the spatial data.
#' @param q The number of clusters.
#' @param platform Spatial transcriptomic platform. Specify 'Visium' for hex 
#'   lattice geometry or 'ST' for square lattice geometry. Specifying this
#'   parameter is optional when analyzing SingleCellExperiments processed using
#'   \code{\link{readVisium}}, \code{\link{spatialPreprocess}}, or
#'   \code{\link{spatialCluster}}, as this information is included in their
#'   metadata.
#' @param use.dimred Name of a reduced dimensionality result in 
#'   \code{reducedDims(sce)}. If provided, cluster on these features directly. 
#' @param d Number of top principal components to use when clustering.
#' @param init Initial cluster assignments for spots.
#' @param init.method If \code{init} is not provided, cluster the top \code{d} 
#'   PCs with this method to obtain initial cluster assignments.
#' @param model Error model. ('normal' or 't')
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
#' @param burn.in Number of iterations to exclude as burn-in period. The MCMC
#'   iterations are currently thinned to every 100; accordingly \code{burn.in}
#'   is rounded down to the nearest multiple of 100.
#' @param jitter_scale Controls the amount of jittering. Small amounts of 
#'   jittering are more likely to be accepted but result in exploring the space
#'   more slowly. We suggest tuning \code{jitter_scale} so that Ychange is on 
#'   average around 30\%.
#' @param jitter_prior Scale factor for the prior variance, parameterized as the
#'   proportion (default = 0.3) of the mean variance of the PCs.
#'   We suggest making \code{jitter_prior} smaller if the jittered values are
#'   not expected to vary much from the overall mean of the spot.
#' @param verbose Log progress to stderr.
#'  
#' @return Returns a new SingleCellExperiment object. By default, the 
#'   \code{assays} of this object are empty, and the enhanced resolution PCs 
#'   are stored as a reduced dimensionality result accessible with
#'   \code{reducedDim(sce, 'PCA')}.
#'   
#' @details 
#' The enhanced \code{SingleCellExperiment} has most of the properties of the
#'   input SCE - \code{rowData}, \code{colData}, \code{reducedDims} - but does
#'   not include expression data in \code{counts} or \code{logcounts}. To impute
#'   enhanced expression vectors, please use [enhanceFeatures()] after
#'   running \code{spatialEnhance}.
#'   
#' The \code{colData} of the enhanced \code{SingleCellExperiment} includes the
#'   following columns to permit referencing the subspots in spatial context and
#'   linking back to the original spots:
#'   \itemize{
#'   \item \code{spot.idx}: Index of the spot this subspot belongs to (with
#'     respect to the input SCE).
#'   \item \code{subspot.idx}: Index of the subspot within its parent spot.
#'   \item \code{spot.row}: Array row of the subspot's parent spot.
#'   \item \code{spot.col}: Array col of the subspot's parent spot.
#'   \item \code{row}: Array row of the subspot. This is the parent spot's row
#'     plus an offset based on the subspot's position within the spot.
#'   \item \code{col}: Array col of the subspot. This is the parent spot's col
#'     plus an offset based on the subspot's position within the spot.
#'   \item \code{imagerow}: Pixel row of the subspot. This is the parent spot's
#'     row plus an offset based on the subspot's position within the spot.
#'   \item \code{imagecol}: Pixel col of the subspot. This is the parent spot's
#'     col plus an offset based on the subspot's position within the spot.
#'   }
#' 
#' @examples
#' set.seed(149)
#' sce <- exampleSCE()
#' sce <- spatialCluster(sce, 7, nrep=200, burn.in=20)
#' enhanced <- spatialEnhance(sce, 7, init=sce$spatial.cluster, nrep=200, burn.in=20)
#' 
#' @seealso \code{\link{spatialCluster}} for clustering at the spot level
#'   before enhancing, \code{\link{clusterPlot}} for visualizing the cluster
#'   assignments, \code{\link{enhanceFeatures}} for imputing enhanced
#'   expression, and \code{\link{mcmcChain}} for examining the full MCMC chain
#'   associated with the enhanced clustering.
#'   .
#' 
#' @name spatialEnhance
NULL

#' Wrapper around C++ \code{iterate_deconv()} function
#' 
#' @return List of enhancement parameter values at each iteration
#' 
#' @keywords internal
#' @importFrom stats cov
deconvolve <- function(Y, positions, xdist, ydist, q, init, nrep = 1000,
    model = "normal", platform = c("Visium", "ST"), verbose = TRUE,
    jitter_scale = 5, jitter_prior = 0.01, mu0 = colMeans(Y), gamma = 2,
    lambda0 = diag(0.01, nrow = ncol(Y)), alpha = 1, beta = 0.01) {
    
    d <- ncol(Y)
    n0 <- nrow(Y)
    Y <- as.matrix(Y)
    c <- jitter_prior * 1 / (2 * mean(diag(cov(Y))))
    
    positions <- as.matrix(positions)
    colnames(positions) <- c("x", "y")
    
    ## If user didn't specify a platform, attempt to parse from SCE metadata
    ## otherwise check against valid options
    if (length(platform) > 1) {
        platform <- .bsData(sce, "platform", match.arg(platform))
    } else {
        platform <- match.arg(platform)
    }

    subspots <- ifelse(platform == "Visium", 6, 9)
    
    init1 <- rep(init, subspots)
    Y2 <- Y[rep(seq_len(n0), subspots), ]  # rbind 6 or 9 times
    positions2 <- positions[rep(seq_len(n0), subspots), ]  # rbind 7 times
    
    shift <- .make_subspot_offsets(subspots)
    shift <- t(t(shift) * c(xdist, ydist))
    dist <- max(rowSums(abs(shift))) * 1.05
    if (platform == "ST") {
        dist <- dist/2
    }
    shift_long <- shift[rep(seq_len(subspots), each=n0), ]
    positions2[, "x"] <- positions2[, "x"] + shift_long[, "Var1"]
    positions2[, "y"] <- positions2[, "y"] + shift_long[, "Var2"]
    n <- nrow(Y2)
    
    if (verbose)
        message("Calculating neighbors...")
    df_j <- find_neighbors(positions2, dist, "manhattan")
    
    if (verbose)
        message("Fitting model...")
    tdist <- (model == "t")
    out <- iterate_deconv(Y=Y2, df_j=df_j, tdist=tdist, nrep=nrep, n=n, n0=n0,
        d=d, gamma=gamma, q=q, init=init1, subspots=subspots, verbose=verbose, 
        jitter_scale=jitter_scale, c=c, mu0=mu0, lambda0=lambda0, alpha=alpha, 
        beta=beta)
    out$positions <- positions2
    out
}

#' Define offsets for each subspot layout.
#' 
#' Hex spots are divided into 6 triangular subspots, square spots are divided
#' into 9 squares. Offsets are relative to the spot center.
#' 
#' @param n_subspots_per Number of subspots per spot
#' @return Matrix of x and y offsets, one row per subspot
#' 
#' @keywords internal
.make_subspot_offsets <- function(n_subspots_per) {
    if (n_subspots_per == 6) {
        rbind(expand.grid(c(1/3, -1/3), c(1/3,-1/3)), expand.grid(c(2/3, -2/3), 0))
    # } else if (n_subspots_per == 7) {
    #     rbind(expand.grid(c(1/3, -1/3), c(1/3, -1/3)), expand.grid(c(2/3, -2/3, 0), 0))
    } else if (n_subspots_per == 9) {
        rbind(expand.grid(c(1/3, -1/3, 0), c(1/3, -1/3, 0)))
    } else {
        stop("Only 6 and 9 subspots currently supported.")
    }
}

#' Add subspot labels and offset row/col locations before making enhanced SCE.
#'
#' Subspots are stored as (1.1, 2.1, 3.1, ..., 1.2, 2.2, 3.2, ...)
#'
#' @param cdata Table of colData (imagerow and imagecol; from deconv$positions)
#' @param sce Original sce (to obtain number of spots and original row/col)
#' @param n_subspots_per Number of subspots per spot
#' 
#' @return Data frame with added subspot names, parent spot indices, and offset
#'   row/column coordinates
#'
#' @keywords internal
#' @importFrom assertthat assert_that
.make_subspot_coldata <- function(positions, sce, n_subspots_per) {
    cdata <- as.data.frame(positions)
    colnames(cdata) <- c("imagecol", "imagerow")
    
    n_spots <- ncol(sce)
    n_subspots <- nrow(cdata)
    assert_that(nrow(cdata) == n_spots * n_subspots_per)
    
    ## Index of parent spot is (subspot % n_spots)
    idxs <- seq_len(n_subspots)
    spot_idxs <- ((idxs - 1) %% n_spots) + 1
    subspot_idxs <- rep(seq_len(n_subspots_per), each=n_spots)
    cdata$spot.idx <- spot_idxs
    cdata$subspot.idx <- subspot_idxs
    rownames(cdata) <- paste0("subspot_", spot_idxs, ".", subspot_idxs)
    
    offsets <- .make_subspot_offsets(n_subspots_per)
    cdata$spot.row <- rep(sce$row, n_subspots_per)
    cdata$spot.col <- rep(sce$col, n_subspots_per)
    cdata$col <- cdata$spot.col + rep(offsets[, 1], each=n_spots)
    cdata$row <- cdata$spot.row + rep(offsets[, 2], each=n_spots)

    cols <- c("spot.idx", "subspot.idx", "spot.row", "spot.col", "row", "col", "imagerow", "imagecol")
    cdata[, cols]
}

#' @export
#' @rdname spatialEnhance
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim<-
#' @importFrom SummarizedExperiment rowData 
#' @importFrom assertthat assert_that
spatialEnhance <- function(sce, q, platform = c("Visium", "ST"),
    use.dimred = "PCA", d = 15,
    init = NULL, init.method = c("spatialCluster", "mclust", "kmeans"),
    model = c("t", "normal"), nrep = 200000, gamma = 3, 
    mu0 = NULL, lambda0 = NULL, alpha = 1, beta = 0.01, 
    save.chain = FALSE, chain.fname = NULL, burn.in=10000,
    jitter_scale = 5, jitter_prior = 0.3, verbose = FALSE) {

    assert_that(nrep >= 100)  # require at least one iteration after thinning
    assert_that(burn.in >= 0)
    if (burn.in >= nrep)
        stop("Please specify a burn-in period shorter than the total number of iterations.")
    
    ## Thinning interval; only every 100 iterations are kept to reduce memory
    ## This is temporarily hard-coded into the C++ code
    thin <- 100

    platform <- match.arg(platform)
    if (platform == "Visium") {
        position.cols <- c("imagecol", "imagerow")
        xdist <- ydist <- NULL  # Compute with .prepare_inputs
    } else if (platform == "ST") {
        position.cols <- c("col", "row")
        xdist <- ydist <- 1
    }
    
    inputs <- .prepare_inputs(sce, use.dimred=use.dimred, d=d,
        positions=NULL, position.cols=position.cols,
        xdist=xdist, ydist=ydist)

    ## Initialize cluster assignments (use spatialCluster by default)
    if (is.null(init)) {
        init.method <- match.arg(init.method)
        if (init.method == "spatialCluster") {
            msg <- paste0("Must run spatialCluster on sce before enhancement ",
                          "if using spatialCluster to initialize.")
            assert_that("spatial.cluster" %in% colnames(colData(sce)), msg=msg)
            init <- sce$spatial.cluster
        } else {
            init <- .init_cluster(inputs$PCs, q, init, init.method)
        }
    }
    
    ## Set model parameters
    model <- match.arg(model)
    if (is.null(mu0))
        mu0 <- colMeans(inputs$PCs)
    if (is.null(lambda0))
        lambda0 <- diag(0.01, ncol(inputs$PCs))
    
    deconv <- deconvolve(inputs$PCs, inputs$positions, nrep=nrep, gamma=gamma, 
        xdist=inputs$xdist, ydist=inputs$ydist, q=q, init=init, model=model, 
        platform=platform, verbose=verbose, jitter_scale=jitter_scale,
        jitter_prior=jitter_prior, mu0=mu0, lambda0=lambda0, alpha=alpha,
        beta=beta)
    
    ## Create enhanced SCE
    n_subspots_per <- ifelse(platform == "Visium", 6, 9)
    cdata <- .make_subspot_coldata(deconv$positions, sce, n_subspots_per)
    enhanced <- SingleCellExperiment(assays=list(), 
        rowData=rowData(sce), colData=cdata)
    
    ## Scale burn.in period to thinned intervals, and
    ## add one to skip initialization values stored before first iteration
    burn.in <- (burn.in %/% thin) + 1
    
    ## Average PCs, excluding burn-in
    deconv_PCs <- Reduce(`+`, deconv$Y[-seq_len(burn.in)]) / (length(deconv$Y) - burn.in)
    colnames(deconv_PCs) <- paste0("PC", seq_len(ncol(deconv_PCs)))
    reducedDim(enhanced, "PCA") <- deconv_PCs
    
    ## Choose modal cluster label, excluding burn-in
    message("Calculating labels using iterations ", (burn.in - 1) * thin,
            " through ", nrep, ".")
    zs <- deconv$z[seq(burn.in, (nrep %/% thin) + 1), ]
    if (burn.in == (nrep %/% thin) + 1)
        labels <- matrix(zs, nrow=1)
    else
        labels <- apply(zs, 2, Mode)

    enhanced$spatial.cluster <- unname(labels)
    
    if (save.chain) {
        deconv <- .clean_chain(deconv, method="enhance")
        params <- c("z", "mu", "lambda", "weights", "Y", "Ychange")
        metadata(enhanced)$chain.h5 <- .write_chain(deconv, chain.fname, params)
    }
    
    ## Add metadata to new SingleCellExperiment object
    metadata(enhanced)$BayesSpace.data <- list()
    metadata(enhanced)$BayesSpace.data$platform <- platform
    metadata(enhanced)$BayesSpace.data$is.enhanced <- TRUE

    enhanced
}
