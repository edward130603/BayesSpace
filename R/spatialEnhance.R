#' Enhance spot resolution
#' 
#' Backend calls iterate_deconv(), written in Rcpp.
#' Inputs are the same as \code{spatialCluster()} except you have to specify
#' xdist and ydist instead of total dist...(maybe would be better to change
#' \code{spatialCluster} to match this)
#' 
#' @param sce A SingleCellExperiment object containing the spatial data.
#' @param q The number of clusters.
#' @param use.dimred Name of a reduced dimensionality result in 
#'   \code{reducedDims(sce)}. If provided, cluster on these features directly. 
#' @param d Number of top principal components to use when clustering.
#' @param positions A matrix or dataframe with two columns (x, y) that gives
#'   the spatial coordinates of each spot.
#' @param position.cols If \code{positions} is not provided, use these columns 
#'   from \code{colData(sce)} as spatial coordinates (x, y).
#' @param init Initial cluster assignments for spots.
#' @param init.method If \code{init} is not provided, cluster the top \code{d} 
#'   PCs with this method to obtain initial cluster assignments.
#' @param xdist The distance along x-axis between neighboring spots. If not 
#'   provided, the distance will be estimated using \code{lm(imagecol ~ col)}.
#' @param ydist The distance along y-axis between neighboring spots. If not
#'   provided, the distance will be estimated using \code{lm(imagerow ~ row)}.
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
#' @param platform Spatial transcriptomic platform. Specify 'Visium' for hex 
#'   lattice geometry or 'ST' for square lattice geometry.
#' @param jitter_scale Controls the amount of jittering. Small amounts of 
#'   jittering are more likely to be accepted but result in exploring the space
#'   more slowly. We suggest tuning \code{jitter_scale} so that Ychange is on 
#'   average around 30\%.
#' @param c Prior precision (1/variance) of the actual spot-level expression 
#'   value in PC space. We suggest making c larger if the jittered values are
#'   not expected to vary much from the overall mean of the spot.
#' @param verbose Log progress to stderr.
#'  
#' @return Returns a new SingleCellExperiment object. By default, the 
#'   \code{assays} of this object are empty, and the enhanced resolution PCs 
#'   are stored as a reduced dimensionality result accessible with
#'   \code{reducedDim(sce, 'PCA')}.
#'   
#' @details TODO: add details on subspot coordinates and linking to parent spots
#' 
#' @name spatialEnhance
NULL

## TODO: re-order arguments so all keyword arguments come after positionals
deconvolve <- function(Y, positions, nrep = 1000, gamma = 2, xdist, ydist, q, 
    init, model = "normal", platform = c("Visium", "ST"), verbose = TRUE, 
    jitter_scale = 5, c = 0.01, mu0 = colMeans(Y), 
    lambda0 = diag(0.01, nrow = ncol(Y)), alpha = 1, beta = 0.01) {
    
    ## TODO: turn progress bar off when verbose=FALSE
    
    d <- ncol(Y)
    n0 <- nrow(Y)
    positions <- as.matrix(positions)
    Y <- as.matrix(Y)
    colnames(positions) <- c("x", "y")
    
    platform <- match.arg(platform)
    subspots <- ifelse(platform == "Visium", 7, 9)  # TODO: parameterize?
    
    init1 <- rep(init, subspots)
    Y2 <- Y[rep(seq_len(n0), subspots), ]  # rbind 7 or 9 times
    positions2 <- positions[rep(seq_len(n0), subspots), ]  # rbind 7 times
    if (platform == "Visium") {
        shift <- rbind(expand.grid(c(1/3, -1/3), c(1/3, -1/3)), 
            expand.grid(c(2/3, -2/3, 0), 0))
    } else {
        shift <- rbind(expand.grid(c(1/3, -1/3, 0), c(1/3, -1/3, 0)))
    }
    
    shift <- t(t(shift) * c(xdist, ydist))
    dist <- max(rowSums(abs(shift))) * 1.05
    if (platform == "ST") {
        dist <- dist/2
    }
    shift_long <- shift[rep(seq_len(subspots), each=n0), ]
    positions2[, "x"] <- positions2[, "x"] + shift_long[, "Var1"]
    positions2[, "y"] <- positions2[, "y"] + shift_long[, "Var2"]
    n <- nrow(Y2)
    
    if (verbose) message("Calculating neighbors...")
    df_j <- find_neighbors(positions2, dist, "manhattan")
    
    if (verbose) message("Fitting model...")
    tdist <- ifelse(model == "t", TRUE, FALSE)
    out <- iterate_deconv(Y=Y2, df_j=df_j, tdist=tdist, nrep=nrep, n=n, n0=n0,
        d=d, gamma=gamma, q=q, init=init1, subspots=subspots, verbose=verbose, 
        jitter_scale=jitter_scale, c=c, mu0=mu0, lambda0=lambda0, alpha=alpha, 
        beta=beta)
    out$positions <- positions2
    
    out
}

#' @export
#' @rdname spatialEnhance
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim<-
#' @importFrom SummarizedExperiment rowData 
spatialEnhance <- function(sce, q, use.dimred = "PCA", d = 15, 
    positions = NULL, position.cols = c("imagecol", "imagerow"), 
    init = NULL, init.method = c("kmeans", "spatialCluster"),
    xdist = NULL, ydist = NULL, 
    model = c("normal", "t"), nrep = 1000, gamma = 2, 
    mu0 = NULL, lambda0 = NULL, alpha = 1, beta = 0.01, 
    save.chain = FALSE, chain.fname = NULL, platform = c("Visium", "ST"), 
    jitter_scale = 5, c = 0.01, verbose = FALSE) {
    
    inputs <- .prepare_inputs(sce, use.dimred=use.dimred, d=d,
        positions=positions, position.cols=position.cols,
        xdist=xdist, ydist=ydist)
    
    ## Initialize cluster assignments (use k-means for now)
    if (is.null(init)) {
        init.method <- match.arg(init.method)
        if (init.method == "kmeans") {
            init <- kmeans(inputs$PCs, centers=q)$cluster
        } else if (init.method == "spatialCluster") {
            ## TODO: check for spatial.cluster in sce, auto-run with same params
            init <- sce$spatial.cluster
        }
    }
    
    ## TODO: pass these through with ...
    model <- match.arg(model)
    platform <- match.arg(platform)
    mu0 <- if (is.null(mu0)) colMeans(inputs$PCs) else mu0
    lambda0 <- if (is.null(lambda0)) diag(0.01, ncol(inputs$PCs)) else lambda0
    
    deconv <- deconvolve(inputs$PCs, inputs$positions, nrep=nrep, gamma=gamma, 
        xdist=inputs$xdist, ydist=inputs$ydist, q=q, init=init, model=model, 
        platform=platform, verbose=verbose, jitter_scale=jitter_scale, c=c, 
        mu0=mu0, lambda0=lambda0, alpha=alpha, beta=beta)
    
    ## Create enhanced SCE
    
    ## TODO: auto-run predictExpression and include as assay 
    ## deconv_PCs <- deconv$Y[[length(deconv$Y)]] 
    ## expr <- predictExpression(sce, deconv_PCs)
    
    cdata <- as.data.frame(deconv$positions)
    colnames(cdata) <- c("imagecol", "imagerow")
    enhanced <- SingleCellExperiment(assays=list(), 
        rowData=rowData(sce), colData=cdata)
    
    reducedDim(enhanced, "PCA") <- deconv$Y[[length(deconv$Y)]]
    ## TODO: fix hard-coding of iterations being used
    enhanced$spatial.cluster <- apply(deconv$z[900:1000, ], 2, Mode)
    
    if (save.chain) {
        deconv <- .clean_chain(deconv, method="enhance")
        params <- c("z", "mu", "lambda", "weights", "Y", "Ychange")
        metadata(enhanced)$chain.h5 <- .write_chain(deconv, chain.fname, params)
    }
    
    enhanced
}
