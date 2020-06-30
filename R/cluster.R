#' Spatial clustering
#'
#' Cluster a spatial expression dataset.
#' 
#' @param sce A SingleCellExperiment object containing the spatial data.
#' @param q The number of clusters.
#' @param use.dimred Name of a reduced dimensionality result in 
#'   \code{reducedDims(sce)}. If provided, cluster on these features directly. 
#' @param pca.method If \code{use.dimred} is not provided, run PCA on the 
#'   specified assay. Can run standard or de-noised PCA.
#' @param d Number of top principal components to use when clustering.
#' @param assay.type If \code{use.dimred} is not provided, use this assay when 
#'   running PCA.
#' @param positions A matrix or dataframe with two columns (x, y) that gives
#'   the spatial coordinates of each spot.
#' @param position.cols If \code{positions} is not provided, use these columns 
#'   from \code{colData(sce)} as spatial coordinates (x, y).
#' @param init Initial cluster assignments for spots.
#' @param init.method If \code{init} is not provided, cluster the top \code{d} 
#'   PCs with this method to obtain initial cluster assignments.
#' @param neighborhood.radius The maximum (L1) distance for two spots to be
#'   considered neighbors. If not provided, the distance will be estimated
#'   using \code{lm(image.coord ~ array.coord)}.
#' @param model Error model. ("normal" or "t")
#' @param precision Covariance structure. ("equal" or "variable" for EEE and 
#'   VVV covariance models, respectively.)
#' @param nrep The number of MCMC iterations.
#' @param gamma Smoothing parameter. (Values in range of 1-3 seem to work well.)
#' @param mu0 Prior mean hyperparameter for mu.
#' @param lambda0 Prior precision hyperparam for mu.
#' @param alpha Hyperparameter for Wishart distributed precision lambda.
#' @param beta Hyperparameter for Wishart distributed precision lambda.
#' @param model Error model. ("normal" or "t")
#' @param precision Covariance structure. ("equal" or "variable" for EEE and 
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
#'   TODO describe method in detail
#'   TODO add details or link to mcmcChain
#'
#' @name spatialCluster
NULL

#' @importFrom purrr map
cluster = function(Y, positions, neighborhood.radius, q, 
                   model = c("normal", "t"), precision = c("equal", "variable"),
                   init = rep(1, nrow(Y)), mu0 = colMeans(Y), lambda0 = diag(0.01, nrow = ncol(Y)),
                   gamma = 2, alpha = 1, beta = 0.01, nrep = 1000) 
{
  
  positions = as.matrix(positions)
  Y = as.matrix(Y)
  d = ncol(Y) 
  n = nrow(Y)
  
  if ((nrow(positions) != n) | (length(init) != n)){
    stop("Dimensions of Y, positions, and init do not match")
  }
  if ((length(mu0) != d) | (ncol(lambda0) != d)) {
    stop("Dimensions of mu0 or lambda0 do not match input data Y")
  }
  
  model <- match.arg(model)
  precision <- match.arg(precision)
  
  if (q==1){
    return(list(z = matrix(rep(1, n), nrow =1)))
  }
  
  # TODO: pass boolean matrix to cpp instead of using sapply?
  df_j <- find_neighbors(positions, neighborhood.radius, "manhattan")
  
  message("Fitting model...")
  if (model == "normal"){
    if (precision == "equal"){
      cluster.FUN <- iterate
    } else if (precision == "variable"){
      cluster.FUN <- iterate_vvv
    } 
  } else if (model == "t"){
    if (precision == "equal"){
      cluster.FUN <- iterate_t
    } else if (precision == "variable"){
      cluster.FUN <- iterate_t_vvv
    }
  } 
  
  cluster.FUN(Y = as.matrix(Y), df_j = df_j, nrep = nrep, n = n, d = d, 
              gamma = gamma, q = q, init = init, mu0 = mu0, lambda0 = lambda0, 
              alpha = alpha, beta = beta)
}

# TODO make generic
#' @importFrom stats kmeans
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors metadata metadata<-
#' @rdname spatialCluster
spatialCluster <- function(sce, q, use.dimred=NULL, 
                           pca.method=c("PCA", "denoised", "geneset"), d=15,
                           assay.type="logcounts",
                           positions=NULL, position.cols=c("imagecol", "imagerow"),
                           init=NULL, init.method=c("kmeans"),
                           neighborhood.radius=NULL,
                           model=c("normal", "t"), precision=c("equal", "variable"),
                           nrep=1000, gamma=2, mu0=NULL, lambda0=NULL, 
                           alpha=1, beta=0.01,
                           save.chain=FALSE, chain.fname=NULL) {

  if (is.null(use.dimred)) {
    use.dimred <- "PCA"
    pca.method <- match.arg(pca.method)
    sce <- addPCA(sce, assay.type, pca.method, d)
  }
  
  PCs <- reducedDim(sce, use.dimred)
  d <- min(d, ncol(PCs))
  PCs <- PCs[, 1:d]
  
  if (is.null(positions)) {
    positions <- as.matrix(colData(sce)[position.cols])
    colnames(positions) <- c("x", "y")
  }
  
  if (is.null(neighborhood.radius)) {
    neighborhood.radius <- compute_neighborhood_radius(sce)
  }
  
  # Initialize cluster assignments (use k-means for now)
  if (is.null(init)) {
    init.method <- match.arg(init.method)
    if (init.method == "kmeans") {
      init <- kmeans(PCs, centers = q)$cluster
    }
  }
  
  # TODO: pass these through with ...
  model <- match.arg(model)
  precision <- match.arg(precision)
  mu0 <- if (is.null(mu0)) colMeans(PCs) else mu0
  lambda0 <- if (is.null(lambda0)) diag(0.01, ncol(PCs)) else lambda0
  
  results <- cluster(PCs, positions, neighborhood.radius, q, 
                     init = init,
                     model = model, precision = precision,
                     mu0 = mu0, lambda0 = lambda0,
                     gamma = gamma, alpha = alpha, beta = beta, nrep = nrep) 
    
  if (save.chain) {
    results <- .clean_chain(results)
    metadata(sce)$chain.h5 <- .write_chain(results, chain.fname)
  }
  
  iter_from <- ifelse(nrep < 2000, max(2, nrep - 1000), 1000)
  # message("Calculating labels using iterations ", iter_from, " through ", nrep, "...")f
  results$labels <- apply(results$z[iter_from:nrep,], 2, Mode)
  
  # TODO: switch to labels computed above (this is just for sanity test)
  colData(sce)$spatial.cluster <- apply(results$z[900:1000, ], 2, Mode)
  
  sce
}
