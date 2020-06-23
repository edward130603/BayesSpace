#' Spatial clustering
#'
#' Backend calls iterate() which is written in Rcpp
#' 
#' @param Y A matrix or dataframe with 1 row per spot and 1 column per outcome 
#'   (e.g. principal components).
#' @param positions A matrix or dataframe with two columns (x, y) that gives
#'   the spatial coordinates of each spot.
#' @param neighborhood.radius The maximum (L1) distance for two spots to be
#'   considered neighbors.
#' @param gamma Smoothing parameter. (Values in range of 1-3 seem to work well.)
#' @param q The number of clusters.
#' @param init Initial cluster assignments (z's). Must be a vector of length
#'   equal to the number of rows of Y and positions.
#' @param model Error model. ("normal" or "t")
#' @param precision Covariance structure. ("equal" or "variable" for EEE and 
#'   VVV covariance models, respectively.)
#' @param nrep The maximum number of MCMC iterations.
#' @param mu0 Prior mean hyperparameter for mu.
#' @param lambda0 Prior precision hyperparam for mu.
#' @param alpha Hyperparameter for Wishart distributed precision lambda.
#' @param beta Hyperparameter for Wishart distributed precision lambda.
#' 
#' @return List of parameter values (`z`, `mu`, `lambda`) and model 
#'         log-likelihoods (`plogLik`) at each MCMC iteration, along with final
#'         cluster labels (`labels`)
#'         
#' @details TODO describe method in detail
cluster = function(Y, positions, neighborhood.radius, q, 
                   model = c("normal", "t"), precision = c("equal", "variable"),
                   z0 = rep(1, nrow(Y)), mu0 = colMeans(Y), lambda0 = diag(0.01, nrow = ncol(Y)),
                   gamma = 2, alpha = 1, beta = 0.01, nrep = 1000) {
  
  positions = as.matrix(positions)
  Y = as.matrix(Y)
  d = ncol(Y) 
  n = nrow(Y)
  
  if ((nrow(positions) != n) | (length(z0) != n)){
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
  
  out <- cluster.FUN(Y = as.matrix(Y), df_j = df_j, nrep = nrep, n = n, d = d, 
                     gamma = gamma, q = q, init = z0, mu0 = mu0, 
                     lambda0 = lambda0, alpha = alpha, beta = beta)
  
  iter_from <- ifelse(nrep < 2000, max(2, nrep - 1000), 1000)
  message("Calculating labels using iterations ", iter_from, " through ", nrep, "...")
  out$labels <- apply(out$z[iter_from:nrep,], 2, Mode)
  
  out
}


#' @importFrom stats kmeans
#' @importFrom SingleCellExperiment reducedDim<-
SpatialCluster <- function(sce, q,
                           init=NULL, init.method=c("kmeans"),
                           positions=NULL, position.cols=c("row", "col"),
                           neighborhood.radius=NULL,
                           assay.type="logcounts") {
  
  # TODO: 
  #   1. Allow specification of reduced dimension (and number of top features)
  #   2. Default to PCA, running PCA if necessary
  # Y <- assay(sce, assay.type)
  # Y <- as.matrix(reducedDim(sce, dim.name)
  
  if (is.null(positions)) {
    positions <- cbind(sce[[x_col]], sce[[y_col]])
    colnames(positions) <- c("x", "y")
  }
  
  if (is.null(neighborhood.radius)) {
    neighborhood.radius <- compute_neighborhood_radius(sce)
  }
  
  # Initialize cluster assignments (use k-means for now)
  if (is.null(init)) {
    init.method <- match.arg(init.method)
    if (init.method == "kmeans") {
      init <- kmeans(PCs$components, centers = q)$cluster
    }
  }
  
  Y <- reducedDim(sce, "PCA")
}