#' Spatial clustering
#'
#' Backend calls iterate() which is written in Rcpp
#' 
#' @param Y A matrix or dataframe with 1 row per spot and 1 column per outcome (e.g. principal component)
#' @param positions A matrix or dataframe with two columns (x, y) that gives the spatial coordinates of the spot
#' @param dist The L1 distance between neighboring spots
#' @param gamma Smoothing parameter. Values in range of 1-3 seem to work well generally
#' @param q The number of clusters
#' @param init Initial cluster assignments (z's). Must be a vector of length equal to the number of rows of Y and positions
#' @param model Error model ("normal" or "t")
#' @param precision Covariance structure ("equal" or "variable" for EEE and VVV covariance models, respectively)
#' @param nrep The maximum number of mcmc iterations
#' @param seed Random state seed
#' @param mu0 Prior mean hyperparameter for mu
#' @param lambda0 Prior precision hyperparam for mu
#' @param alpha Hyperparameter for Wishart distributed precision lambda
#' @param beta Hyperparameter for Wishart distributed precision lambda
#' 
#' @return List of parameter values (`z`, `mu`, `lambda`) and model log-likelihoods (`plogLik`) at each MCMC iteration, along with final cluster labels (`labels`)
cluster = function(Y, positions, dist, gamma = 2, q, init = rep(1, nrow(Y)), model = "normal", precision = "equal", nrep = 1000, seed = 100, mu0 = colMeans(Y), lambda0 = diag(0.01, nrow = ncol(Y)), alpha = 1, beta = 0.01){
  
  set.seed(seed)
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
  if (q==1){
    return(list(z = matrix(rep(1, n), nrow =1)))
  }
  colnames(positions) = c("x", "y")
  message("Calculating neighbors...")
  df_j = sapply(1:n, function(x){which((abs(positions[,1] -positions[x,1]) + abs(positions[,2] - positions[x,2])) <= dist &  #L1 distance
                                         (abs(positions[,1] -positions[x,1]) + abs(positions[,2] - positions[x,2])) > 0)-1})
  num_neighbors = sapply(df_j, length)
  num_message = paste0("Neighbors were identified for ", sum(num_neighbors >0), " out of ", n, " spots.")
  message(num_message)
  
  message("Fitting model...")
  if (model == "normal"){
    if (precision == "equal"){
      out = iterate(Y = as.matrix(Y), df_j = df_j, nrep = nrep, n = n, d = d, gamma = gamma, q = q, init = init, mu0 = mu0, lambda0 = lambda0, alpha = alpha, beta = beta)
    } else if (precision == "variable"){
      out = iterate_vvv(Y = as.matrix(Y), df_j = df_j, nrep = nrep, n = n, d = d, gamma = gamma, q = q, init = init, mu0 = mu0, lambda0 = lambda0, alpha = alpha, beta = beta)
    } else {
      stop("precision should be either 'equal' or 'variable'")
    }
  } else if (model == "t"){
    if (precision == "equal"){
      out = iterate_t(Y = as.matrix(Y), df_j = df_j, nrep = nrep, n = n, d = d, gamma = gamma, q = q, init = init, mu0 = mu0, lambda0 = lambda0, alpha = alpha, beta = beta)
    } else if (precision == "variable"){
      out = iterate_t_vvv(Y = as.matrix(Y), df_j = df_j, nrep = nrep, n = n, d = d, gamma = gamma, q = q, init = init, mu0 = mu0, lambda0 = lambda0, alpha = alpha, beta = beta)
    } else {
      stop("precision should be either 'equal' or 'variable'")
    }
  } else {
    stop("model should be either 'normal' or 't'")
  }
  out$labels = apply(out$z[max(nrep-1000, 2):nrep,], 2, Mode)
  out
}
