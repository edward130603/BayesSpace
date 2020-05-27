#cluster() is used to perform spatial clustering
#cluster() calls iterate() which is written in Rcpp

#Y is a matrix or dataframe with 1 row per spot and 1 column per outcome (e.g. principal component)
#positions is a matrix or dataframe with two columns (x, y) that gives the spatial coordinates of the spot
#dist is the L1 distance between neighboring spots
#model is either "normal" or "t"
#precision is either "equal" or "variable" for EEE and VVV covariance models respectively
#q is the number of clusters
#nrep is the number of mcmc iterations
#gamma is the smoothing parameter, values in range of 1-3 seem to work well generally
#init is the initial states (z's), vector of length equal to the number of rows of Y and positions
#mu0 is prior mean hyperparam for mu
#lambda0 is prior precision hyperparam for mu
#alpha,beta are additional hyperparameters Wishart distributed precision lambda

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
