#deconvolve
#deconvolve calls iterate_deconv(), written in Rcpp
#inputs are the same as cluster() except you have to specify xdist and ydist instead of total dist...(maybe would be better to change cluster() to match this)
deconvolve = function(Y, positions, nrep = 1000, every = 1, gamma = 2, xdist, ydist, q, init, seed = 100, mu0 = colMeans(Y), lambda0 = diag(0.01, nrow = ncol(Y)), alpha = 1, beta = 0.01){
  set.seed(seed)
  
  d = ncol(Y)
  n0 = nrow(Y)
  positions = as.matrix(positions)
  Y = as.matrix(Y)
  colnames(positions) = c("x", "y")
  
  init1 = rep(init, 7)
  Y2 = Y[rep(seq_len(n0), 7), ] #rbind 7 times
  positions2 = positions[rep(seq_len(n0), 7), ] #rbind 7 times
  shift = rbind(expand.grid(c(1/3, -1/3), c(1/3,-1/3)), expand.grid(c(2/3, -2/3,0), 0))
  shift = t(t(shift)*c(xdist, ydist))
  dist = max(rowSums(abs(shift)))*1.05
  shift_long = shift[rep(seq_len(7), each = n0), ]
  positions2[,"x"] = positions2[,"x"] + shift_long[,"Var1"]
  positions2[,"y"] = positions2[,"y"] + shift_long[,"Var2"]
  n = nrow(Y2)
  print("Calculating neighbors...")
  df_j = sapply(1:n, function(x){which((abs(positions2[,1] -positions2[x,1]) + abs(positions2[,2] - positions2[x,2])) <= dist &  #L1 distance
                                         (abs(positions2[,1] -positions2[x,1]) + abs(positions2[,2] - positions2[x,2])) > 0)-1})
  print("Fitting model...")
  iterate_deconv(Y = Y2, df_j = df_j, nrep = nrep, n = n, n0 = n0, d = d, gamma = gamma, q = q, init = init1, mu0 = mu0, lambda0 = lambda0, alpha = alpha, beta = beta)
}