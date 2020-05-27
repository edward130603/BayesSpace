#deconvolution as described by Ripley 1991
deconvolveRipley = function(Y, positions, nrep = 1000, every = 1, gamma = 2, dist, q, init, seed = 100, mu0 = colMeans(Y), lambda0 = diag(0.01, nrow = ncol(Y)), alpha = 1, beta = 0.01){
  set.seed(seed)
  d = ncol(Y)
  n0 = nrow(Y)
  positions = as.matrix(positions)
  Y = as.matrix(Y)
  colnames(positions) = c("x", "y")
  
  Y = Y[rep(seq_len(n0), 9), ] #rbind 9 times
  positions = positions[rep(seq_len(n0), 9), ] #rbind 9 times
  shift = rbind(expand.grid(c(1/3, -1/3, 0), c(1/3,-1/3,0)))
  shift_long = shift[rep(seq_len(9), each = n0), ]
  positions[,"x"] = positions[,"x"] + shift_long[,"Var1"]
  positions[,"y"] = positions[,"y"] + shift_long[,"Var2"]
  n = nrow(Y)
  df_j = sapply(1:n, function(x){which((abs(positions[,1] -positions[x,1]) + abs(positions[,2] - positions[x,2])) <= dist &  
                                         (abs(positions[,1] -positions[x,1]) + abs(positions[,2] - positions[x,2])) > 0)-1})
  
  iterate2(Y = Y, df_j = df_j, nrep = nrep, n = n, n0 = n0, d = d, gamma = gamma, q = q, init = init, mu0 = mu0, lambda0 = lambda0, alpha = alpha, beta = beta)
}