library(mvtnorm)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

run_mcmc_multi = function(df, nrep = 1000, q = 3, d = 2, mu0 = colMeans(df[,grep("Y",colnames(df))]), lambda0 = diag(0.01, nrow = d), alpha = 1, beta = 0.01, gamma = 2, seed = 100){
  set.seed(seed)
  n = nrow(df)
  df[,"j"] = 1:n
  df_j = sapply(1:n, function(x){df[(abs(df[,"x"] -df[x,"x"]) + abs(df[,"y"] - df[x,"y"])) <= 220 & #new coordinates
                                      (abs(df[,"x"] -df[x,"x"]) + abs(df[,"y"] - df[x,"y"])) > 0,"j"]})
  #df_j = sapply(1:n, function(x){df[(abs(df[,"x"] -df[x,"x"]) + abs(df[,"y"] - df[x,"y"])) == 2,"j"]}) #old coordinates
  
  #Initialize matrices storing iterations
  df_sim_z = matrix(NA, nrow = nrep, ncol = n)
  colnames(df_sim_z) = c(paste0("z",1:n))
  df_sim_mu = matrix(NA, nrow = nrep, ncol = q*d)
  qd = expand.grid(d = 1:d, q = 1:q)
  colnames(df_sim_mu) = sapply(1:(q*d), function(x){paste0("muq", qd[x,2],"_d",qd[x,1])})
  df_sim_lambda = replicate(nrep, matrix(NA, nrow = 2, ncol = 2), simplify = F)

  lambdamu0 = lambda0*mu0
  alpha_n = alpha + n/2
  
  #Initialize parameters
  for (pc in 1:d){
    df_sim_mu[1, grep(paste0("d", pc), colnames(df_sim_mu))] = mu0[pc]
  }
  df_sim_lambda[[1]] = solve(cov(df[,grep("Y", colnames(df))]))
  df_sim_z[1,] = 1 #initialize to all 1
  #df_sim_z[1,] = df$z #initialize to all truth (testing only)
  #df_sim_z[1,] = sample(3, n, replace = T) #initialize to random
  
  Y = as.matrix(df[,grep("Y", colnames(df))])
  #Iterate
  for (i in 2:nrep){
    #mu
    mu_i = matrix(nrow = q, ncol = d)
    index_1 = matrix(nrow = q, ncol = n)
    lambda_prev = df_sim_lambda[[i-1]]
    for(k in 1:q){
      index_1[k,] = df_sim_z[i-1,]==k
      n_i = sum(index_1[k,])
      mean_i = solve(lambda0 + n_i * lambda_prev) %*% (lambda0 %*% mu0 + lambda_prev %*% colSums(Y[index_1[k,],, drop = F]))
      var_i = solve(lambda0 + n_i * lambda_prev)
      mu_i[k,] = rmvnorm(1, mean = mean_i, sigma = var_i)
    }
    df_sim_mu[i,] = as.vector(t(mu_i))
    
    #lambda
    mu_i_long = mu_i[df_sim_z[i-1,],]
    sumofsq = crossprod(Y-mu_i_long)
    Vinv = diag(beta, d)
    lambda_i = rWishart(1, df = n + alpha, Sigma = solve(Vinv + sumofsq))[,,1]
    df_sim_lambda[[i]] = lambda_i
    sigma_i = solve(lambda_i)
    #beta_i = beta + sum(sapply(1:q, function(x){sum((df$Y[index_1[x,]]-mu_i[x])^2)}))/2
    #lambda_i = rgamma(1, alpha_n, beta_i)
    #df_sim[i,q+2] = lambda_i
    
    #z
    df_sim_z[i,] = df_sim_z[i-1, ]
    for (j in 1:n){
      z_j_prev = df_sim_z[i,j]
      qlessk = setdiff(1:q, z_j_prev)
      z_j_new = sample(qlessk, 1)
      j_vector = df_j[[j]]
      h_z_prev = gamma/length(j_vector)* 2*sum(((z_j_prev == df_sim_z[i, j_vector])-0.5)) + dmvnorm(Y[j,], mean = mu_i[z_j_prev,], sigma = sigma_i, log = T)
      h_z_new = gamma/length(j_vector) * 2*sum(((z_j_new  == df_sim_z[i, j_vector])-0.5)) + dmvnorm(Y[j,], mean = mu_i[z_j_new, ] , sigma = sigma_i, log = T)
      prob_j = min(exp(h_z_new - h_z_prev),1)
      df_sim_z[i, j] = sample(x = c(z_j_prev, z_j_new), size = 1, prob = c(1-prob_j, prob_j))
    }
  }
  return(list(z = df_sim_z, mu = df_sim_mu, lambda = df_sim_lambda))
}


run_mcmc_deconv = function(df, nrep = 1000, q = 3, d = 2, mu0 = colMeans(df[,grep("Y",colnames(df))]), lambda0 = diag(0.01, nrow = d), alpha = 1, beta = 0.01, gamma = 2, seed = 100, prev){
  set.seed(seed)
  n = nrow(df)
  df[,"j"] = 1:n
  df_j = sapply(1:n, function(x){df[(abs(df[,"x"] -df[x,"x"]) + abs(df[,"y"] - df[x,"y"])) == 2,"j"]})
  
  df2 = df[rep(seq_len(n), 4), ] #rbind 4 times
  df2[,"j2"] = NA
  df2[,"j2"] = rep(1:(4*n))
  
  shift = expand.grid(c(1/3, -1/3), c(1/3,-1/3))
  shift_long = shift[rep(seq_len(4), each = n), ]
  
  df2[,"x"] = df2[,"x"] + shift_long$Var1
  df2[,"y"] = df2[,"y"] + shift_long$Var2
  df2_j = sapply(1:(4*n), function(x){df2[(abs(df2[,"x"] -df2[x,"x"]) + abs(df2[,"y"] - df2[x,"y"])) <= 2/3 + 0.01 &
                                            (abs(df2[,"x"] -df2[x,"x"]) + abs(df2[,"y"] - df2[x,"y"])) > 0,"j2"]})
  n = nrow(df2)
  n0 = nrow(df)
  
  #Initialize matrices storing iterations
  df_sim_z = matrix(NA, nrow = nrep, ncol = n)
  colnames(df_sim_z) = c(paste0("z",1:n))
  df_sim_mu = matrix(NA, nrow = nrep, ncol = q*d)
  qd = expand.grid(d = 1:d, q = 1:q)
  colnames(df_sim_mu) = sapply(1:(q*d), function(x){paste0("muq", qd[x,2],"_d",qd[x,1])})
  df_sim_lambda = replicate(nrep, matrix(NA, nrow = d, ncol = d), simplify = F)
  df_sim_Y = replicate(nrep, matrix(NA, nrow = n, ncol = d), simplify = F)
  
  lambdamu0 = lambda0*mu0
  alpha_n = alpha + n/2
  
  #Initialize parameters
  df_sim_mu[1, ] = prev$mu[nrow(prev$mu),]/4
  df_sim_lambda[[1]] = prev$lambda[[length(prev$lambda)]]*16
  df_sim_z[1,] = rep(prev$z[nrow(prev$z),], 4)
  df_sim_Y[[1]] = as.matrix(df2[,grep("Y", colnames(df2))]/4)
  
  Y = as.matrix(df[,grep("Y", colnames(df))])
  
  testY = rep(NA, nrep)
  testY[1] = 0
  #Iterate
  for (i in 2:nrep){
    #mu
    mu_i = matrix(nrow = q, ncol = d)
    index_1 = matrix(nrow = q, ncol = n)
    lambda_prev = df_sim_lambda[[i-1]]
    for(k in 1:q){
      index_1[k,] = df_sim_z[i-1,]==k
      n_i = sum(index_1[k,])
      mean_i = solve(lambda0 + n_i * lambda_prev) %*% (lambda0 %*% mu0 + lambda_prev %*% colSums(df_sim_Y[[i-1]][index_1[k,],, drop = F]))
      var_i = solve(lambda0 + n_i * lambda_prev)
      mu_i[k,] = rmvnorm(1, mean = mean_i, sigma = var_i)
    }
    df_sim_mu[i,] = as.vector(t(mu_i))
    
    #lambda
    mu_i_long = mu_i[df_sim_z[i-1,],]
    sumofsq = crossprod(df_sim_Y[[i-1]]-mu_i_long)
    Vinv = diag(beta, d)
    lambda_i = rWishart(1, df = n + alpha, Sigma = solve(Vinv + sumofsq))[,,1]
    df_sim_lambda[[i]] = lambda_i
    sigma_i = solve(lambda_i)
    
    #Y
    four_map = lapply(1:n0, function(x){which(df2[,"j"] == x)})
    test = rep(NA, n0) #testing purpose only!!
    for (j in 1:n0){
      Y_j_prev = df_sim_Y[[i-1]][four_map[[j]],]
      error = scale(rmvnorm(n = 4, rep(0, d), sigma = diag(d)/5), scale = F)
      Y_j_new = Y_j_prev + error
      mu_i_four = mu_i[df_sim_z[i-1,j + 0:3 * n0],]
      p_prev = prod(sapply(1:4, function(x){dmvnorm(Y_j_prev[x,], mu_i_four[x,], sigma_i)}))
      p_new = prod(sapply(1:4, function(x){dmvnorm(Y_j_new[x,], mu_i_four[x,], sigma_i)}))
      #probY_j = min(p_new/p_prev, 1)
      probY_j = min(p_new/p_prev * exp(-0.1*(sum(diag(crossprod(df_sim_Y[[1]][c(j,j,j,j),] - Y_j_new))) -
                                           sum(diag(crossprod(df_sim_Y[[1]][c(j,j,j,j),] - Y_j_prev))))), 1)
      test[j] = sample(x = 0:1, size = 1, prob = c(1-probY_j, probY_j)) #(testing purpose only!!)
      if (sample(x = 0:1, size = 1, prob = c(1-probY_j, probY_j))){
        df_sim_Y[[i]][four_map[[j]],] = Y_j_new
      } else {
        df_sim_Y[[i]][four_map[[j]],] = Y_j_prev
      }
    }
    testY[i] = prop.table(table(test))["1"]
    #z [CHANGE j to j2!!]
    df_sim_z[i,] = df_sim_z[i-1, ]
    for (j2 in 1:n){
      z_j_prev = df_sim_z[i,j2]
      qlessk = setdiff(1:q, z_j_prev)
      z_j_new = sample(qlessk, 1)
      j_vector = df2_j[[j2]]
      h_z_prev = gamma/length(j_vector)* 2*sum(((z_j_prev == df_sim_z[i, j_vector])-0.5)) + dmvnorm(df_sim_Y[[i]][j2,], mean = mu_i[z_j_prev,], sigma = sigma_i, log = T)
      h_z_new = gamma/length(j_vector) * 2*sum(((z_j_new  == df_sim_z[i, j_vector])-0.5)) + dmvnorm(df_sim_Y[[i]][j2,], mean = mu_i[z_j_new, ], sigma = sigma_i, log = T)
      prob_j = min(exp(h_z_new - h_z_prev),1)
      df_sim_z[i, j2] = sample(x = c(z_j_prev, z_j_new), size = 1, prob = c(1-prob_j, prob_j))
    }
  }
  return(list(z = df_sim_z, mu = df_sim_mu, lambda = df_sim_lambda, Y = df_sim_Y))
}

run_mcmc_hexdeconv = function(df, nrep = 1000, q = 3, d = 2, mu0 = colMeans(df[,grep("Y",colnames(df))]), lambda0 = diag(0.01, nrow = d), alpha = 1, beta = 0.01, gamma = 2, seed = 100, prev){
  set.seed(seed)
  n = nrow(df)
  df[,"j"] = 1:n
  df_j = sapply(1:n, function(x){df[(abs(df[,"x"] -df[x,"x"]) + abs(df[,"y"] - df[x,"y"])) == 2,"j"]})
  
  df2 = as.data.frame(df[rep(seq_len(n), 7), ]) #rbind 7 times
  df2$j2 = rep(1:(7*n))
  
  shift = rbind(expand.grid(c(1/2, -1/2), c(1/4,-1/4)), expand.grid(0, c(0, 1/2, -1/2)))
  shift_long = shift[rep(seq_len(7), each = n), ]
  
  df2$x = df2$x + shift_long$Var1
  df2$y = df2$y + shift_long$Var2
  df2 = as.matrix(df2)
  df2_j = sapply(1:(7*n), function(x){df2[(abs(df2[,"x"] -df2[x,"x"]) + abs(df2[,"y"] - df2[x,"y"])) <= 0.75 &
                                            (abs(df2[,"x"] -df2[x,"x"]) + abs(df2[,"y"] - df2[x,"y"])) > 0,"j2"]})
  n = nrow(df2)
  n0 = nrow(df)
  
  #Initialize matrices storing iterations
  df_sim_z = matrix(NA, nrow = nrep, ncol = n)
  colnames(df_sim_z) = c(paste0("z",1:n))
  df_sim_mu = matrix(NA, nrow = nrep, ncol = q*d)
  qd = expand.grid(d = 1:d, q = 1:q)
  colnames(df_sim_mu) = sapply(1:(q*d), function(x){paste0("muq", qd[x,2],"_d",qd[x,1])})
  df_sim_lambda = replicate(nrep, matrix(NA, nrow = d, ncol = d), simplify = F)
  df_sim_Y = replicate(nrep, matrix(NA, nrow = n, ncol = d), simplify = F)
  
  lambdamu0 = lambda0*mu0
  alpha_n = alpha + n/2
  
  #Initialize parameters
  df_sim_mu[1, ] = prev$mu[nrow(prev$mu),]/7
  df_sim_lambda[[1]] = prev$lambda[[length(prev$lambda)]]*(7^2)
  df_sim_z[1,] = rep(prev$z[nrow(prev$z),], 7)
  df_sim_Y[[1]] = as.matrix(df2[,grep("Y", colnames(df2))]/7)
  
  Y = as.matrix(df[,grep("Y", colnames(df))])
  
  testY = rep(NA, nrep)
  testY[1] = 0
  #Iterate
  for (i in 2:nrep){
    #mu
    mu_i = matrix(nrow = q, ncol = d)
    index_1 = matrix(nrow = q, ncol = n)
    lambda_prev = df_sim_lambda[[i-1]]
    for(k in 1:q){
      index_1[k,] = df_sim_z[i-1,]==k
      n_i = sum(index_1[k,])
      mean_i = solve(lambda0 + n_i * lambda_prev) %*% (lambda0 %*% mu0 + lambda_prev %*% colSums(df_sim_Y[[i-1]][index_1[k,],, drop = F]))
      var_i = solve(lambda0 + n_i * lambda_prev)
      mu_i[k,] = rmvnorm(1, mean = mean_i, sigma = var_i)
    }
    df_sim_mu[i,] = as.vector(t(mu_i))
    
    #lambda
    mu_i_long = mu_i[df_sim_z[i-1,],]
    sumofsq = crossprod(df_sim_Y[[i-1]]-mu_i_long)
    Vinv = diag(beta, d)
    lambda_i = rWishart(1, df = n + alpha, Sigma = solve(Vinv + sumofsq))[,,1]
    df_sim_lambda[[i]] = lambda_i
    sigma_i = solve(lambda_i)
    
    #Y
    four_map = lapply(1:n0, function(x){which(df2[,"j"] == x)})
    test = rep(NA, n0) #testing purpose only!!
    for (j in 1:n0){
      Y_j_prev = df_sim_Y[[i-1]][four_map[[j]],]
      error = scale(rmvnorm(n = 7, rep(0, d), sigma = diag(d)/5), scale = F)
      Y_j_new = Y_j_prev + error
      mu_i_four = mu_i[df_sim_z[i-1,j + 0:6 * n0],]
      p_prev = prod(sapply(1:7, function(x){dmvnorm(Y_j_prev[x,], mu_i_four[x,], sigma_i)}))
      p_new = prod(sapply(1:7, function(x){dmvnorm(Y_j_new[x,], mu_i_four[x,], sigma_i)}))
      #probY_j = min(p_new/p_prev, 1)
      probY_j = min(p_new/p_prev * exp(-0.1*(sum(diag(crossprod(df_sim_Y[[1]][rep(j,7),] - Y_j_new))) -
                                               sum(diag(crossprod(df_sim_Y[[1]][rep(j,7),] - Y_j_prev))))), 1)
      test[j] = sample(x = 0:1, size = 1, prob = c(1-probY_j, probY_j)) #(testing purpose only!!)
      if (sample(x = 0:1, size = 1, prob = c(1-probY_j, probY_j))){
        df_sim_Y[[i]][four_map[[j]],] = Y_j_new
      } else {
        df_sim_Y[[i]][four_map[[j]],] = Y_j_prev
      }
    }
    testY[i] = prop.table(table(test))["1"]
    #z [CHANGE j to j2!!]
    df_sim_z[i,] = df_sim_z[i-1, ]
    for (j2 in 1:n){
      z_j_prev = df_sim_z[i,j2]
      qlessk = setdiff(1:q, z_j_prev)
      z_j_new = sample(qlessk, 1)
      j_vector = df2_j[[j2]]
      h_z_prev = gamma/length(j_vector)* 2*sum(((z_j_prev == df_sim_z[i, j_vector])-0.5)) + dmvnorm(df_sim_Y[[i]][j2,], mean = mu_i[z_j_prev,], sigma = sigma_i, log = T)
      h_z_new = gamma/length(j_vector) * 2*sum(((z_j_new  == df_sim_z[i, j_vector])-0.5)) + dmvnorm(df_sim_Y[[i]][j2,], mean = mu_i[z_j_new, ], sigma = sigma_i, log = T)
      prob_j = min(exp(h_z_new - h_z_prev),1)
      df_sim_z[i, j2] = sample(x = c(z_j_prev, z_j_new), size = 1, prob = c(1-prob_j, prob_j))
    }
  }
  return(list(z = df_sim_z, mu = df_sim_mu, lambda = df_sim_lambda, Y = df_sim_Y))
}
