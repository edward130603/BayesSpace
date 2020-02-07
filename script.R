library(mvtnorm)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

run_mcmc = function(df, nrep = 1000, mu0 = mean(df[,"Y"]), lambda0 = 1/100, alpha = 1, beta = 0.01, gamma = 2, seed = 100){
  set.seed(seed)
  
  n = nrow(df)
  df$j = 1:nrow(df)+3
  df_j = sapply(1:n, function(x){df[abs(df[,"x"] -df[x,"x"]) + abs(df[,"y"] - df[x,"y"]) == 1,"j"]})
  
  #Initialize matrix storing iterations
  df_sim = matrix(NA, nrow = nrep, ncol = nrow(df)+3)
  colnames(df_sim) = c("i", "mu", "lambda", paste0("z",1:n))
  df_sim[,"i"] = 1:nrep
  lambdamu0 = lambda0*mu0
  alpha_n = alpha + n/2
  
  #Initialize parameters
  df_sim[1,"mu"] = mu0
  df_sim[1,"lambda"] = 1/var(df[,"Y"])
  df_sim[1,4:ncol(df_sim)] = -1
  
  #Iterate
  for (i in 2:nrep){
    #mu
    lambda_prev = df_sim[i-1,3]
    index_1 = df_sim[i-1,4:ncol(df_sim)]==1
    n_i = sum(index_1)
    mean_i = (lambdamu0 + lambda_prev*sum(df$Y[index_1]))/(lambda0 + n_i*lambda_prev)
    var_i = 1/(lambda0 + n_i * lambda_prev)
    mu_i = rnorm(1, mean = mean_i, sd = sqrt(var_i))
    df_sim[i, 2] = mu_i
    
    #lambda
    beta_i = beta+ sum((df$Y-mu_i*index_1)^2)/2
    lambda_i = rgamma(1, alpha_n, beta_i)
    df_sim[i,3] = lambda_i
    
    #z
    df_sim[i, 4:ncol(df_sim)] = df_sim[i-1, 4:ncol(df_sim)]
    for (j in 1:n){
      z_j_prev = df_sim[i,j+3]
      z_j_new = z_j_prev*-1
      x_j = df$x[j]
      y_j = df$y[j]
      j_vector = df_j[[j]]
      h_z_prev = gamma/length(j_vector)* sum(z_j_prev * df_sim[i, j_vector]) + dnorm(df$Y[j], mean = mu_i * (z_j_prev == 1), sd = 1/sqrt(lambda_i), log = T)
      h_z_new = gamma/length(j_vector) * sum(z_j_new * df_sim[i, j_vector]) + dnorm(df$Y[j], mean = mu_i * (z_j_new == 1), sd = 1/sqrt(lambda_i), log = T)
      prob_j = min(exp(h_z_new - h_z_prev),1)
      df_sim[i, j+3] = sample(x = c(z_j_new, z_j_prev), size = 1, prob = c(prob_j, 1-prob_j))
    }
  }
  return(df_sim)
}

run_mcmc_potts = function(df, nrep = 1000, q = 3, mu0 = mean(df[,"Y"]), lambda0 = 1/100, alpha = 1, beta = 0.01, gamma = 2, seed = 100){
  set.seed(seed)
  
  n = nrow(df)
  df$j = 1:nrow(df)+q+2
  df_j = sapply(1:n, function(x){df[(abs(df[,"x"] -df[x,"x"]) + abs(df[,"y"] - df[x,"y"])) == 1,"j"]})
  
  #Initialize matrix storing iterations
  df_sim = matrix(NA, nrow = nrep, ncol = nrow(df)+2+q)
  colnames(df_sim) = c("i", paste0("mu", 1:q), "lambda", paste0("z",1:n))
  df_sim[,"i"] = 1:nrep
  lambdamu0 = lambda0*mu0
  alpha_n = alpha + n/2
  
  #Initialize parameters
  df_sim[1,2:(q+1)] = mu0
  df_sim[1,"lambda"] = 1/var(df[,"Y"])
  df_sim[1,(q+3):ncol(df_sim)] = 1 #initialize to all 1
  #df_sim[1,(q+3):ncol(df_sim)] = df$z #initialize with truth (testing only)
  #df_sim[1,(q+3):ncol(df_sim)] = sample(1:3, n, replace = T) #random init
  
  #Iterate
  for (i in 2:nrep){
    #mu
    mu_i = numeric(q)
    index_1 = matrix(nrow = q, ncol = n)
    for(k in 1:q){
      lambda_prev = df_sim[i-1,q+2]
      index_1[k,] = df_sim[i-1,(q+3):ncol(df_sim)]==k
      n_i = sum(index_1[k,])
      mean_i = (lambdamu0 + lambda_prev*sum(df$Y[index_1[k,]]))/(lambda0 + n_i*lambda_prev)
      var_i = 1/(lambda0 + n_i * lambda_prev)
      mu_i[k] = rnorm(1, mean = mean_i, sd = sqrt(var_i))
      df_sim[i, 1+k] = mu_i[k]
    }
    
    #lambda
    beta_i = beta + sum(sapply(1:q, function(x){sum((df$Y[index_1[x,]]-mu_i[x])^2)}))/2
    lambda_i = rgamma(1, alpha_n, beta_i)
    df_sim[i,q+2] = lambda_i
    
    #z
    df_sim_z[i,] = df_sim_z[i-1,]
    for (j in 1:n){
      z_j_prev = df_sim_z[i,j]
      qlessk = setdiff(1:q, z_j_prev)
      z_j_new = sample(qlessk, 1)
      j_vector = df_j[[j]]
      h_z_prev = gamma/length(j_vector)* 2*sum(((z_j_prev == df_sim[i, j_vector])-0.5)) + dnorm(df$Y[j], mean = mu_i[z_j_prev], sd = 1/sqrt(lambda_i), log = T)
      h_z_new = gamma/length(j_vector) * 2*sum(((z_j_new  == df_sim[i, j_vector])-0.5)) + dnorm(df$Y[j], mean = mu_i[z_j_new] , sd = 1/sqrt(lambda_i), log = T)
      prob_j = min(exp(h_z_new - h_z_prev),1)
      df_sim[i, j+2+q] = sample(x = c(z_j_prev, z_j_new), size = 1, prob = c(1-prob_j, prob_j))
    }
  }
  return(df_sim)
}

run_mcmc_multi = function(df, nrep = 1000, q = 3, d = 2, mu0 = colMeans(df[,grep("Y",colnames(df))]), lambda0 = diag(0.01, nrow = d), alpha = 1, beta = 0.01, gamma = 2, seed = 100){
  set.seed(seed)
  n = nrow(df)
  df$j = 1:n
  df_j = sapply(1:n, function(x){df[(abs(df[,"x"] -df[x,"x"]) + abs(df[,"y"] - df[x,"y"])) == 2,"j"]})
  
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
  #df_sim_z[1,] = sample(3, n, replace = T) #initialize to all truth (testing only)
  
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
