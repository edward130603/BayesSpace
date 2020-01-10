run_mcmc = function(df, nrep = 1000, mu0 = mean(df[,"Y"]), lambda0 = 1/100, alpha = 1, beta = 0.01, gamma = 2, seed = 100){
  set.seed(seed)
  
  n = nrow(df)
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
  df_j = sapply(1:n, function(x){df[abs(df[,"x"] -df[x,"x"]) + abs(df[,"y"] - df[x,"y"]) == 1,"j"]})
  
  #Initialize matrix storing iterations
  df_sim = matrix(NA, nrow = nrep, ncol = nrow(df)+2+q)
  colnames(df_sim) = c("i", paste0("mu", 1:q), "lambda", paste0("z",1:n))
  df_sim[,"i"] = 1:nrep
  lambdamu0 = lambda0*mu0
  alpha_n = alpha + n/2
  
  #Initialize parameters
  df_sim[1,2:(q+1)] = mu0
  df_sim[1,"lambda"] = 1/var(df[,"Y"])
  df_sim[1,(q+3):ncol(df_sim)] = 1
  
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
    sapply(1:q, function(x){sum((df$Y[index_1[x,]]-mu_i[x])^2)/2})
    beta_i = beta + sum(sapply(1:q, function(x){sum((df$Y[index_1[x,]]-mu_i[x])^2)/2}))
    lambda_i = rgamma(1, alpha_n, beta_i)
    df_sim[i,q+2] = lambda_i
    
    #z
    df_sim[i, (q+3):ncol(df_sim)] = df_sim[i-1, (q+3):ncol(df_sim)]
    for (j in 1:n){
      z_j_prev = df_sim[i,j+2+q]
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