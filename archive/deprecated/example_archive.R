library(tidyverse)
library(patchwork)
set.seed(100)

#mat = matrix(0, nrow = 100, ncol = 100)
df = expand.grid(x = 1:100, y = 1:100)
df$z = -1

rect_center = sample(10:90, 2, replace = T)
df$z[df$x %in% seq(rect_center[1] - 5, rect_center[1]+5) & df$y %in% seq(rect_center[2] - 5, rect_center[2]+5)] = 1

circ_center = sample(10:90, 2, replace = T)
r = 10 #radius
s = seq(0,r) #sagitta
l = sqrt(2*(r-s)*r-(r-s)^2) #half chord length
for (i in s){
  df$z[(df$x %in% c(circ_center[1] + c(i, -i))) & 
         (df$y > circ_center[2] - l[i+1]) & 
         (df$y < circ_center[2] + l[i+1])] = 1
}

p1 = ggplot(df, aes(x,y,fill = factor(z, levels = c(1,-1)))) + 
  geom_tile() +
  scale_fill_manual(values = c("#b2182b", "#2166ac")) +
  theme_classic() +
  labs(fill = "State", x = NULL, y = NULL)

set.seed(100)  
mu = 5
lambda = 0.5
df$Y = rnorm(n = nrow(df), 
             mean = (df$z==1)*mu, 
             sd = 1/sqrt(lambda))

df$j = 1:nrow(df)+3
p2 = ggplot(df, aes(x,y,fill = ifelse(Y>5, 5, ifelse(Y<0,0, Y)))) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", labels = c(expression(""<= 0), 1:4, expression("">=5)))+
  theme_classic()+
  labs(fill = "Expression", x = NULL, y = NULL)


#mcmc
ptm = proc.time()
set.seed(100)
nrep = 10
df_sim = matrix(NA, nrow = nrep, ncol = nrow(df) + 3)
colnames(df_sim) = c("i", "mu", "lambda", paste0("z",1:nrow(df)))
mu0 = mean(df$Y)
n = nrow(df)
df_j = sapply(1:n, function(x){df[abs(df[,"x"] -df[x,"x"]) + abs(df[,"y"] - df[x,"y"]) == 1,"j"]})
lambda0 = 1/100
alpha = 1
beta = 0.01
gamma = 2
df_sim[,"i"] = 1:nrep
df_sim[1,"mu"] = mu0
df_sim[1,"lambda"] = 1/var(df$Y)
df_sim[1,4:ncol(df_sim)] = -1
lambdamu0 = lambda0*mu0
alpha_n = alpha + n/2
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
proc.time()-ptm

p1+p2 + plot_annotation(tag_levels = "A")

df$test = df_sim[i,-(1:3)]
df$test2 = test2[10,-(1:3)]

df$z_gamma2 = colMeans(df_sim_gamma2[-(1:100),-(1:3)])
df$z_gamma4 = colMeans(df_sim_gamma4[-(1:100),-(1:3)])
df$z_gamma6 = colMeans(df_sim_gamma6[-(1:100),-(1:3)])

ggplot(df, aes(x,y,fill = factor(test))) + 
  geom_tile() +
  scale_fill_discrete(breaks = c(1, -1)) +
  theme_classic()
ggplot(df, aes(x,y,fill = factor(test2))) + 
  geom_tile() +
  scale_fill_discrete(breaks = c(1, -1)) +
  theme_classic()

z_2 = ggplot(df, aes(x,y,fill = z_gamma2)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", breaks = c(-1,1))+
  theme_classic() +
  labs(x = NULL, y = NULL, fill = "State")
z_4 = ggplot(df, aes(x,y,fill = z_gamma4)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", breaks = c(-1,1))+
  theme_classic() +
  labs(x = NULL, y = NULL, fill = "State")
z_6 = ggplot(df, aes(x,y,fill = z_gamma6)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", breaks = c(-1,1))+
  theme_classic() +
  labs(x = NULL, y = NULL, fill = "State")

mu_2 = ggplot(as_tibble(df_sim_gamma2[-(1:100),]), aes(x = i, y = mu)) + geom_line() +
  theme_classic() +
  labs(x = NULL, y = expression(mu)) + 
  geom_hline(yintercept = 5, color = "red") +
  coord_cartesian(ylim = c(4.6, 5.3))
mu_4 = ggplot(as_tibble(df_sim_gamma4[-(1:100),]), aes(x = i, y = mu)) + geom_line() +
  theme_classic() +
  labs(x = NULL, y = expression(mu)) + 
  geom_hline(yintercept = 5, color = "red") +
  coord_cartesian(ylim = c(4.6, 5.3))
mu_6 = ggplot(as_tibble(df_sim_gamma6[-(1:100),]), aes(x = i, y = mu)) + geom_line() +
  theme_classic() +
  labs(x = NULL, y = expression(mu)) + 
  geom_hline(yintercept = 5, color = "red") +
  coord_cartesian(ylim = c(4.6, 5.3))

lambda_2 = ggplot(as_tibble(df_sim_gamma2[-(1:100),]), aes(x = i, y = lambda)) + geom_line() +
  theme_classic() +
  labs(x = NULL, y = expression(lambda)) + 
  geom_hline(yintercept = 0.5, color = "red") +
  coord_cartesian(ylim = c(0.48, 0.55))

lambda_4 = ggplot(as_tibble(df_sim_gamma4[-(1:100),]), aes(x = i, y = lambda)) + geom_line() +
  theme_classic() +
  labs(x = NULL, y = expression(lambda)) + 
  geom_hline(yintercept = 0.5, color = "red") +
  coord_cartesian(ylim = c(0.48, 0.55))

lambda_6 = ggplot(as_tibble(df_sim_gamma6[-(1:100),]), aes(x = i, y = lambda)) + geom_line() +
  theme_classic() +
  labs(x = NULL, y = expression(lambda)) + 
  geom_hline(yintercept = 0.5, color = "red") +
  coord_cartesian(ylim = c(0.48, 0.55))

((z_2+z_4+z_6) +plot_layout(guides = 'collect')) / (mu_2 + mu_4 + mu_6) / (lambda_2 + lambda_4 + lambda_6) + 
  plot_annotation(tag_levels = "A")
df_sim_gamma2[,"gamma"] = 2
df_sim_gamma4[,"gamma"] = 4
df_sim_gamma6[,"gamma"] = 6
sim_mat_tall = as_tibble(rbind(df_sim_gamma2, df_sim_gamma4, df_sim_gamma6))
sim_mat_tall %>%
  gather
