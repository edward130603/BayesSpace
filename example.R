library(tidyverse)
library(patchwork)
set.seed(100)

#Generate truth
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

#Plot truth
p1 = ggplot(df, aes(x,y,fill = factor(z, levels = c(1,-1)))) + 
  geom_tile() +
  scale_fill_manual(values = c("#b2182b", "#2166ac")) +
  theme_classic() +
  labs(fill = "State", x = NULL, y = NULL)

#Generate simulated expression
set.seed(100)  
mu = 5
lambda = 0.5
df$Y = rnorm(n = nrow(df), 
             mean = (df$z==1)*mu, 
             sd = 1/sqrt(lambda))

df$j = 1:nrow(df)+3

#Plots
p2 = ggplot(df, aes(x,y,fill = ifelse(Y>5, 5, ifelse(Y<0,0, Y)))) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", labels = c(expression(""<= 0), 1:4, expression("">=5)))+
  theme_classic()+
  labs(fill = "Expression", x = NULL, y = NULL)

p1+p2 + plot_annotation(tag_levels = "A")

#Run mcmc
df_sim_gamma2 = run_mcmc(df = df, gamma = 2)
df_sim_gamma4 = run_mcmc(df = df, gamma = 4)
df_sim_gamma6 = run_mcmc(df = df, gamma = 6)

#Plot mcmc
df$z_gamma2 = colMeans(df_sim_gamma2[-(1:100),-(1:3)])
df$z_gamma4 = colMeans(df_sim_gamma4[-(1:100),-(1:3)])
df$z_gamma6 = colMeans(df_sim_gamma6[-(1:100),-(1:3)])

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

reg1 = lm(data = df, Y~I(z == 1) - 1)
mu_mle = reg1$coefficients[2]
lambda_mle = 1/summary(reg1)$sigma^2

mu_2 = ggplot(as_tibble(df_sim_gamma2[-(1:100),]), aes(x = i, y = mu)) + geom_line() +
  theme_classic() +
  labs(x = NULL, y = expression(mu)) + 
  geom_hline(yintercept = 5, color = "red") +
  geom_hline(yintercept = mu_mle, color = "blue") +
  coord_cartesian(ylim = c(4.6, 5.3))
mu_4 = ggplot(as_tibble(df_sim_gamma4[-(1:100),]), aes(x = i, y = mu)) + geom_line() +
  theme_classic() +
  labs(x = NULL, y = expression(mu)) + 
  geom_hline(yintercept = 5, color = "red") +
  geom_hline(yintercept = mu_mle, color = "blue") +
  coord_cartesian(ylim = c(4.6, 5.3))
mu_6 = ggplot(as_tibble(df_sim_gamma6[-(1:100),]), aes(x = i, y = mu)) + geom_line() +
  theme_classic() +
  labs(x = NULL, y = expression(mu)) + 
  geom_hline(yintercept = 5, color = "red") +
  geom_hline(yintercept = mu_mle, color = "blue") +
  coord_cartesian(ylim = c(4.6, 5.3))

lambda_2 = ggplot(as_tibble(df_sim_gamma2[-(1:100),]), aes(x = i, y = lambda)) + geom_line() +
  theme_classic() +
  labs(x = NULL, y = expression(lambda)) + 
  geom_hline(yintercept = 0.5, color = "red") +
  geom_hline(yintercept = lambda_mle, color = "blue") +
  coord_cartesian(ylim = c(0.48, 0.55))

lambda_4 = ggplot(as_tibble(df_sim_gamma4[-(1:100),]), aes(x = i, y = lambda)) + geom_line() +
  theme_classic() +
  labs(x = NULL, y = expression(lambda)) + 
  geom_hline(yintercept = 0.5, color = "red") +
  geom_hline(yintercept = lambda_mle, color = "blue") +
  coord_cartesian(ylim = c(0.48, 0.55))

lambda_6 = ggplot(as_tibble(df_sim_gamma6[-(1:100),]), aes(x = i, y = lambda)) + geom_line() +
  theme_classic() +
  labs(x = NULL, y = expression(lambda)) + 
  geom_hline(yintercept = 0.5, color = "red") +
  geom_hline(yintercept = lambda_mle, color = "blue") +
  coord_cartesian(ylim = c(0.48, 0.55))

((z_2+z_4+z_6) +plot_layout(guides = 'collect')) / (mu_2 + mu_4 + mu_6) / (lambda_2 + lambda_4 + lambda_6) + 
  plot_annotation(tag_levels = "A")
df_sim_gamma2 = cbind(df_sim_gamma2, gamma = 2)
df_sim_gamma4 = cbind(df_sim_gamma4, gamma = 4)
df_sim_gamma6 = cbind(df_sim_gamma6, gamma = 6)

sim_mat_tall = as_tibble(rbind(df_sim_gamma2, df_sim_gamma4, df_sim_gamma6))
sim_mat_tall %>%
  filter(i >=100) %>%
  ggplot(aes(x = factor(gamma), y = mu)) + 
  geom_hline(aes(yintercept = 5.0), color = "red") +
  geom_violin(draw_quantiles = 0.5, alpha = 0.5) +
  labs(x = expression(gamma), y = expression(mu))+
  coord_flip()

sim_mat_tall %>%
  filter(i >=100) %>%
  ggplot(aes(x = factor(gamma), y = lambda)) + 
  geom_hline(aes(yintercept = 0.5), color = "red") +
  geom_violin(draw_quantiles = 0.5, alpha = 0.5) +
  labs(x = expression(gamma), y = expression(lambda))+
  coord_flip()
