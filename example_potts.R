library(tidyverse)
library(patchwork)
set.seed(101)

#Generate truth
df = expand.grid(x = 1:100, y = 1:100)
df$z = 1

rect_center = sample(10:90, 2, replace = T)
df$z[df$x %in% seq(rect_center[1] - 5, rect_center[1]+5) & df$y %in% seq(rect_center[2] - 5, rect_center[2]+5)] = 2

circ_center = sample(10:90, 2, replace = T)
r = 10 #radius
s = seq(0,r) #sagitta
l = sqrt(2*(r-s)*r-(r-s)^2) #half chord length
for (i in s){
  df$z[(df$x %in% c(circ_center[1] + c(i, -i))) & 
         (df$y > circ_center[2] - l[i+1]) & 
         (df$y < circ_center[2] + l[i+1])] = 3
}

#Plot truth
p1 = ggplot(df, aes(x,y,fill = factor(z))) + 
  geom_tile() +
  #scale_fill_manual(values = c("#f6f6f6", "#b2182b", "#2166ac")) +
  theme_classic() +
  labs(fill = "State", x = NULL, y = NULL)

#Generate simulated expression
set.seed(101)  
mu2 = 5
mu3 = -5
lambda = 0.5
df$Y = rnorm(n = nrow(df), 
             mean = (df$z==2)*mu2 + (df$z==3)*mu3, 
             sd = 1/sqrt(lambda))

#Plots
p2 = ggplot(df, aes(x,y,fill = ifelse(Y>5, 5, ifelse(Y< -5,-5, Y)))) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", labels = c(expression(""<= -5), seq(-2.5, 2.5, 2.5), expression("">=5)))+
  theme_classic()+
  labs(fill = "Expression", x = NULL, y = NULL)

p1+p2 + plot_annotation(tag_levels = "A")

#Run mcmc

df_sim_gamma1 = run_mcmc_potts(df = df, gamma = 1, q = 3)
df_sim_gamma2 = run_mcmc_potts(df = df, gamma = 2, q = 3)
df_sim_gamma3 = run_mcmc_potts(df = df, gamma = 3, q = 3)
df_sim_gamma4 = run_mcmc_potts(df = df, gamma = 4, q = 3)
df_sim_gamma6 = run_mcmc_potts(df = df, gamma = 6, q = 3)

 #Plot mcmc means
df$z_gamma1 = apply(df_sim_gamma1[-(1:100),-(1:5)], 2, Mode)
df$z_gamma1_alpha = pmax(colMeans(df_sim_gamma1[-(1:100),-(1:5)]==1),
                         colMeans(df_sim_gamma1[-(1:100),-(1:5)]==2),
                         colMeans(df_sim_gamma1[-(1:100),-(1:5)]==3))
df$z_gamma3 = apply(df_sim_gamma3[-(1:100),-(1:5)], 2, Mode)
df$z_gamma3_alpha = pmax(colMeans(df_sim_gamma3[-(1:100),-(1:5)]==1),
                         colMeans(df_sim_gamma3[-(1:100),-(1:5)]==2),
                         colMeans(df_sim_gamma3[-(1:100),-(1:5)]==3))

df$z_gamma2 = apply(df_sim_gamma2[-(1:100),-(1:5)], 2, Mode)
df$z_gamma2_alpha = pmax(colMeans(df_sim_gamma2[-(1:100),-(1:5)]==1),
                         colMeans(df_sim_gamma2[-(1:100),-(1:5)]==2),
                         colMeans(df_sim_gamma2[-(1:100),-(1:5)]==3))
df$z_gamma3 = apply(df_sim_gamma3[-(1:100),-(1:5)], 2, Mode)
df$z_gamma3_alpha = pmax(colMeans(df_sim_gamma3[-(1:100),-(1:5)]==1),
                         colMeans(df_sim_gamma3[-(1:100),-(1:5)]==2),
                         colMeans(df_sim_gamma3[-(1:100),-(1:5)]==3))
df$z_gamma4 = apply(df_sim_gamma4[-(1:100),-(1:5)], 2, Mode)
df$z_gamma4_alpha = pmax(colMeans(df_sim_gamma4[-(1:100),-(1:5)]==1),
                         colMeans(df_sim_gamma4[-(1:100),-(1:5)]==2),
                         colMeans(df_sim_gamma4[-(1:100),-(1:5)]==3))
df$z_gamma6 = apply(df_sim_gamma6[-(1:100),-(1:5)], 2, Mode)
df$z_gamma6_alpha = pmax(colMeans(df_sim_gamma6[-(1:100),-(1:5)]==1),
                         colMeans(df_sim_gamma6[-(1:100),-(1:5)]==2),
                         colMeans(df_sim_gamma6[-(1:100),-(1:5)]==3))

potts2 = ggplot(df, aes(x, y)) +
  geom_tile(aes(fill = factor(z_gamma2), alpha = z_gamma2_alpha)) +
  labs(fill = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+
  theme_classic()
potts3 = ggplot(df, aes(x, y)) +
  geom_tile(aes(fill = factor(z_gamma3), alpha = z_gamma3_alpha)) +
  labs(fill = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+
  theme_classic()
potts4 = ggplot(df, aes(x, y)) +
  geom_tile(aes(fill = factor(z_gamma4), alpha = z_gamma4_alpha)) +
  labs(fill = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+
  theme_classic()

#Calculate MLEs
reg1 = lm(data = df, Y~factor(z))
mu1_mle = sum(reg1$coefficients[c(1)])
mu2_mle = sum(reg1$coefficients[c(1,2)])
mu3_mle = sum(reg1$coefficients[c(1,3)])
lambda_mle = 1/summary(reg1)$sigma^2

#Trace plots
mu2_2 = ggplot(as_tibble(df_sim_gamma2[-(1:100),]), aes(x = i, y = mu2)) + geom_line() +
  theme_classic() +
  labs(x = NULL, y = expression(mu[2])) + 
  geom_hline(yintercept = 5, color = "red") +
  geom_hline(yintercept = mu2_mle, color = "blue") +
  coord_cartesian(ylim = c(-5, 5))
mu2_3 = ggplot(as_tibble(df_sim_gamma3[-(1:100),]), aes(x = i, y = mu2)) + geom_line() +
  theme_classic() +
  labs(x = NULL, y = expression(mu[2])) + 
  geom_hline(yintercept = 5, color = "red") +
  geom_hline(yintercept = mu2_mle, color = "blue") +
  coord_cartesian(ylim = c(-5, 5))
mu2_4 = ggplot(as_tibble(df_sim_gamma4[-(1:100),]), aes(x = i, y = mu3)) + geom_line() +
  theme_classic() +
  labs(x = NULL, y = expression(mu[2])) + 
  geom_hline(yintercept = 5, color = "red") +
  geom_hline(yintercept = mu2_mle, color = "blue") +
  coord_cartesian(ylim = c(-5, 5))

((potts2+potts3+potts4)) /
  (mu2_2 + mu2_3 + mu2_4) +
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "A")
