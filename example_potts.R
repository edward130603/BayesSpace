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
  scale_fill_manual(values = c("#f6f6f6", "#b2182b", "#2166ac")) +
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
df_sim_gamma1 = run_mcmc_potts(df = df, gamma = 2, q = 3, nrep = 1000)

df$test2 = df_sim_gamma1[nrow(df_sim_gamma1), 6:ncol(df_sim_gamma1)]
ggplot(df, aes(x,y,fill = factor(test2))) + 
  geom_tile()
