library(tidyverse)
library(mvtnorm)
library(patchwork)
set.seed(102)

#Generate truth
df = expand.grid(x = 1:100, y = 1:100)
df = df[df$x %% 2 + df$y %% 2 == 1,]
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
p1 = ggplot(df, aes(x,y,color = factor(z))) + 
  geom_point() +
  #scale_fill_manual(values = c("#f6f6f6", "#b2182b", "#2166ac")) +
  theme_classic() +
  theme(legend.position = "bottom")+
  labs(color = "State", x = NULL, y = NULL)

#Generate simulated expression
set.seed(101)  
mu1 = c(0,0)
mu2 = c(5,0)
mu3 = c(5,5)
Sigma = diag(1, nrow = 2) + 1
Lambda = solve(Sigma)
rmvnorm(5, mean = mu1, sigma = Sigma)
df$Y1 = NA
df$Y2 = NA
df[df$z == 1,c("Y1", "Y2")] = rmvnorm(n = sum(df$z == 1), 
                                      mean = mu1, 
                                      sigma = Sigma)
df[df$z == 2,c("Y1", "Y2")] = rmvnorm(n = sum(df$z == 2), 
                                      mean = mu2, 
                                      sigma = Sigma)
df[df$z == 3,c("Y1", "Y2")] = rmvnorm(n = sum(df$z == 3), 
                                      mean = mu3, 
                                      sigma = Sigma)


#Plots
p2 = ggplot(df, aes(x,y, color = Y1)) + 
  geom_point() +
  scale_color_distiller(palette = "Blues", direction = 1, breaks = c(0,5)) +
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(color = expression(Y[1]), x = NULL, y = NULL)
p3 = ggplot(df, aes(x,y, color = Y2)) + 
  geom_point() +
  scale_color_distiller(palette = "Blues", direction = 1, breaks = c(0,5)) +
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(color = expression(Y[2]), x = NULL, y = NULL)

p1+p2+p3 + plot_annotation(tag_levels = "A")


