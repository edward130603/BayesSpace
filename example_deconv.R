library(tidyverse)
library(mvtnorm)
library(patchwork)
set.seed(102)

#Generate truth
df = expand.grid(x = 1:50, y = 1:50)
df = df[df$x %% 2 + df$y %% 2 == 1,]
df$z = 1
n = nrow(df)
df$j = 1:n

df2 = df[rep(seq_len(n), 4), ] #rbind 4 times
df2$j2 = rep(1:(4*n))

shift = expand.grid(c(1/3, -1/3), c(1/3,-1/3))
shift_long = shift[rep(seq_len(4), each = n), ]

df2$x = df2$x + shift_long$Var1
df2$y = df2$y + shift_long$Var2

rect_center = sapply(1:10,function(x)df2[sample(4*n, 1) ,c("x","y")])
for(i in 1:10){
  df2$z[df2$x > unlist(rect_center[,i])[1] - 1.1 & df2$x < unlist(rect_center[,i])[1]+1.1 & df2$y > unlist(rect_center[,i])[2] - 1.1 & df2$y < unlist(rect_center[,i])[2]+1.1] = 2
}


circ_center = sample(5:45, 2, replace = T)
r = 2 #radius
s = seq(0,r,1/3) #sagitta
l = sqrt(2*(r-s)*r-(r-s)^2) #half chord length
for (i in 1:7){
  df2$z[(df2$x %in% c(circ_center[1] + c(s[i], -s[i]))) & 
         (df2$y > circ_center[2] - l[i+1]) & 
         (df2$y < circ_center[2] + l[i+1])] = 3
}

#Plot truth
p1 = ggplot(data = df2, aes(x = x, y = y))+
  geom_text(aes(color = factor(z), label = rep(c("\u25E3", "\u25E2", "\u25E4", "\u25E5"), each = 1250)), 
            size = 4, family = "Arial Unicode MS")+
  labs(color = "State")+theme_classic()+
  guides(shape = F, size = F, alpha = F)


#Generate simulated expression
set.seed(101)  
mu1 = c(0,0)
mu2 = c(5,0)
mu3 = c(5,5)
Sigma = diag(1, nrow = 2) + 1
Lambda = solve(Sigma)

df2$Y1 = NA
df2$Y2 = NA
df2[df2$z == 1,c("Y1", "Y2")] = rmvnorm(n = sum(df2$z == 1), 
                                      mean = mu1, 
                                      sigma = Sigma)
df2[df2$z == 2,c("Y1", "Y2")] = rmvnorm(n = sum(df2$z == 2), 
                                      mean = mu2, 
                                      sigma = Sigma)
df2[df2$z == 3,c("Y1", "Y2")] = rmvnorm(n = sum(df2$z == 3), 
                                      mean = mu3, 
                                      sigma = Sigma)


#Plots
p2 = ggplot(data = df2, aes(x = x, y = y))+
  geom_text(aes(color = Y1, label = rep(c("\u25E3", "\u25E2", "\u25E4", "\u25E5"), each = 1250)), 
            size = 4, family = "Arial Unicode MS")+
  labs(color = expression(Y[1]))+theme_classic()+
  scale_color_distiller(palette = "Blues", direction = 1, breaks = c(0,5)) +
  guides(shape = F, size = F)

p3 = ggplot(data = df2, aes(x = x, y = y))+
  geom_text(aes(color = Y2, label = rep(c("\u25E3", "\u25E2", "\u25E4", "\u25E5"), each = 1250)), 
            size = 4, family = "Arial Unicode MS")+
  labs(color = expression(Y[2]))+theme_classic()+
  scale_color_distiller(palette = "Blues", direction = 1, breaks = c(0,5)) +
  guides(shape = F, size = F)

p1+p2+p3 + plot_annotation(tag_levels = "A")

df2 %>%
  as.data.frame() %>%
  group_by(j) %>%
  summarise(x = mean(x), y = mean(y),
            Y1= sum(Y1), Y2 = sum(Y2)) %>%
  mutate(j2 = NA) ->
  df2sum
df2sum = as.matrix(df2sum)

#Run mcmc

deconv_sim = run_mcmc_multi(df = df2sum, gamma = 2, q = 3, d = 2,nrep = 1500)
saveRDS(deconv_sim, "data/deconv_sim_3-1.RDS")

df2sum = as.data.frame(df2sum)
df2sum$z_deconv1 = apply(deconv_sim$z[100:1500,], 2, Mode)
df2sum$z_deconv1_alpha = pmax(colMeans(deconv_sim$z[100:1500,]==1),
                           colMeans(deconv_sim$z[100:1500,]==2),
                           colMeans(deconv_sim$z[100:1500,]==3))

pnodeconv = ggplot(as.data.frame(df2sum), aes(x, y)) +
  geom_point(aes(color = factor(z_deconv1), alpha = z_deconv1_alpha), size = 4, shape = 18) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+
  theme_classic()

#Deconvolution
deconv_sim_deconv = run_mcmc_deconv(df = df2sum, gamma = 2, q = 3, d = 2,nrep = 5000, prev = deconv_sim)
saveRDS(deconv_sim_deconv, "data/deconv2_sim_3-2.RDS")

df2 = as.data.frame(df2)
df2$z_deconv2 = apply(deconv_sim_deconv$z[1000:5000,], 2, Mode)
df2$z_deconv2_alpha = pmax(colMeans(deconv_sim_deconv$z[1000:5000,]==1),
                              colMeans(deconv_sim_deconv$z[1000:5000,]==2),
                              colMeans(deconv_sim_deconv$z[1000:5000,]==3))


ggplot(data = as.data.frame(df2), aes(x = x, y = y))+
  geom_text(aes(color = factor(z_deconv2), label = rep(c("\u25E3", "\u25E2", "\u25E4", "\u25E5"), each = 1250),
                alpha = z_deconv2_alpha), size = 4, family = "Arial Unicode MS")+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
  labs(color = "State")+theme_classic()+
  guides(shape = F, size = F, alpha = F)

plot(deconv_sim_deconv$mu[,1], type = "l", ylab = "mu1")
plot(deconv_sim_deconv$mu[,3], type = "l", ylab = "mu1")
plot(sapply(df_sim_Y,"[[", 8), type = "l")
lines(sapply(df_sim_Y,"[[", 8+1250*1), type = "l", col = "red")
lines(sapply(df_sim_Y,"[[", 8+1250*2), type = "l", col = "blue")
lines(sapply(df_sim_Y,"[[", 8+1250*3), type = "l", col = "green")
pdeconv = ggplot(data = as.data.frame(df2), aes(x = x, y = y))+
  geom_text(aes(color = factor(zmode), label = rep(c("\u25E3", "\u25E2", "\u25E4", "\u25E5"), each = 1250),
                alpha = zalpha), size = 4, family = "Arial Unicode MS")+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
  labs(color = "State")+theme_classic()+
  guides(shape = F, size = F, alpha = F)
plot(df_sim_mu[,4], type = "l")
plot(sapply(df_sim_lambda,"[[", 1), type = "l")
plot(sapply(df_sim_Y,"[[", 1), type = "l")
zmode=apply(df_sim_z[100:500,], 2, Mode)
zalpha=pmax(colMeans(df_sim_z[100:500,]==1),
            colMeans(df_sim_z[100:500,]==2),
            colMeans(df_sim_z[100:050,]==3))
plotY1 = ggplot(data = as.data.frame(df2), aes(x = x, y = y))+
  geom_text(aes(color = ifelse(Y1>5,5,ifelse(Y1<0,0, Y1)), label = rep(c("\u25E3", "\u25E2", "\u25E4", "\u25E5"), each = 1250),
                alpha = 1), size = 4, family = "Arial Unicode MS")+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
  labs(color = "Y1")+theme_classic()+
  guides(shape = F, size = F, alpha = F)

plotY2 = ggplot(data = as.data.frame(df2), aes(x = x, y = y))+
  geom_text(aes(color = ifelse(Y2>5,5,ifelse(Y2<0,0, Y2)), label = rep(c("\u25E3", "\u25E2", "\u25E4", "\u25E5"), each = 1250),
                alpha = 1), size = 4, family = "Arial Unicode MS")+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
  labs(color = "Y2")+theme_classic()+
  guides(shape = F, size = F, alpha = F)

plotY1 = ggplot(data = as.data.frame(df2), aes(x = x, y = y))+
  geom_text(aes(color = Y1, label = rep(c("\u25E3", "\u25E2", "\u25E4", "\u25E5"), each = 1250),
                alpha = 1), size = 4, family = "Arial Unicode MS")+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
  labs(color = "Y1")+theme_classic()+
  guides(shape = F, size = F, alpha = F)

plotY2 = ggplot(data = as.data.frame(df2), aes(x = x, y = y))+
  geom_text(aes(color = Y2, label = rep(c("\u25E3", "\u25E2", "\u25E4", "\u25E5"), each = 1250),
                alpha = 1), size = 4, family = "Arial Unicode MS")+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
  labs(color = "Y2")+theme_classic()+
  guides(shape = F, size = F, alpha = F)
