library(tidyverse)
library(SingleCellExperiment)
library(viridis)
library(patchwork)
melanoma1.1 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel1_rep1.RDS")
melanoma1.2 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel1_rep2.RDS")
melanoma2.1 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel2_rep1.RDS")
melanoma2.2 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel2_rep2.RDS")
library(scater)
library(mclust)
source("script.R")
melanoma1.1 = runPCA(melanoma1.1)
melanoma1.2 = runPCA(melanoma1.2)
melanoma2.1 = runPCA(melanoma2.1)
melanoma2.2 = runPCA(melanoma2.2)
positions1.1 = cbind(x = melanoma1.1$col, y = - melanoma1.1$row)
positions1.2 = cbind(x = melanoma1.2$col, y = melanoma1.2$row)
positions2.1 = cbind(y = melanoma2.1$col, x = - melanoma2.1$row)
positions2.2 = cbind(y = - melanoma2.2$col, x = melanoma2.2$row)
Y1.1 = reducedDim(melanoma1.1)[,1:10]
Y1.2 = reducedDim(melanoma1.2)[,1:10]
Y2.1 = reducedDim(melanoma2.1)[,1:10]
Y2.2 = reducedDim(melanoma2.2)[,1:10]

plot1.2 = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = librarySizeFactors(melanoma1.2))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Library Size") +
  scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed()
plot2.1 = ggplot(data.frame(positions2.1), aes(x = x, y = y, fill = librarySizeFactors(melanoma2.1))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Library Size") +
  scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed()

set.seed(100)
km1.1 = kmeans(Y1.1, 3)$cluster
km1.2 = kmeans(Y1.2, 3)$cluster
km2.1 = kmeans(Y2.1, 3)$cluster
km2.2 = kmeans(Y2.2, 3)$cluster
plot_km1.2 = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = factor(km1.2))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("purple", "yellow", "red", "grey"))+
  theme_void()+ coord_fixed()
plot_km2.1 = ggplot(data.frame(positions2.1), aes(x = x, y = y, fill = factor(km2.1))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("orange", "black", "grey"))+
  theme_void()+ coord_fixed()

#Cluster
clust1.2 = cluster(Y= Y1.2, positions = positions1.2, q = 3, init = km1.2, nrep = 10000, gamma = 2, dist = 1)
clust2.1 = cluster(Y= Y2.1, positions = positions2.1, q = 3, init = km2.1, nrep = 10000, gamma = 3, dist = 1)

plot_clust1.2 = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = factor(clust1.2$labels))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("purple", "yellow", "red", "grey"))+
  theme_void()+ coord_fixed()
plot_clust2.1 = ggplot(data.frame(positions2.1), aes(x = x, y = y, fill = factor(clust2.1$labels))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("orange", "black", "grey"))+
  theme_void()+ coord_fixed()

#Deconvolution
deconv1.2 = deconvolve(Y= Y1.2, positions = positions1.2, q = 3, init = clust1.2$labels, nrep = 10000, gamma = 1, xdist = 1, ydist = 1, platform = "ST")
deconv1.2 = deconvolve(Y= test[,1:9], positions = positions1.2, q = 3, init = km1.2, nrep = 50000, gamma = 1, xdist = 1, ydist = 1, platform = "ST")
deconv2.1 = deconvolve(Y= Y2.1, positions = positions2.1, q = 3, init = clust1.2$labels, nrep = 10000, gamma = 2, xdist = 1, ydist = 1, platform = "ST")
deconv1.2 = readRDS("data-raw/deconv1.2_g1_div9.RDS")
deconv1.2 = readRDS("data-raw/deconv1.2_g3.RDS")
deconv1.2_old = readRDS("data-raw/deconv3_km3.RDS")
deconv1.2_oldcol= apply(deconv1.2_old$z[seq(10000,60000,10),],2,Mode)
deconv1.2_alpha = pmax(colMeans(deconv1.2_old$z[seq(10000,60000,10),]==1),
                       colMeans(deconv1.2_old$z[seq(10000,60000,10),]==2),
                       colMeans(deconv1.2_old$z[seq(10000,60000,10),]==3))
deconv1.2colv2 = apply(deconv1.2v2$z[seq(10000,50000,10),], 2, Mode)
deconv1.2colg1 = apply(deconv1.2$obj$z[seq(10000,100000,10),], 2, Mode)
deconv1.2colg2 = apply(deconv1.2$obj$z[seq(10000,100000,10),], 2, Mode)
deconv1.2colg3 = apply(deconv1.2$obj$z[seq(10000,100000,10),], 2, Mode)
deconv1.2_alpha = pmax(colMeans(deconv1.2$obj$z[seq(10000,100000,10),]==1),
                       colMeans(deconv1.2$obj$z[seq(10000,100000,10),]==2),
                       colMeans(deconv1.2$obj$z[seq(10000,100000,10),]==3))
deconv1.2_alpha = pmax(colMeans(deconv1.2$z[9000:10000,]==1),
                       colMeans(deconv1.2$z[9000:10000,]==2),
                       colMeans(deconv1.2$z[9000:10000,]==3))
ggplot(as.data.frame(deconv1.2$obj$positions), 
       aes(x = x, y = y, alpha = 1, fill = factor(deconv1.2colg1))) + 
  geom_point(size = 2.1, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+ 
  labs(x = NULL, y = NULL, fill = "Cluster") +
  guides(alpha = F) +
  scale_fill_manual(values = c("purple", "yellow", "red", "grey"))+
  theme_void()+ coord_fixed()

ggplot(as.data.frame(deconv1.2$obj$positions), 
       aes(x = x, y = y, alpha = deconv1.2_alpha, fill = factor(deconv1.2_oldcol))) + 
  geom_point(size = 2.1, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+ 
  labs(x = NULL, y = NULL, fill = "Cluster") +
  guides(alpha = F) +
  scale_fill_manual(values = c("purple", "yellow", "red", "grey"))+
  theme_void()+ coord_fixed()

test = readRDS("data-raw/df_melanomaA.RDS")
plot(deconv1$mu[,1+15*0], type = "l")
plot(deconv1$plogLik, type = "l")
plot(sapply(deconv1$lambda,"[", 3,3), type = "l")
positions2 = positions[rep(seq_len(nrow(positions)), 9), ] #rbind 9 times
shift = rbind(expand.grid(c(1/3, -1/3, 0), c(1/3,-1/3,0)))
shift_long = shift[rep(seq_len(9), each = nrow(positions)), ]
positions2[,"x"] = positions2[,"x"] + shift_long[,"Var1"]
positions2[,"y"] = positions2[,"y"] + shift_long[,"Var2"]

deconv1_col = apply(deconv1$z[5000:10000,], 2, Mode)
deconv1_alpha = pmax(colMeans(deconv1$z[4000:5000,]==1),
                     colMeans(deconv1$z[4000:5000,]==2),
                     colMeans(deconv1$z[4000:5000,]==3),
                     colMeans(deconv1$z[4000:5000,]==4))
ggplot(positions2, aes(x = y, y = x, fill = factor(deconv1_col), alpha = 1)) + 
  geom_point(size = 2.1, pch = 22)+
  #geom_point(data = as.data.frame(colData(sce_1.2)), aes(x = x, y = -y), size = 5.5, pch = 22, inherit.aes = F, stroke = 1)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("purple", "yellow", "red", "blue"))+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+ 
  guides(alpha= F) +
  #scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()
deconv4_alpha = pmax(colMeans(deconv1$z[40000:50000,]==1),
                     colMeans(deconv1$z[40000:50000,]==2),
                     colMeans(deconv1$z[40000:50000,]==3),
                     colMeans(deconv1$z[40000:50000,]==4))
#deconv1 = run_mcmc_squaredeconv(df, nrep = 100, q = 4, d = 10, seed = 100, prev = clust4)
deconv2 = readRDS("data-raw/deconv3.RDS")
plot(deconv2$mu[,1+0*9], type = "l")
plot(deconv2$mu[,1+2*9], type = "l")
plot(sapply(deconv2$Y,"[[", 2000), type = "l")
deconv2_col = apply(deconv2$z[seq(10000,60000,100),], 2, Mode)
deconv2_alpha = pmax(colMeans(deconv2$z[seq(10000,60000,100),]==1),
                     colMeans(deconv2$z[seq(10000,60000,100),]==2),
                     colMeans(deconv2$z[seq(10000,60000,100),]==3))
Ymeans = Reduce('+', deconv1$Y[250:292])/length(252:292)
ggplot(as.data.frame(positions2), aes(x = y, y = x, 
                               alpha = deconv2_alpha,
                               fill = factor(deconv2_col))) + 
  geom_point(size = 2.1, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+ 
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #scale_fill_viridis_d(option = "A")+
  guides(alpha = F) +
  scale_fill_manual(values = c("purple", "yellow", "red", "blue"))+
  theme_void()+ coord_fixed()
deconv3 = readRDS("data-raw/deconv3_km3.RDS")
deconv3_col = apply(deconv3$z[seq(10000,60000,100),], 2, Mode)
deconv3_alpha = pmax(colMeans(deconv3$z[seq(10000,60000,100),]==1),
                     colMeans(deconv3$z[seq(10000,60000,100),]==2),
                     colMeans(deconv3$z[seq(10000,60000,100),]==3))

deconv4 = readRDS("data-raw/deconv3_mclust3.RDS")
deconv4_col = apply(deconv4$z[seq(10000,60000,100),], 2, Mode)
deconv4_alpha = pmax(colMeans(deconv4$z[seq(10000,60000,100),]==1),
                     colMeans(deconv4$z[seq(10000,60000,100),]==2),
                     colMeans(deconv4$z[seq(10000,60000,100),]==3))

p1 = ggplot(as.data.frame(positions2), aes(x = y, y = x, 
                                           alpha = deconv2_alpha,
                                           fill = factor(deconv2_col))) + 
  geom_point(size = 2.1, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+ 
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #scale_fill_viridis_d(option = "A")+
  guides(alpha = F, fill = F) +
  scale_fill_manual(values = c("purple", "yellow", "red", "blue"))+
  theme_void()+ coord_fixed()
p2 = ggplot(as.data.frame(positions2), aes(x = y, y = x, 
                                           alpha = deconv3_alpha,
                                           fill = factor(deconv3_col))) + 
  geom_point(size = 2.1, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+ 
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #scale_fill_viridis_d(option = "A")+
  guides(alpha = F, fill = F) +
  scale_fill_manual(values = c("purple", "red", "yellow", "blue"))+
  theme_void()+ coord_fixed()
p3 = ggplot(as.data.frame(positions2), aes(x = y, y = x, 
                                           alpha = deconv4_alpha,
                                           fill = factor(deconv4_col))) + 
  geom_point(size = 2.1, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+ 
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #scale_fill_viridis_d(option = "A")+
  guides(alpha = F, fill = F) +
  scale_fill_manual(values = c("yellow", "red", "purple", "blue"))+
  theme_void()+ coord_fixed()
p1+p2+p3+
##mclust
library(mclust)
test = Mclust(PCs$components, G = 1:14, modelNames = c("EEE"))
plot(test, legendArgs = list(x = "topright", ncol = 1, cex = 0.7), what = "BIC")
summary(test, param = T)
ggplot(as.data.frame(colData(sce_A)), aes(x = y, y = x, fill = factor(test$classification))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("plum", "pink", "black", "yellow", "red", "purple", "darkorchid"))+
  #scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()

plot(test)
library(scamp)
test2 = scamp(PCs$components, numberIterations = 10, minimumClusterSize = 10, getVotingHistory = T, clusterOutputString = "test.txt")
ggplot(as.data.frame(PCs$components), aes(x = PCs$components[,1], y = PCs$components[,2], fill = factor(clust3col))) + 
  geom_point(size = 5, pch = 21, alpha = 0.7)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("plum", "pink", "black", "yellow", "red", "purple", "darkorchid"))+
  #scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()
deconv4_col = apply(deconv4$z[seq(10000,60000,100),], 2, Mode)
deconv5_col = apply(deconv4$z[seq(10000,35000,100),], 2, Mode)
deconv6_col = apply(deconv4$z[seq(35000,60000,100),], 2, Mode)
ggplot(as.data.frame(colData(sce_A)), aes(x = y, y = x, fill = PCs$components[,2])) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "PC1") +
  scale_fill_distiller(palette = "RdBu", limit = max(abs(PCs$components[,2]))*c(-1,1))+
  #scale_fill_manual(values = c("plum", "pink", "black", "yellow", "red", "purple", "darkorchid"))+
  #scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()
ggplot(as.data.frame(colData(sce_A)), aes(x = y, y = x, fill = ifelse(PCs$components[,2]>15,15,PCs$components[,2]))) + 
  geom_point(size = 7, pch = 22)+
  scale_fill_distiller(palette = "RdBu", limit =  15*c(-1,1))+
  labs(x = NULL, y = NULL, fill = "PC2") +
  #scale_fill_manual(values = c("plum", "pink", "black", "yellow", "red", "purple", "darkorchid"))+
  #scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()
ggplot(as.data.frame(colData(sce_A)), aes(x = y, y = x, fill = PCs$components[,2])) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "PC1") +
  scale_fill_distiller(palette = "RdBu", limit = max(abs(PCs$components[,2]))*c(-1,1))+
  #scale_fill_manual(values = c("plum", "pink", "black", "yellow", "red", "purple", "darkorchid"))+
  #scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()
simdata = as.data.frame(expand.grid(1:60,1:60)/3)
simdata$z = 1
ggplot(as.data.frame(expand.grid(1:20,1:20)), aes(x = Var1, y = Var2)) + 
  geom_point(size = 12, pch = 22)+
  labs(x = NULL, y = NULL, fill = "PC1") +
  guides(fill = F)+
  theme_void()+ coord_fixed()
cormat = cov2cor(solve(Reduce("+", clust3$lambda)/ length(clust3$lambda)))
colnames(cormat) = 1:9
as_tibble(cormat) %>% mutate(PC1 = 1:9) %>% gather("PC2", "Corr", `1`:`9`) %>% mutate(PC1 = as.character(PC1)) %>%
  ggplot(aes(x = PC1, y = PC2, fill = Corr))+
  geom_tile()+scale_fill_distiller(palette = "RdBu", limit = c(-1,1), direction = 1) +labs(x = NULL, y = NULL)
