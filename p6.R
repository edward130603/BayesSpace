source('script.R')
library(tidyverse)
library(scater)
library(patchwork)
data1 <- read.csv("data-raw/data_visium_p6.csv", stringsAsFactors=FALSE)
ggplot(data1,aes(x=spot_col,y=spot_row*sqrt(3),col=CD4)) +
  geom_point(size=5) + 
  coord_fixed() +
  theme_bw()    
positions = cbind(data1$spot_col, data1$spot_row * sqrt(3))
logcounts = t(data1[,-c(1:11,13,51,52)])
sce = SingleCellExperiment(assays = list(logcounts = logcounts))
colnames(sce) = data1$spot
set.seed(100)
sce = runPCA(sce, ncomponents = 4, ntop = 38)

Y_p6 = reducedDim(sce)
set.seed(101)
km_init = kmeans(Y_p6, centers = 5)$cluster

clust_p6 = cluster(Y= Y_p6, positions = positions, q = 5, init = km_init, nrep = 10000, gamma = 2, 
                   dist = (sqrt(3) + 1)*1.05) #dist is l1 distance
plot_km_p6 =ggplot(data1,aes(x=spot_col,y=spot_row*sqrt(3),col=factor(km_init))) +
  geom_point(size=5) + 
  coord_fixed() +
  labs(color = "Cluster", title = "k means")+
  theme_void()
plot_clust_p6 = ggplot(data1,aes(x=spot_col,y=spot_row*sqrt(3),col=factor(clust_p6$labels))) +
  geom_point(size=5) + 
  coord_fixed() +
  labs(color = "Cluster", title = "spatial clustering")+
  theme_void()
plot_km_p6+plot_clust_p6

load("data-raw/p6_Visium_deconv.RData")

deconv2 = deconvolve(Y= deconv1$Y[[1]][1:576,], 
                     positions = deconv1$positions[1:576,], 
                     q = 5, 
                     init = deconv1$z[1,1:576], 
                     nrep = 1000, 
                     gamma = 0, 
                     xdist = 1, ydist = 1, platform = "visium", c = 0.003, jitter_scale = 0.0005)
cols = apply(deconv2$z[500:1000,], 2, Mode)
ggplot(data.frame(deconv2$positions), aes(x = x, y =y ))+
  geom_point(col = factor(cols))+coord_fixed(ratio = sqrt(3))
