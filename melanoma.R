library(tidyverse)
library(SingleCellExperiment)
library(viridis)
library(patchwork)
melanoma1.1 = readRDS("ST_mel1_rep1.RDS")
melanoma1.2 = readRDS("ST_mel1_rep2.RDS")
melanoma2.1 = readRDS("ST_mel2_rep1.RDS")
melanoma2.2 = readRDS("ST_mel2_rep2.RDS")
# melanoma1.1 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel1_rep1.RDS")
# melanoma1.2 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel1_rep2.RDS")
# melanoma2.1 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel2_rep1.RDS")
# melanoma2.2 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel2_rep2.RDS")
library(scater)
library(mclust)
source("script.R")
set.seed(100)
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
  #scale_fill_manual(values = c("orange", "black", "grey", "red", "green"))+
  theme_void()+ coord_fixed()

#Cluster
# clust1.2 = cluster(Y= Y1.2, positions = positions1.2, q = 3, init = km1.2, nrep = 1000, gamma = 2, dist = 1)
# clust2.1 = cluster(Y= Y2.1, positions = positions2.1, q = 10, init = km2.1, nrep = 5000, gamma = 3, dist = 1, model = "t")

# plot_clust1.2 = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = factor(clust1.2$labels))) + 
#   geom_point(size = 7, pch = 22)+
#   labs(x = NULL, y = NULL, fill = "Cluster") +
#   scale_fill_manual(values = c("purple", "yellow", "red", "grey"))+
#   theme_void()+ coord_fixed()
# plot_clust2.1 = ggplot(data.frame(positions2.1), aes(x = x, y = y, fill = factor(clust2.1$labels))) + 
#   geom_point(size = 7, pch = 22)+
#   labs(x = NULL, y = NULL, fill = "Cluster") +
#   #scale_fill_manual(values = c("orange", "black", "grey"))+
#   theme_void()+ coord_fixed()

#Deconvolution
deconv1.2 = deconvolve(Y= Y1.2, positions = positions1.2, q = 3, init = km1.2, nrep = 100000, gamma = 1, xdist = 1, ydist = 1, platform = "ST", c = 0.01)
deconv1.2col = apply(deconv1.2$z[seq(10000,100000,10),], 2, Mode)
deconv1.2_alpha = pmax(colMeans(deconv1.2$z[seq(10000,100000,10),]==1),
                       colMeans(deconv1.2$z[seq(10000,100000,10),]==2),
                       colMeans(deconv1.2$z[seq(10000,100000,10),]==3))
saveRDS(list(obj = deconv1.2, cols = deconv1.2col, alpha = deconv1.2_alpha), "deconv1.2_g1_c0.01.RDS")

deconv2.1 = deconvolve(Y= Y2.1, positions = positions2.1, q = 3, init = km2.1, nrep = 100000, gamma = 1, xdist = 1, ydist = 1, platform = "ST", c = 0.01)
deconv2.1col = apply(deconv2.1$z[seq(10000,100000,10),], 2, Mode)
deconv2.1_alpha = pmax(colMeans(deconv2.1$z[seq(10000,100000,10),]==1),
                       colMeans(deconv2.1$z[seq(10000,100000,10),]==2),
                       colMeans(deconv2.1$z[seq(10000,100000,10),]==3))
saveRDS(list(obj = deconv2.1, cols = deconv2.1col, alpha = deconv2.1_alpha), "deconv2.1_g1_c0.01.RDS")

