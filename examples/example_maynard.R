#library(devtools)
#install_github("LieberInstitute/spatialLIBD")
library(tidyverse)
library(spatialLIBD)
library(scater)
library(scran)
library(extrafont)
library(patchwork)


sce <- fetch_data(type = 'sce')
sce_image_clus(
  sce = sce,
  clustervar = 'layer_guess_reordered',
  sampleid = '151673',
  colors = libd_layer_colors,
  ... = ' DLPFC Human Brain Layers\nMade with github.com/LieberInstitute/spatialLIBD'
)

sce <- fetch_data(type = 'sce')
sce = sce[,sce$sample_name == "151673"]

#Dim reduction and clustering
set.seed(101)
dec <- modelGeneVarByPoisson(sce)
top <- getTopHVGs(dec, prop=0.1)
plot(dec$mean, dec$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)

set.seed(102)
sce <- denoisePCA(sce, technical=dec, subset.row=top)
PCs = getDenoisedPCs(sce, technical=dec, subset.row=top, min.rank = 5, max.rank = 50)
sce <- runTSNE(sce, dimred="PCA")


#Spatial clustering

#choose input
km = sapply(1:10,function(k){kmeans(PCs$components, centers = k)$cluster})
positions = cbind(sce$imagecol, sce$imagerow) #save positions as df
colnames(positions) = c("x", "y") 
xdist = coef(lm(sce$imagecol~sce$col))[2] #x distance between neighbors
ydist = coef(lm(sce$imagerow~sce$row))[2] #y distance between neighbors
dist = xdist + ydist + 0.2
#Normal model
clust1 = cluster(Y= PCs$components, positions = as.matrix(positions), q = 7, init = km[,7], nrep = 10000, gamma = 1.5, dist = dist)
clust1col = apply(clust1$z[9000:10000,], 2, Mode)
clust1alpha = pmax(colMeans(clust1$z[9000:10000,]==1),
                   colMeans(clust1$z[9000:10000,]==2),
                   colMeans(clust1$z[9000:10000,]==3),
                   colMeans(clust1$z[9000:10000,]==4),
                   colMeans(clust1$z[9000:10000,]==5),
                   colMeans(clust1$z[9000:10000,]==6),
                   colMeans(clust1$z[9000:10000,]==7))
clust1_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(clust1col, levels = c(7,1,4,3,6,2,5))), alpha = clust1alpha,
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F, col = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()

truth_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(sce$layer_guess_reordered)),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()

km_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(km[,7], c(7,1,4,3,6,2,5))),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F, col = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()

maynardclust_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(sce$HVG_PCA_spatial, levels = c(5,1,3,2,7,8,6,4))),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F, col = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()

truth_out + clust_out #Truth vs spatial clustering
truth_out + maynardclust_out #Truth vs best clustering implementation from Maynard

#t model
clust2 = cluster(Y= PCs$components, positions = as.matrix(positions), model = "t", q = 7, init = km[,7], nrep = 10000, gamma = 1.5, dist = dist)
clust2col = apply(clust2$z[9000:10000,], 2, Mode)
clust2alpha = pmax(colMeans(clust2$z[9000:10000,]==1),
                   colMeans(clust2$z[9000:10000,]==2),
                   colMeans(clust2$z[9000:10000,]==3),
                   colMeans(clust2$z[9000:10000,]==4),
                   colMeans(clust2$z[9000:10000,]==5),
                   colMeans(clust2$z[9000:10000,]==6),
                   colMeans(clust2$z[9000:10000,]==7))
clust2_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(clust2col, levels = c(7,1,4,3,6,2,5))), alpha = clust2alpha,
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F, col = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()

###variable normal model
clust3 = cluster(Y= PCs$components, positions = as.matrix(positions), precision = "variable", q = 7, init = km[,7], nrep = 10000, gamma = 1.5, dist = dist)
clust3col = apply(clust3$z[9000:10000,], 2, Mode)
clust3alpha = pmax(colMeans(clust3$z[9000:10000,]==1),
                   colMeans(clust3$z[9000:10000,]==2),
                   colMeans(clust3$z[9000:10000,]==3),
                   colMeans(clust3$z[9000:10000,]==4),
                   colMeans(clust3$z[9000:10000,]==5),
                   colMeans(clust3$z[9000:10000,]==6),
                   colMeans(clust3$z[9000:10000,]==7))
clust3_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(clust3col, levels = c(7,1,4,3,6,2,5))), alpha = clust3alpha,
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F, col = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()

###variable t model
clust4 = cluster(Y= PCs$components, positions = as.matrix(positions), model = "t", precision = "variable", q = 7, init = km[,7], nrep = 10000, gamma = 1.5, dist = dist)
clust4col = apply(clust4$z[9000:10000,], 2, Mode)
clust4alpha = pmax(colMeans(clust4$z[9000:10000,]==1),
                   colMeans(clust4$z[9000:10000,]==2),
                   colMeans(clust4$z[9000:10000,]==3),
                   colMeans(clust4$z[9000:10000,]==4),
                   colMeans(clust4$z[9000:10000,]==5),
                   colMeans(clust4$z[9000:10000,]==6),
                   colMeans(clust4$z[9000:10000,]==7))
clust4_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(clust4col, levels = c(7,1,4,3,6,2,5))), alpha = clust4alpha,
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F, col = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()

((truth_out | maynardclust_out | km_out) + plot_layout(guides = "collect"))/(clust1_out | clust2_out|clust3_out|clust4_out) 

#calculate ARI
mclust::adjustedRandIndex(sce$layer_guess, clust1col)
mclust::adjustedRandIndex(sce$layer_guess, clust2col)
mclust::adjustedRandIndex(sce$layer_guess, clust3col)
mclust::adjustedRandIndex(sce$layer_guess, clust4col)
mclust::adjustedRandIndex(clust3col, clust4col)
mclust::adjustedRandIndex(sce$layer_guess, sce$HVG_PCA_spatial)
mclust::adjustedRandIndex(sce$layer_guess, km[,7])
mclust::adjustedRandIndex(clust2col, clust1col)

#Spatial deconvolution
ptm = proc.time()
deconv1 = deconvolve(Y = PCs$components, positions = positions, nrep = 300, gamma = 2, xdist = xdist, ydist = ydist, init = clust1col, q = 7)
proc.time()-ptm #147.75 seconds, ~25x faster than non-Rcpp code
deconv1col = apply(deconv1$z[100:300,], 2, Mode)

