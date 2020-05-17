#library(devtools)
#install_github("LieberInstitute/spatialLIBD")
library(tidyverse)
library(spatialLIBD)
library(scater)
library(scran)
library(mvtnorm)
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
positions = cbind(sce$imagerow, sce$imagecol) #save positions as df
colnames(positions) = c("x", "y") 
xdist = coef(lm(sce$imagecol~sce$col))[2] #x distance between neighbors
ydist = coef(lm(sce$imagerow~sce$row))[2] #y distance between neighbors
dist = xdist + ydist + 0.2
clust1 = cluster(Y= PCs$components, positions = as.matrix(positions), q = 7, init = km[,7], nrep = 10000, gamma = 1.5, dist = dist)
clust1col = apply(clust1$z[9000:10000,], 2, Mode)
clust1alpha = pmax(colMeans(clust1$z[9000:10000,]==1),
                   colMeans(clust1$z[9000:10000,]==2),
                   colMeans(clust1$z[9000:10000,]==3),
                   colMeans(clust1$z[9000:10000,]==4),
                   colMeans(clust1$z[9000:10000,]==5),
                   colMeans(clust1$z[9000:10000,]==6),
                   colMeans(clust1$z[9000:10000,]==7))
clust_out = ggplot(data.frame(positions), aes(y, -x)) +
  geom_text(aes(color = factor(clust1col)), alpha = clust1alpha,
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  #geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

truth_out = ggplot(data.frame(positions), aes(y, -x)) +
  geom_text(aes(color = factor(sce$layer_guess)),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  #geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

clust_out + truth_out

#calculate ARI
mclust::adjustedRandIndex(sce$layer_guess, clust1col)
mclust::adjustedRandIndex(sce$layer_guess, km[,7])
