#library(devtools)
#install_github("LieberInstitute/spatialLIBD")
library(tidyverse)
library(spatialLIBD)
library(scater)
library(scran)
library(mvtnorm)
library(extrafont)


sce <- fetch_data(type = 'sce')
sce_image_clus(
  sce = sce,
  clustervar = 'layer_guess_reordered',
  sampleid = '151673',
  colors = libd_layer_colors,
  ... = ' DLPFC Human Brain Layers\nMade with github.com/LieberInstitute/spatialLIBD'
)

sce_151673 = sce[,sce$sample_name == "151673"]
ggplot(as.data.frame(colData(sce_151673)), aes(imagecol, -imagerow, color = layer_guess_reordered)) +
  geom_text(label = "\u2B22", size = 6, family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = NULL, x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  theme_classic() + coord_fixed()

ggplot(as.data.frame(colData(sce_151673)), aes(imagecol, -imagerow, color = reducedDim(sce_151673v2,"TSNE")[,2])) +
  geom_text(label = "\u2B22", size = 6, family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = NULL, x = NULL, y = NULL) +
  scale_color_viridis(option = "A")+
  theme_classic() + coord_fixed()
coef(lm(sce_151673$imagecol~sce_151673$col))
coef(lm(sce_151673$imagerow~sce_151673$row))

#rerun dim reduc
sce_151673v2 = sce_151673
set.seed(101)
dec <- modelGeneVarByPoisson(sce_151673v2)
top <- getTopHVGs(dec, prop=0.1)
plot(dec$mean, dec$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)

#Dim reduction and clustering
set.seed(102)
sce_151673v2 <- denoisePCA(sce_151673v2, technical=dec, subset.row=top)
PCs = getDenoisedPCs(sce_151673v2, technical=dec, subset.row=top, min.rank = 5, max.rank = 50)
sce_151673v2 <- runTSNE(sce_151673v2, dimred="PCA")
sce_151673v2 <- runTSNE(sce_151673v2, dimred="PCA", ncomponents = 3, name = "TSNE3")

snn.gr <- buildSNNGraph(sce_A, use.dimred="PCA", k=40)
sce_A$cluster <- factor(igraph::cluster_walktrap(snn.gr)$membership)
table(sce_A$cluster)
sce_A$cluster = factor(sce_A$cluster, levels = c(1,2,4,5,3))

#mcmc method

#choose input
df_151673 = data.frame(reducedDim(sce_151673, "PCA"))[,1:9]
df_151673 = data.frame(reducedDim(sce_151673, "TSNE_perplexity5"))
df_151673 = data.frame(reducedDim(sce_151673, "UMAP_neighbors15"))
df_151673 = data.frame(reducedDim(sce_151673v2, "PCA"))[,1:9]
df_151673 = data.frame(reducedDim(sce_151673v2, "TSNE"))[,1:2]
df_151673 = data.frame(reducedDim(sce_151673v2, "TSNE3"))[,1:3]


km5 = kmeans(reducedDim(sce_151673v2, "PCA"), centers = 5)$cluster #k means initialization
positions = cbind(sce_151673$imagerow, sce_151673$imagecol) #save positions as df
colnames(positions) = c("x", "y") 
xdist = coef(lm(sce_151673$imagecol~sce_151673$col))[2] #x distance between neighbors
ydist = coef(lm(sce_151673$imagerow~sce_151673$row))[2] #y distance between neighbors
dist = xdist + ydist + 0.2
test = cluster(Y= PCs$components, positions = as.matrix(positions), q = 5, init = km5, nrep = 1000, gamma = 1.5, dist = dist)

colnames(df_151673) = paste0("Y", 1:9)
df_151673$kmeans = kmeans(reducedDim(sce_151673, "PCA"), centers = 7)$cluster
df_151673$kmeans8 = kmeans(reducedDim(sce_151673, "PCA"), centers = 8)$cluster
df_151673$kmeans6 = kmeans(reducedDim(sce_151673, "PCA"), centers = 6)$cluster
df_151673$kmeans5 = kmeans(reducedDim(sce_151673, "PCA"), centers = 5)$cluster
df_151673$kmeans4 = kmeans(reducedDim(sce_151673, "PCA"), centers = 4)$cluster
df_151673$x = sce_151673$imagecol
df_151673$y = sce_151673$imagerow
df_151673$j = 1:nrow(df_151673)
num_neighbors = sapply(sapply(1:nrow(df_151673), function(x){df_151673[(abs(df_151673[,"x"] -df_151673[x,"x"]) + abs(df_151673[,"y"] - df_151673[x,"y"])) <= 9 &
                                                                         (abs(df_151673[,"x"] -df_151673[x,"x"]) + abs(df_151673[,"y"] - df_151673[x,"y"])) > 0,"j"]}), length)
df_151673 = df_151673[num_neighbors>0,]
clust1_151673 = run_mcmc_multi(df = df_151673 %>% select(-(Y6:Y9)), gamma = 2, q = 7, d = 5,nrep = 2000,
                               init = df_151673$kmeans)
df_151673$z = apply(clust1_151673$z[-(1:500),], 2, Mode)
df_151673$z_gamma = pmax(colMeans(clust1_151673$z[-(1:500),]==1),
                              colMeans(clust1_151673$z[-(1:500),]==2),
                              colMeans(clust1_151673$z[-(1:500),]==3),
                              colMeans(clust1_151673$z[-(1:500),]==4),
                              colMeans(clust1_151673$z[-(1:500),]==5),
                              colMeans(clust1_151673$z[-(1:500),]==6),
                              colMeans(clust1_151673$z[-(1:500),]==7))
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z, levels = c(5,6,2,3,7,1,4))),#, alpha = z_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
    labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

clust2_151673 = run_mcmc_multi(df = df_151673 %>% select(-(Y6:Y9)), gamma = 4, q = 7, d = 5,nrep = 2000,
                               init = df_151673$kmeans)
df_151673$z2 = apply(clust2_151673$z[-(1:500),], 2, Mode)
df_151673$z2_gamma = pmax(colMeans(clust2_151673$z[-(1:500),]==1),
                         colMeans(clust2_151673$z[-(1:500),]==2),
                         colMeans(clust2_151673$z[-(1:500),]==3),
                         colMeans(clust2_151673$z[-(1:500),]==4),
                         colMeans(clust2_151673$z[-(1:500),]==5),
                         colMeans(clust2_151673$z[-(1:500),]==6),
                         colMeans(clust2_151673$z[-(1:500),]==7))
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z2, levels = c(5,6,2,3,7,1,4))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

clust3_151673 = run_mcmc_multi(df = df_151673, gamma = 2, q = 7, d = 9,nrep = 2000, 
                               init = df_151673$kmeans)
df_151673$z3 = apply(clust3_151673$z[500:2000,], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z3, levels = c(5,6,2,3,7,1,4))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()
clust4_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 7, d = 9,nrep = 2000,
                               init = df_151673$kmeans)
df_151673$z4 = apply(clust4_151673$z[500:2000,], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z4, levels = c(5,6,2,3,7,1,4))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

df_151673_v2 = df_151673
df_151673_v2$Y1 = reducedDim(sce_151673,"TSNE_perplexity5")[num_neighbors>0,1]
df_151673_v2$Y2 = reducedDim(sce_151673,"TSNE_perplexity5")[num_neighbors>0,2]
clust5_151673 = run_mcmc_multi(df = df_151673_v2, gamma = 4, q = 7, d = 2,nrep = 2000,
                               init = df_151673_v2$kmeans)
df_151673$z5 = apply(clust5_151673$z[1000:2000,], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z5, levels = c(5,6,2,3,7,1,4))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

clust6_151673 = run_mcmc_multi(df = df_151673_v2, gamma = 6, q = 7, d = 2,nrep = 2000,
                               init = df_151673_v2$kmeans)
df_151673$z6 = apply(clust6_151673$z[500:2000,], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z6, levels = c(5,6,2,3,7,1,4))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

clust7_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 8, d = 9,nrep = 2000,
                               init = df_151673$kmeans8)
df_151673$z7 = apply(clust7_151673$z[500:2000,], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z7, levels = c(3,4,7,8,6,5,1,2))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

clust8_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 6, d = 9,nrep = 2000,
                               init = df_151673$kmeans6)
df_151673$z8 = apply(clust8_151673$z[500:2000,], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z8, levels = c(6,1,2,5,3,4))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

clust9_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 6, d = 2,nrep = 2000,
                               init = df_151673$kmeans6)
df_151673$z9 = apply(clust9_151673$z[500:2000,], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z9, levels = c(6,1,2,5,3,4))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

clust10_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 7, d = 2,nrep = 2000,
                               init = df_151673$kmeans)
df_151673$z10 = apply(clust10_151673$z[500:2000,], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z10 , levels = c(4,6,1,2,5,7,3))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

clust11_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 6, d = 2,nrep = 2000,
                                init = df_151673$kmeans6)
df_151673$z11 = apply(clust11_151673$z[500:2000,], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z11 , levels = c(4,6,1,2,5,3))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

truth_init = as.numeric(sce_151673$layer_guess[num_neighbors>0])
truth_init[is.na(truth_init)] = 6
clust12_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 7, d = 2,nrep = 2000,
                                init = truth_init)

df_151673$z12 = apply(clust12_151673$z[500:2000,], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z12, levels = c(4,6,1,2,5,3,7))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

clust13_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 7, d = 2,nrep = 20000,
                                init = truth_init)

df_151673$z13 = apply(clust13_151673$z[10000:20000,], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z13, levels = c(4,6,1,2,5,3,7))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

clust14_151673 = run_mcmc_multi(df = df_151673 %>% select(-(Y6:Y9)), gamma = 4, q = 7, d = 5,nrep = 2000,
                               init = df_151673$kmeans)
df_151673$z14 = apply(clust14_151673$z[-(1:500),], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z14, levels = c(5,6,2,3,7,1,4))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  theme_classic() + coord_fixed()


clust15_151673 = run_mcmc_multi(df = df_151673 %>% select(-(Y6:Y9)), gamma = 4, q = 6, d = 5,nrep = 2000,
                                init = df_151673$kmeans6)
df_151673$z15 = apply(clust15_151673$z[-(1:500),], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z15, levels = c(6,1,5,4,3,2))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  theme_classic() + coord_fixed()

clust16_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 6, d = 9,nrep = 2000,
                                init = df_151673$kmeans6)
df_151673$z16 = apply(clust16_151673$z[-(1:500),], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z16, levels = c(6,1,5,4,3,2))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  theme_classic() + coord_fixed()
z16means = colMeans(clust16_151673$mu[-(1:500),])
z16means1 = z16means[seq(1,46,9)]

ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = df_151673$Y1-z16means1[df_151673$z16]), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "PC1 residuals", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_distiller(palette = "RdBu")+
  theme_classic() + coord_fixed()



clust17_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 7, d = 9,nrep = 2000,
                                init = df_151673$kmeans)
df_151673$z17 = apply(clust17_151673$z[-(1:700),], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z17, levels = c(1,6,3,7,5,4,2))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  theme_classic() + coord_fixed()

clust18_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 7, d = 2,nrep = 2000,
                                init = df_151673$kmeans)
df_151673$z18 = apply(clust18_151673$z[-(1:500),], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z18, levels = c(1,6,3,7,5,4,2))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  theme_classic() + coord_fixed()

clust19_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 6, d = 2,nrep = 2000,
                                init = df_151673$kmeans6)
df_151673$z19 = apply(clust19_151673$z[-(1:500),], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z19, levels = c(1,6,3,5,4,2))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  theme_classic() + coord_fixed()

clust20_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 6, d = 3,nrep = 2000,
                                init = df_151673$kmeans6)
df_151673$z20 = apply(clust20_151673$z[-(1:500),], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z20, levels = c(1,6,3,5,4,2))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  theme_classic() + coord_fixed()

clust21_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 7, d = 3,nrep = 2000,
                                init = df_151673$kmeans)
df_151673$z21 = apply(clust21_151673$z[-(1:1000),], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z21, levels = c(1,6,3,5,4,2,7))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  theme_classic() + coord_fixed()

clust22_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 5, d = 9,nrep = 2000,
                                init = df_151673$kmeans5)
df_151673$z22 = apply(clust22_151673$z[-(1:500),], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z22, levels = c(1,5,4,3,2))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  theme_classic() + coord_fixed()

clust23_151673 = run_mcmc_multi(df = df_151673, gamma = 4, q = 4, d = 9,nrep = 2000,
                                init = df_151673$kmeans4)
df_151673$z23 = apply(clust23_151673$z[-(1:500),], 2, Mode)
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(z23, levels = c(1,4,3,2))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  theme_classic() + coord_fixed()


#spatialLIBD clustering methods 
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(sce_151673$HVG_PCA_spatial[num_neighbors>0], 
                               levels = c(5,1,3,2,7,8,6,4))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(sce_151673$markers_PCA_spatial[num_neighbors>0], 
                               levels = c(4,6,1,5,7,2,3,8))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()
ggplot(df_151673, aes(x, -y)) +
  geom_text(aes(color = factor(sce_151673$pseudobulk_PCA_spatial[num_neighbors>0], 
                               levels = c(8,2,3,1,4,7,6,5))), #alpha = z2_gamma),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A")+
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

#calculate ARI
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], sce_151673$HVG_PCA_spatial[num_neighbors>0])
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], sce_151673$markers_PCA_spatial[num_neighbors>0])
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], sce_151673$pseudobulk_PCA_spatial[num_neighbors>0])
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z2)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z3)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z4)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z5)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z6)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z7)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z8)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z9)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z10)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z11)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z12)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z13)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z14)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z15)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z16)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z17)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z18)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z19)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z20)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z21)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z22)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$z23)
mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], df_151673$kmeans8)
#mclust::adjustedRandIndex(sce_151673$layer_guess[num_neighbors>0], kmeans(reducedDim(sce_151673v2, "PCA"), centers = 6)$cluster[num_neighbors>0])
mclust::adjustedRandIndex(df_151673$z12, df_151673$z13)
mclust::adjustedRandIndex(df_151673$z10, df_151673$z13)
mclust::adjustedRandIndex(df_151673$z10, df_151673$z12)

results = read.csv("output/results.csv")
ggplot(results, aes(x = Input, y = ARI, color = factor(Dimensions))) +
  geom_point(size = 2) +
  labs(color = "Dimensions", linetype = "References")+
  geom_hline(aes(yintercept = 0.274, linetype = "HVG"))+
  geom_hline(aes(yintercept = 0.197, linetype = "Markers"))+
  geom_hline(aes(yintercept = 0.380, linetype = "Semi-supervised"))+
  theme_light()

ggplot(results, aes(x = Input, y = ARI, color = factor(Clusters))) +
  geom_point(size = 2) +
  labs(color = "Clusters", linetype = "References")+
  geom_hline(aes(yintercept = 0.274, linetype = "HVG"))+
  geom_hline(aes(yintercept = 0.197, linetype = "Markers"))+
  geom_hline(aes(yintercept = 0.380, linetype = "Semi-supervised"))+  theme_light()

ggplot(results, aes(x = Input, y = ARI, color = (Notes==""))) +
  geom_point(size = 2) +
  labs(color = "Single sample?", linetype = "References")+
  geom_hline(aes(yintercept = 0.274, linetype = "HVG"))+
  geom_hline(aes(yintercept = 0.197, linetype = "Markers"))+
  geom_hline(aes(yintercept = 0.380, linetype = "Semi-supervised"))+  theme_light()

# plotdim + plotclust + plotsingle +
#   plot_layout(guides = "collect") & theme(legend.position = 'bottom')
#plot t-SNE
tsne1 = ggplot(as.data.frame(colData(sce_151673)), aes(imagecol, -imagerow, color = reducedDim(sce_151673,"TSNE_perplexity5")[,1])) +
  geom_text(label = "\u2B22", size = 4, family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 4, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = NULL, x = NULL, y = NULL) +
  scale_color_viridis(option = "A")+
  theme_classic() + coord_fixed()
tsne2 = ggplot(as.data.frame(colData(sce_151673)), aes(imagecol, -imagerow, color = reducedDim(sce_151673,"TSNE_perplexity5")[,2])) +
  geom_text(label = "\u2B22", size = 4, family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 4, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = NULL, x = NULL, y = NULL) +
  scale_color_viridis(option = "A")+
  theme_classic() + coord_fixed()
#plot PCs
pc1 = ggplot(as.data.frame(colData(sce_151673)), aes(imagecol, -imagerow, color = reducedDim(sce_151673)[,"PC1"])) +
  geom_text(label = "\u2B22", size = 2.7, family = "Lucida Sans Unicode") +
  labs(color = "PC1", x = NULL, y = NULL) +
  scale_color_viridis(option = "A")+
  theme_classic() + coord_fixed()
pc2 = ggplot(as.data.frame(colData(sce_151673)), aes(imagecol, -imagerow, color = reducedDim(sce_151673)[,"PC2"])) +
  geom_text(label = "\u2B22", size = 2.7, family = "Lucida Sans Unicode") +
  labs(color = "PC2", x = NULL, y = NULL) +
  scale_color_viridis(option = "A")+
  theme_classic() + coord_fixed()
pc3 = ggplot(as.data.frame(colData(sce_151673)), aes(imagecol, -imagerow, color = reducedDim(sce_151673)[,"PC3"])) +
  geom_text(label = "\u2B22", size = 2.7, family = "Lucida Sans Unicode") +
  labs(color = "PC3", x = NULL, y = NULL) +
  scale_color_viridis(option = "A")+
  theme_classic() + coord_fixed()
pc4 = ggplot(as.data.frame(colData(sce_151673)), aes(imagecol, -imagerow, color = reducedDim(sce_151673)[,"PC4"])) +
  geom_text(label = "\u2B22", size = 2.7, family = "Lucida Sans Unicode") +
  labs(color = "PC4", x = NULL, y = NULL) +
  scale_color_viridis(option = "A")+
  theme_classic() + coord_fixed()
pc1+pc2+pc3+pc4
