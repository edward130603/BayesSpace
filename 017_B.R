library(tidyverse)
library(DropletUtils)
library(scater)
library(ggbeeswarm)
library(scran)
library(viridis)
library(extrafont)
library(patchwork)

#Make sce
sce_17B = read10xCounts("data-raw/017_B")
sce_17B$Sample = "17B"
rownames(sce_17B) = uniquifyFeatureNames(rowData(sce_17B)$ID, rowData(sce_17B)$Symbol)
colnames(sce_17B) = paste0(sce_17B$Sample, '.', sce_17B$Barcode)

#Add position
pos = read.csv("data-raw/017_B/tissue_positions_list.csv", header=FALSE)
colnames(pos) = c("Barcode", "tissue", "Y1", "X1", "Y2", "X2")
colData(sce_17B) = merge(colData(sce_17B), pos, by = "Barcode")

#Preprocess
is.mito = grep("_MT-", rownames(sce_17B))
is.RBC = grep("_HB", rownames(sce_17B))
df_qc = perCellQCMetrics(sce_17B, subsets = list(Mito = is.mito,
                                          RBC = is.RBC))
qc = quickPerCellQC(df_qc, percent_subsets = "subsets_Mito_percent")
sce_17B$discard = qc$discard
coord1=ggplot(as.data.frame(colData(sce_17B)), aes(x = X2, y = Y2, col = clusters)) + 
  geom_text(label = "\u2B22", size = 6, family = "Lucida Sans Unicode")+
  theme_classic()+coord_fixed()+guides(col = F)
sce_17B = sce_17B[,!qc$discard]

#Normalization
set.seed(100)
clusters = quickCluster(sce_17B)
sce_17B = computeSumFactors(sce_17B, clusters = clusters)
sce_17B = logNormCounts(sce_17B)
plot(librarySizeFactors(sce_17B), sizeFactors(sce_17B), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", log="xy", 
     col = ifelse(sce_17B$discard, "red", "black"))

set.seed(101)
dec <- modelGeneVarByPoisson(sce_17B)
top <- getTopHVGs(dec, prop=0.1)
plot(dec$mean, dec$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)

#Dim reduction and clustering
set.seed(102)
sce_17B <- denoisePCA(sce_17B, technical=dec, subset.row=top)
PCs = getDenoisedPCs(sce_17B, technical=dec, subset.row=top, min.rank = 5, max.rank = 50)
sce_17B <- runTSNE(sce_17B, dimred="PCA")

snn.gr <- buildSNNGraph(sce_17B, use.dimred="PCA", k=40)
sce_17B$cluster <- factor(igraph::cluster_walktrap(snn.gr)$membership)
table(sce_17B$cluster)

tsne1 = ggplot(as.data.frame(reducedDim(sce_17B, "TSNE")), aes(x = V1, y = V2, col = clusters))+
  geom_point(size = 2) +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()
  
spatial1 = ggplot(as.data.frame(colData(sce_17B)), aes(x = X2, y = Y2, col = clusters)) + 
  geom_text(label = "\u2B22", 
            size = 7, family = "Lucida Sans Unicode")+
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()+
  theme(legend.position = c(0.05,0.8)) + coord_fixed()

tsne1 + spatial1 +
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "A")
  

#Distribution of principal components
pc1_spatial = ggplot(as.data.frame(colData(sce_17B)), aes(x = X2, y = Y2, col = reducedDim(sce_17B, "PCA")[,"PC1"])) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "PC1")+scale_color_viridis(option = "A")+theme_classic() + coord_fixed()
pc2_spatial = ggplot(as.data.frame(colData(sce_17B)), aes(x = X2, y = Y2, col = reducedDim(sce_17B, "PCA")[,"PC2"])) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "PC2")+scale_color_viridis(option = "A")+theme_classic()  + coord_fixed()
pc3_spatial = ggplot(as.data.frame(colData(sce_17B)), aes(x = X2, y = Y2, col = reducedDim(sce_17B, "PCA")[,"PC3"])) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "PC3")+scale_color_viridis(option = "A")+theme_classic()  + coord_fixed()
pc4_spatial = ggplot(as.data.frame(colData(sce_17B)), aes(x = X2, y = Y2, col = reducedDim(sce_17B, "PCA")[,"PC4"])) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "PC4")+scale_color_viridis(option = "A")+theme_classic()  + coord_fixed()
pc1_spatial + pc2_spatial + pc3_spatial + pc4_spatial

#Run mcmc
df_17B = data.frame(reducedDim(sce_17B, "PCA"))[,1:4]
colnames(df_17B) = paste0("Y", 1:4)
df_17B$x = sce_17B$X2
df_17B$y = sce_17B$Y2
df_17B$j = 1:nrow(df_17B)
num_neighbors = sapply(sapply(1:nrow(df_17B), function(x){df_17B[(abs(df_17B[,"x"] -df_17B[x,"x"]) + abs(df_17B[,"y"] - df_17B[x,"y"])) <= 40 &
                                                               (abs(df_17B[,"x"] -df_17B[x,"x"]) + abs(df_17B[,"y"] - df_17B[x,"y"])) > 0,"j"]}), length)
df_17B = df_17B[num_neighbors >0, ] #all have at least 1 neighbor (in fact at least 3)

df_17B_gamma2_q5 = run_mcmc_multi(df = df_17B, gamma = 2, q = 5, d = 4,nrep = 2000)
saveRDS(df_17B_gamma2_q5, "mcmc_d5_q5_A.RDS")
df_17B_gamma2_q5 = readRDS("data/mcmc_d5_q5_A.RDS")

#5 states, gamma = 2
df_17B$z_gamma2_q5 = apply(df_17B_gamma2_q5$z[-(1:500),], 2, Mode)
df_17B$z_gamma2_q5_alpha = pmax(colMeans(df_17B_gamma2_q5$z[-(1:500),]==1),
                               colMeans(df_17B_gamma2_q5$z[-(1:500),]==2),
                               colMeans(df_17B_gamma2_q5$z[-(1:500),]==3),
                               colMeans(df_17B_gamma2_q5$z[-(1:500),]==4),
                               colMeans(df_17B_gamma2_q5$z[-(1:500),]==5))

ggplot(df_17B, aes(x, y)) +
  geom_text(aes(color = factor(z_gamma2_q5), alpha = z_gamma2_q5_alpha),
            size = 9, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

ggplot(as.data.frame(colData(sce_17B)), aes(x = X1, y = Y1, col = sce_17B$cluster)) + 
  geom_text(label = "\u2B22", size = 13, family = "Lucida Sans Unicode") +
  theme_classic() + coord_fixed() + labs(col = "State")

ggplot(df_17B, aes(x, y)) +
  geom_point(aes(color = factor(z_gamma2_q5), alpha = z_gamma2_q5_alpha), size = 4, shape = 18) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
  theme_classic()
ggplot(as.data.frame(reducedDim(sce_17B, "TSNE")[num_neighbors>0,]), aes(x = V1, y = V2, col = factor(df_17B$z_gamma2_q5)))+
  geom_point() +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()


ggplot(df_17B, aes(x, y, color = logcounts(sce_17B)[grep("Human_GRCh38_CD3D", rownames(sce_17B)),num_neighbors>0] +
                   logcounts(sce_17B)[grep("Human_GRCh38_CD3E$", rownames(sce_17B)),num_neighbors>0]+
                   logcounts(sce_17B)[grep("Human_GRCh38_CD3G", rownames(sce_17B)),num_neighbors>0] )) +
  geom_text(label = "\u2B22", size = 9, family = "Lucida Sans Unicode") +
  labs(color = "CD3", x = NULL, y = NULL) +
  scale_color_viridis(option = "A")+
  theme_classic() + coord_fixed()

ggplot(df_17B, aes(x, -y, color = logcounts(sce_17B)[grep("CHGA", rownames(sce_17B)),num_neighbors>0])) +
  geom_text(label = "\u2B22", size = 9, family = "Lucida Sans Unicode") +
  labs(color = "CHGA", x = NULL, y = NULL) +
  scale_color_viridis(option = "A")+
  theme_classic() + coord_fixed()

ggplot(df_17B, aes(x, y, color = logcounts(sce_17B)[grep("MCPyV________MCPyV_gp3", rownames(sce_17B)),num_neighbors>0])) +
  geom_text(label = "\u2B22", size = 9, family = "Lucida Sans Unicode") +
  labs(color = "CHGA", x = NULL, y = NULL) +
  scale_color_viridis()+
  theme_classic() + coord_fixed()

ggplot(as.data.frame(reducedDim(sce_17B, "TSNE")[num_neighbors>0,]), aes(x = V1, y = V2, col = logcounts(sce_17B)[grep("KRT20", rownames(sce_17B)),]))+
  geom_point() +
  scale_color_viridis() +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()

##Deconvolution
deconv1 = run_mcmc_hexdeconv(df = df_17B, gamma = 2, q = 5, d = 5,nrep = 2000, prev = df_17B_gamma2_q5)
save(deconv1, file = "data/deconv1.RDS")
df_17B2 = df_17B[rep(seq_len(nrow(df_17B)), 7), ] #rbind 7 times
df_17B2$j2 = 1:(7*nrow(df_17B))

shift = rbind(expand.grid(c(1/2, -1/2), c(1/4,-1/4)), expand.grid(0, c(0, 1/2, -1/2)))
shift_long = shift[rep(seq_len(7), each = nrow(df_17B)), ]

df_17B2$x = df_17B2$x + shift_long$Var1
df_17B2$y = df_17B2$y + shift_long$Var2

df_17B2$z_deconv1 = apply(deconv1$z[500:2000,], 2, Mode)
df_17B2$z_deconv1_alpha = pmax(colMeans(deconv1$z[500:2000,]==1),
                           colMeans(deconv1$z[500:2000,]==2),
                           colMeans(deconv1$z[500:2000,]==3),
                           colMeans(deconv1$z[500:2000,]==4),
                           colMeans(deconv1$z[500:2000,]==5))


plot(deconv1$mu[,1], type = "l", ylab = "mu")
plot(deconv1$mu[,6], type = "l", ylab = "mu")
plot(sapply(deconv1$lambda,"[[",1), type = "l", ylab = "lambda")
plot(sapply(deconv1$lambda,"[[",2), type = "l", ylab = "lambda")


ggplot(data = df_17B2, aes(x = x, y = y))+
  geom_text(aes(color = factor(z_deconv1), alpha = z_deconv1_alpha), label = "\u2B23",
            size = 4, family = "Lucida Sans Unicode")+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+
  labs(color = "State")+theme_classic()+
  guides(shape = F, size = F, alpha = F)+
  geom_text(data = df_17B, aes(x = x+0.03, y = y+0.2), size = 13, label = "\u2B21", alpha = 0.3, family = "Lucida Sans Unicode")+
  coord_fixed()

##New coordinates
#cluster
df_17Bcoord2 = df_17B  
df_17Bcoord2$x = sce_17B$X2
df_17Bcoord2$y = sce_17B$Y2
# index = 2
# test = rep(1, n)
# test[index] = 2
# test[df_j[[index]]] = 3
# ggplot(df_17Bcoord2, aes(x, y, color = factor(test))) +
#   geom_text(label = "\u2B22", size = 10, family = "Lucida Sans Unicode") +
#   labs(x = NULL, y = NULL) +
#   scale_color_viridis_d()+
#   theme_classic() + coord_fixed()
dfv2_A_gamma2_q5 = run_mcmc_multi(df = df_17Bcoord2, gamma = 2, q = 5, d = 5,nrep = 2000)
df_17Bcoord2$z_gamma2_q5 = apply(dfv2_A_gamma2_q5$z[-(1:500),], 2, Mode)
df_17Bcoord2$z_gamma2_q5 = factor(df_17Bcoord2$z_gamma2_q5, levels = c(1,2, 5, 3, 4))

df_17Bcoord2$z_gamma2_q5_alpha = pmax(colMeans(dfv2_A_gamma2_q5$z[-(1:500),]==1),
                              colMeans(dfv2_A_gamma2_q5$z[-(1:500),]==2),
                              colMeans(dfv2_A_gamma2_q5$z[-(1:500),]==3), 
                              colMeans(dfv2_A_gamma2_q5$z[-(1:500),]==4),
                              colMeans(dfv2_A_gamma2_q5$z[-(1:500),]==5))

pnodeconv2=ggplot(df_17Bcoord2, aes(x, y)) +
  geom_text(aes(color = factor(z_gamma2_q5), alpha = z_gamma2_q5_alpha),
            size = 10, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()
#deconv
df_17B2coord2 = df_17Bcoord2[rep(seq_len(nrow(df_17Bcoord2)), 7), ] #rbind 7 times
df_17B2coord2$j2 = 1:(7*nrow(df_17Bcoord2))

shift = rbind(expand.grid(c(1/3, -1/3), c(1/3,-1/3)), expand.grid(c(2/3, -2/3,0), 0))
coord_scale = c(77.65, 135.5)
shift = t(t(shift)*coord_scale)
shift_long = shift[rep(seq_len(7), each = nrow(df_17Bcoord2)), ]

df_17B2coord2$x = df_17B2coord2$x + shift_long[,"Var1"]
df_17B2coord2$y = df_17B2coord2$y + shift_long[,"Var2"]
ggplot(df_17B2coord2, aes(x, y)) +
  geom_text(aes(color = rep(logcounts(sce_17B)[grep("CHGA", rownames(sce_17B)),],7)), label = "\u2B22", size = 4, family = "Lucida Sans Unicode") +
  labs(color = "CHGA", x = NULL, y = NULL) +
  scale_color_viridis()+
  geom_text(data = df_17Bcoord2, aes(x = x, y = y+20), size = 11, label = "\u2B21", alpha = 0.3, family = "Lucida Sans Unicode")+
  theme_classic() + coord_fixed()

deconv2 = run_mcmc_hexdeconv(df = df_17Bcoord2, gamma = 2, q = 5, d = 5,nrep = 2000, prev = dfv2_A_gamma2_q5)
plot(deconv2$mu[,21], type = "l", ylab = "mu")
plot(sapply(deconv2$lambda,"[[",2), type = "l", ylab = "lambda")
plot(sapply(deconv2$Y,"[[",12), type = "l", ylab = "lambda")


df_17B2coord2$z_deconv2 = apply(deconv2$z[500:2000,], 2, Mode)
df_17B2coord2$z_deconv2_alpha = pmax(colMeans(deconv2$z[500:2000,]==1),
                             colMeans(deconv2$z[500:2000,]==2),
                             colMeans(deconv2$z[500:2000,]==3),
                             colMeans(deconv2$z[500:2000,]==4),
                             colMeans(deconv2$z[500:2000,]==5))
pdeconv2= ggplot(data = df_17B2coord2, aes(x = x, y = y))+
  geom_text(aes(color = factor(z_deconv2, levels = c(1,2,5,3,4)), alpha = z_deconv2_alpha), label = "\u2B23",
            size = 2.8, family = "Lucida Sans Unicode")+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+
  labs(color = "State", x = NULL, y = NULL)+theme_classic()+
  guides(shape = F, size = F, alpha = F)+
  geom_text(data = df_17Bcoord2, aes(x = x+2, y = y+20), size = 10, label = "\u2B21", alpha = 0.3, family = "Lucida Sans Unicode")+
  coord_fixed()
pnodeconv2+pdeconv2
df_17B2coord2$z_deconv3 = apply(deconv2$z[500:600,], 2, Mode)
df_17B2coord2$z_deconv3_alpha = pmax(colMeans(deconv2$z[500:600,]==1),
                                   colMeans(deconv2$z[500:600,]==2),
                                   colMeans(deconv2$z[500:600,]==3),
                                   colMeans(deconv2$z[500:600,]==4),
                                   colMeans(deconv2$z[500:600,]==5))
#deconv 3
#deconv3 = list(z = df_sim_z, mu = df_sim_mu, lambda = df_sim_lambda, Y = df_sim_Y, Ychange = testY)
deconv3 = run_mcmc_hexdeconv(df = df_17Bcoord2, gamma = 2, q = 5, d = 5,nrep = 2000, prev = dfv2_A_gamma2_q5)

plot(deconv3$mu[,21], type = "l", ylab = "mu")
plot(sapply(deconv3$Y,"[[",13), type = "l", ylab = "lambda")
plot(sapply(deconv3$lambda,"[[",13), type = "l", ylab = "lambda")
plot(deconv3$Ychange, type = "l")
df_17B2coord2$z_deconv3 = apply(deconv3$z[1000:2000,], 2, Mode)
df_17B2coord2$z_deconv3_alpha = pmax(colMeans(deconv3$z[1000:2000,]==1),
                                   colMeans(deconv3$z[1000:2000,]==2),
                                   colMeans(deconv3$z[1000:2000,]==3),
                                   colMeans(deconv3$z[1000:2000,]==4),
                                   colMeans(deconv3$z[1000:2000,]==5))
pdeconv3= ggplot(data = df_17B2coord2, aes(x = x, y = y))+
  geom_text(aes(color = factor(z_deconv3, levels = c(1,2,5,3,4)), alpha = z_deconv3_alpha), label = "\u2B23",
            size = 2.8, family = "Lucida Sans Unicode")+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+
  labs(color = "State", x = NULL, y = NULL)+theme_classic()+
  guides(shape = F, size = F, alpha = F)+
  geom_text(data = df_17Bcoord2, aes(x = x+2, y = y+20), size = 10, label = "\u2B21", alpha = 0.3, family = "Lucida Sans Unicode")+
  coord_fixed()

pnodeconv2+pdeconv3

simYmean = as.data.frame(Reduce(`+`, deconv3$Y[100:200])/101)
colnames(simYmean) = paste0("PC", 1:5)
pc1 = ggplot(data = df_17B2coord2, aes(x = x, y = y))+
  geom_text(aes(color = simYmean$PC1), label = "\u2B23",
            size = 2, family = "Lucida Sans Unicode")+
  scale_color_viridis(option = "A")+
  labs(color = "PC1", x = NULL, y = NULL)+theme_classic()+
  guides(shape = F, size = F, alpha = F)+
  #geom_text(data = df_17Bcoord2, aes(x = x+2, y = y+20), size = 11, label = "\u2B21", alpha = 0.3, family = "Lucida Sans Unicode")+
  coord_fixed()

pc2 = ggplot(data = df_17B2coord2, aes(x = x, y = y))+
  geom_text(aes(color = simYmean$PC2), label = "\u2B23",
            size = 2, family = "Lucida Sans Unicode")+
  scale_color_viridis(option = "A")+
  labs(color = "PC2", x = NULL, y = NULL)+theme_classic()+
  guides(shape = F, size = F, alpha = F)+
  #geom_text(data = df_17Bcoord2, aes(x = x+2, y = y+20), size = 11, label = "\u2B21", alpha = 0.3, family = "Lucida Sans Unicode")+
  coord_fixed()
pc3 = ggplot(data = df_17B2coord2, aes(x = x, y = y))+
  geom_text(aes(color = simYmean$PC3), label = "\u2B23",
            size = 2, family = "Lucida Sans Unicode")+
  scale_color_viridis(option = "A")+
  labs(color = "PC3", x = NULL, y = NULL)+theme_classic()+
  guides(shape = F, size = F, alpha = F)+
  #geom_text(data = df_17Bcoord2, aes(x = x+2, y = y+20), size = 11, label = "\u2B21", alpha = 0.3, family = "Lucida Sans Unicode")+
  coord_fixed()
pc4 = ggplot(data = df_17B2coord2, aes(x = x, y = y))+
  geom_text(aes(color = simYmean$PC4), label = "\u2B23",
            size = 2, family = "Lucida Sans Unicode")+
  scale_color_viridis(option = "A")+
  labs(color = "PC4", x = NULL, y = NULL)+theme_classic()+
  guides(shape = F, size = F, alpha = F)+
  #geom_text(data = df_17Bcoord2, aes(x = x+2, y = y+20), size = 11, label = "\u2B21", alpha = 0.3, family = "Lucida Sans Unicode")+
  coord_fixed()
(pc1+pc2)/(pc3+pc4)
pc1_0 = ggplot(df_17Bcoord2, aes(x, y)) +
  geom_text(aes(color = df_17Bcoord2$Y1),
            size = 5, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "PC1", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_color_viridis(option = "A")+
  theme_classic() + coord_fixed()
pc2_0 = ggplot(df_17Bcoord2, aes(x, y)) +
  geom_text(aes(color = df_17Bcoord2$Y2),
            size = 5, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "PC2", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_color_viridis(option = "A")+
  theme_classic() + coord_fixed()
pc3_0 = ggplot(df_17Bcoord2, aes(x, y)) +
  geom_text(aes(color = df_17Bcoord2$Y3),
            size = 5, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "PC3", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_color_viridis(option = "A")+
  theme_classic() + coord_fixed()
pc4_0 = ggplot(df_17Bcoord2, aes(x, y)) +
  geom_text(aes(color = df_17Bcoord2$Y4),
            size = 5, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "PC4", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_color_viridis(option = "A")+
  theme_classic() + coord_fixed()
(pc1_0+pc2_0)/(pc3_0+pc4_0)


Xhat = PCs$components[,1:5] %*% t(PCs$rotation[,1:5])
Xpca = prcomp(t(logcounts(sce_17B)))
Xhat2 =as.matrix(simYmean) %*% t(PCs$rotation[,1:5])
testdata = data.frame(x = Xhat[,"CD3E"], y = logcounts(sce_17B)[grep("CD3E$", rownames(sce_17B)),])
ggplot(testdata, aes(x = x, y = y+logcounts(sce_17B)[grep("CD3D", rownames(sce_17B)),]))+
  geom_point()+
  geom_smooth(se = F)+theme_classic()
reg_CD3 = lm(rep(logcounts(sce_17B)[grep("CD3D", rownames(sce_17B)),],7) ~ .,  data = simYmean)
#plot(predict(reg1), rep(logcounts(sce_17B)[grep("CD3D", rownames(sce_17B)),],7))
predCD3 = ggplot(data = df_17B2coord2, aes(x = x, y = y))+
  geom_text(aes(color = predict(reg_CD3)), label = "\u2B23",
            size = 2, family = "Lucida Sans Unicode")+
  scale_color_viridis(option = "A")+
  labs(color = "CD3", x = NULL, y = NULL)+theme_classic()+
  guides(shape = F, size = F, alpha = F)+
  coord_fixed()

reg_CHGA = lm(rep(logcounts(sce_17B)[grep("CHGA", rownames(sce_17B)),],7) ~ .,  data = simYmean)
predCHGA = ggplot(data = df_17B2coord2, aes(x = x, y = y))+
  geom_text(aes(color = predict(reg_CHGA)), label = "\u2B23",
            size = 2, family = "Lucida Sans Unicode")+
  scale_color_viridis(option = "A")+
  labs(color = "CHGA", x = NULL, y = NULL)+theme_classic()+
  guides(shape = F, size = F, alpha = F)+
  coord_fixed()
reg_CD8 = lm(rep(logcounts(sce_17B)[grep("CD8A", rownames(sce_17B)),],7) ~ .,  data = simYmean)
predCD8 = ggplot(data = df_17B2coord2, aes(x = x, y = y))+
  geom_text(aes(color = predict(reg_CD8)), label = "\u2B23",
            size = 2, family = "Lucida Sans Unicode")+
  scale_color_viridis(option = "A")+
  labs(color = "CD8", x = NULL, y = NULL)+theme_classic()+
  guides(shape = F, size = F, alpha = F)+
  coord_fixed()
predCHGA + predCD3 + predCD8
