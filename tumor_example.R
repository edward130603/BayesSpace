library(tidyverse)
library(DropletUtils)
library(scater)
library(ggbeeswarm)
library(scran)
library(viridis)
library(patchwork)

#Make sce
sce_A = read10xCounts("data-raw/044_A")
sce_A$Sample = "A"
rownames(sce_A) = uniquifyFeatureNames(rowData(sce_A)$ID, rowData(sce_A)$Symbol)
colnames(sce_A) = paste0(sce_A$Sample, '.', sce_A$Barcode)

#Add position
pos = read.csv("data-raw/044_A/179566_ID044_A_spatial__tissue_positions_list.csv", header=FALSE)
colnames(pos) = c("Barcode", "tissue", "Y1", "X1", "Y2", "X2")
colData(sce_A) = merge(colData(sce_A), pos, by = "Barcode")

#Preprocess
is.mito = grep("^MT", rownames(sce_A))
is.RBC = grep("^HB", rownames(sce_A))
df_qc = perCellQCMetrics(sce_A, subsets = list(Mito = is.mito,
                                          RBC = is.RBC))
qc = quickPerCellQC(df_qc, percent_subsets = "subsets_Mito_percent")
sce_A$discard = qc$discard
coord1=ggplot(as.data.frame(colData(sce_A)), aes(x = X2, y = Y2, col = clusters)) + 
  geom_text(label = "\u2B22", size = 6, family = "Lucida Sans Unicode")+
  theme_classic()+coord_fixed()+guides(col = F)
coord2=ggplot(as.data.frame(colData(sce_A)), aes(x = X1, y = Y1, col = sce_A$cluster)) + 
  geom_text(label = "\u2B25", size = 12, family = "Lucida Sans Unicode") +
  theme_classic() + coord_fixed()#+guides(col = F)
coord1+coord2
#Normalization
set.seed(100)
clusters = quickCluster(sce_A)
sce_A = computeSumFactors(sce_A, clusters = clusters)
sce_A = logNormCounts(sce_A)
plot(librarySizeFactors(sce_A), sizeFactors(sce_A), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", log="xy", 
     col = ifelse(sce_A$discard, "red", "black"))

set.seed(101)
dec <- modelGeneVarByPoisson(sce_A)
top <- getTopHVGs(dec, prop=0.1)
plot(dec$mean, dec$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)

#Dim reduction and clustering
set.seed(102)
sce_A <- denoisePCA(sce_A, technical=dec, subset.row=top)
sce_A <- runTSNE(sce_A, dimred="PCA")

snn.gr <- buildSNNGraph(sce_A, use.dimred="PCA", k=40)
sce_A$cluster <- factor(igraph::cluster_walktrap(snn.gr)$membership)
table(sce_A$cluster)
sce_A$cluster = factor(sce_A$cluster, levels = c(1,2,4,5,3))




tsne1 = ggplot(as.data.frame(reducedDim(sce_A, "TSNE")), aes(x = V1, y = V2, col = sce_A$cluster))+
  geom_point(size = 2) +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()
  
spatial1 = ggplot(as.data.frame(colData(sce_A)), aes(x = X1, y = Y1, col = cluster)) + 
  geom_text(label = "\u2B22", 
            size = 7, family = "Lucida Sans Unicode")+
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()+
  theme(legend.position = c(0.05,0.8)) + coord_fixed()

tsne1 + spatial1 +
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "A")
  

#Distribution of principal components
pc1_spatial = ggplot(as.data.frame(colData(sce_A)), aes(x = X1, y = Y1, col = reducedDim(sce_A, "PCA")[,"PC1"])) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "PC1")+scale_color_viridis()+theme_classic()
pc2_spatial = ggplot(as.data.frame(colData(sce_A)), aes(x = X1, y = Y1, col = reducedDim(sce_A, "PCA")[,"PC2"])) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "PC2")+scale_color_viridis()+theme_classic()
pc3_spatial = ggplot(as.data.frame(colData(sce_A)), aes(x = X1, y = Y1, col = reducedDim(sce_A, "PCA")[,"PC3"])) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "PC3")+scale_color_viridis()+theme_classic()
pc4_spatial = ggplot(as.data.frame(colData(sce_A)), aes(x = X1, y = Y1, col = reducedDim(sce_A, "PCA")[,"PC4"])) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "PC4")+scale_color_viridis()+theme_classic()
pc5_spatial = ggplot(as.data.frame(colData(sce_A)), aes(x = X1, y = Y1, col = reducedDim(sce_A, "PCA")[,"PC5"])) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "PC5")+scale_color_viridis()+theme_classic()
pc1_spatial + pc2_spatial 
pc3_spatial + pc4_spatial + pc5_spatial

#Run mcmc
df_A = data.frame(reducedDim(sce_A, "PCA"))[,1:5]
colnames(df_A) = paste0("Y", 1:5)
df_A$x = sce_A$X1
df_A$y = sce_A$Y1
df_A$j = 1:nrow(df_A)
num_neighbors = sapply(sapply(1:nrow(df_A), function(x){df_A[(abs(df_A[,"x"] -df_A[x,"x"]) + abs(df_A[,"y"] - df_A[x,"y"])) == 2,"j"]}), length)
df_A = df_A[num_neighbors >0, ] #all have at least 1 neighbor (in fact at least 3)

df_A_gamma2_q5 = run_mcmc_multi(df = df_A, gamma = 4, q = 5, d = 5,nrep = 2000)
df_A_gamma2_q5v2 = df_A_gamma2_q5
saveRDS(df_A_gamma2_q5, "mcmc_d5_q5_A.RDS")
df_A_gamma2_q5 = readRDS("data/mcmc_d5_q5_A.RDS")

#5 states, gamma = 2
df_A$z_gamma2_q5 = apply(df_A_gamma2_q5$z[-(1:500),], 2, Mode)
df_A$z_gamma2_q5[df_A$z_gamma2_q5 == 5] = 99
df_A$z_gamma2_q5[df_A$z_gamma2_q5 == 3] = 5
df_A$z_gamma2_q5[df_A$z_gamma2_q5 == 99] = 3

df_A$z_gamma2_q5_alpha = pmax(colMeans(df_A_gamma2_q5$z[-(1:500),]==1),
                               colMeans(df_A_gamma2_q5$z[-(1:500),]==2),
                               colMeans(df_A_gamma2_q5$z[-(1:500),]==3),
                               colMeans(df_A_gamma2_q5$z[-(1:500),]==4),
                               colMeans(df_A_gamma2_q5$z[-(1:500),]==5))

ggplot(df_A, aes(x, y)) +
  geom_text(aes(color = factor(z_gamma2_q5), alpha = z_gamma2_q5_alpha),
            size = 13, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()

ggplot(as.data.frame(colData(sce_A)), aes(x = X1, y = Y1, col = sce_A$cluster)) + 
  geom_text(label = "\u2B22", size = 13, family = "Lucida Sans Unicode") +
  theme_classic() + coord_fixed() + labs(col = "State")

ggplot(df_A, aes(x, y)) +
  geom_point(aes(color = factor(z_gamma2_q5), alpha = z_gamma2_q5_alpha), size = 4, shape = 18) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
  theme_classic()
ggplot(as.data.frame(reducedDim(sce_A, "TSNE")[num_neighbors>0,]), aes(x = V1, y = V2, col = factor(df_A$z_gamma2_q5)))+
  geom_point() +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()


ggplot(df_A, aes(x, y, color = logcounts(sce_A)[grep("CD3D", rownames(sce_A)),] +
                   logcounts(sce_A)[grep("CD3E$", rownames(sce_A)),]+
                   logcounts(sce_A)[grep("CD3G", rownames(sce_A)),] > 0.1)) +
  geom_text(label = "\u2B25", size = 21, family = "Lucida Sans Unicode") +
  labs(color = "CD3", x = NULL, y = NULL) +
  scale_color_viridis_d()+
  theme_classic() + coord_fixed()

ggplot(df_A, aes(x, y, color = logcounts(sce_A)[grep("CHGA", rownames(sce_A)),])) +
  geom_text(label = "\u2B25", size = 21, family = "Lucida Sans Unicode") +
  labs(color = "CHGA", x = NULL, y = NULL) +
  scale_color_viridis()+
  theme_classic() + coord_fixed()

ggplot(df_A, aes(x, y, color = logcounts(sce_A)[grep("KRT20", rownames(sce_A)),])) +
  geom_text(label = "\u2B25", size = 21, family = "Lucida Sans Unicode") +
  labs(color = "CHGA", x = NULL, y = NULL) +
  scale_color_viridis()+
  theme_classic() + coord_fixed()

ggplot(as.data.frame(reducedDim(sce_A, "TSNE")[num_neighbors>0,]), aes(x = V1, y = V2, col = logcounts(sce_A)[grep("KRT20", rownames(sce_A)),]))+
  geom_point() +
  scale_color_viridis() +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()

##Deconvolution
deconv1 = run_mcmc_hexdeconv(df = df_A, gamma = 2, q = 5, d = 5,nrep = 2000, prev = df_A_gamma2_q5)
save(deconv1, file = "data/deconv1.RDS")
df_A2 = df_A[rep(seq_len(nrow(df_A)), 7), ] #rbind 7 times
df_A2$j2 = 1:(7*nrow(df_A))

shift = rbind(expand.grid(c(1/2, -1/2), c(1/4,-1/4)), expand.grid(0, c(0, 1/2, -1/2)))
shift_long = shift[rep(seq_len(7), each = nrow(df_A)), ]

df_A2$x = df_A2$x + shift_long$Var1
df_A2$y = df_A2$y + shift_long$Var2

df_A2$z_deconv1 = apply(deconv1$z[500:2000,], 2, Mode)
df_A2$z_deconv1_alpha = pmax(colMeans(deconv1$z[500:2000,]==1),
                           colMeans(deconv1$z[500:2000,]==2),
                           colMeans(deconv1$z[500:2000,]==3),
                           colMeans(deconv1$z[500:2000,]==4),
                           colMeans(deconv1$z[500:2000,]==5))


plot(deconv1$mu[,1], type = "l", ylab = "mu")
plot(deconv1$mu[,6], type = "l", ylab = "mu")
plot(sapply(deconv1$lambda,"[[",1), type = "l", ylab = "lambda")
plot(sapply(deconv1$lambda,"[[",2), type = "l", ylab = "lambda")


ggplot(data = df_A2, aes(x = x, y = y))+
  geom_text(aes(color = factor(z_deconv1), alpha = z_deconv1_alpha), label = "\u2B23",
            size = 4, family = "Lucida Sans Unicode")+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+
  labs(color = "State")+theme_classic()+
  guides(shape = F, size = F, alpha = F)+
  geom_text(data = df_A, aes(x = x+0.03, y = y+0.2), size = 13, label = "\u2B21", alpha = 0.3, family = "Lucida Sans Unicode")+
  coord_fixed()

##New coordinates
df_Acoord2 = df_A  
df_Acoord2$x = sce_A$X2
df_Acoord2$y = sce_A$Y2
# index = 2
# test = rep(1, n)
# test[index] = 2
# test[df_j[[index]]] = 3
# ggplot(df_Acoord2, aes(x, y, color = factor(test))) +
#   geom_text(label = "\u2B22", size = 10, family = "Lucida Sans Unicode") +
#   labs(x = NULL, y = NULL) +
#   scale_color_viridis_d()+
#   theme_classic() + coord_fixed()
dfv2_A_gamma2_q5 = run_mcmc_multi(df = df_Acoord2, gamma = 2, q = 5, d = 5,nrep = 2000)
df_Acoord2$z_gamma2_q5 = apply(dfv2_A_gamma2_q5$z[-(1:500),], 2, Mode)
df_Acoord2$z_gamma2_q5 = factor(df_Acoord2$z_gamma2_q5, levels = c(1,2, 5, 3, 4))

df_Acoord2$z_gamma2_q5_alpha = pmax(colMeans(dfv2_A_gamma2_q5$z[-(1:500),]==1),
                              colMeans(dfv2_A_gamma2_q5$z[-(1:500),]==2),
                              colMeans(dfv2_A_gamma2_q5$z[-(1:500),]==3),
                              colMeans(dfv2_A_gamma2_q5$z[-(1:500),]==4),
                              colMeans(dfv2_A_gamma2_q5$z[-(1:500),]==5))

ggplot(df_Acoord2, aes(x, y)) +
  geom_text(aes(color = factor(z_gamma2_q5), alpha = z_gamma2_q5_alpha),
            size = 10, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
  theme_classic() + coord_fixed()
