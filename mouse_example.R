library(tidyverse)
library(DropletUtils)
library(scater)
library(ggbeeswarm)
library(AnnotationHub)
library(scran)
library(viridis)
library(patchwork)


#Make sce
sce_A2 = read10xCounts("data-raw/A2")
sce_A2$Sample = "A2"
rownames(sce_A2) = uniquifyFeatureNames(rowData(sce_A2)$ID, rowData(sce_A2)$Symbol)
colnames(sce_A2) = paste0(sce_A2$Sample, '.', sce_A2$Barcode)

#Add position
pos = read.csv("data-raw/A2/tissue_positions_list.txt", header=FALSE)
colnames(pos) = c("Barcode", "tissue", "Y1", "X1", "Y2", "X2")
colData(sce_A2) = merge(colData(sce_A2), pos, by = "Barcode")

#Preprocess
ens.mm.v97 = AnnotationHub()[["AH73905"]]
rowData(sce_A2)$SEQNAME = mapIds(ens.mm.v97, keys=rowData(sce_A2)$ID,
                                   keytype="GENEID", column="SEQNAME")

is.mito = which(rowData(sce_A2)$SEQNAME == "MT")
is.RBC = grep("^Hb", rownames(sce_A2))
df_qc = perCellQCMetrics(sce_A2, subsets = list(Mito = is.mito,
                                          RBC = is.RBC))
qc = quickPerCellQC(df_qc, percent_subsets = "subsets_Mito_percent")
sce_A2$discard = qc$discard
ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = discard)) + 
  geom_point(size = 4)

#Normalization
set.seed(100)
clusters = quickCluster(sce_A2)
sce_A2 = computeSumFactors(sce_A2, clusters = clusters)
sce_A2 = logNormCounts(sce_A2)
plot(librarySizeFactors(sce_A2), sizeFactors(sce_A2), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", log="xy", 
     col = ifelse(sce_A2$discard, "red", "black"))

set.seed(101)
dec <- modelGeneVarByPoisson(sce_A2)
top <- getTopHVGs(dec, prop=0.1)
plot(dec$mean, dec$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)

#Dim reduction and clustering
set.seed(102)
sce_A2 <- denoisePCA(sce_A2, technical=dec, subset.row=top)
sce_A2 <- runTSNE(sce_A2, dimred="PCA")

snn.gr <- buildSNNGraph(sce_A2, use.dimred="PCA", k=100)
sce_A2$cluster <- factor(igraph::cluster_walktrap(snn.gr)$membership)
table(sce_A2$cluster)

tsne1 = ggplot(as.data.frame(reducedDim(sce_A2, "TSNE")), aes(x = V1, y = V2, col = sce_A2$cluster))+
  geom_point() +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()
  
spatial1 = ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = cluster)) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()

tsne1 + spatial1 +
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "A")
  

#Distribution of principal components
pc1 = ggplot(as.data.frame(reducedDim(sce_A2, "PCA")), aes(x = PC1))+geom_histogram(bins = 20)+
  facet_wrap(~sce_A2$cluster, scales = "free") +
  labs(y = NULL) +
  theme_classic()
pc2 = ggplot(as.data.frame(reducedDim(sce_A2, "PCA")), aes(x = PC2))+geom_histogram(bins = 20)+
  facet_wrap(~sce_A2$cluster, scales = "free") +
  labs(y = NULL) +
  theme_classic()
pc1 / pc2 + plot_annotation(tag_levels = "A")
pc1_spatial = ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = reducedDim(sce_A2, "PCA")[,"PC1"])) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "PC1")+scale_color_viridis()+theme_classic()
pc2_spatial = ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = reducedDim(sce_A2, "PCA")[,"PC2"])) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "PC2")+scale_color_viridis()+theme_classic()
pc3_spatial = ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = reducedDim(sce_A2, "PCA")[,"PC3"])) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "PC3")+scale_color_viridis()+theme_classic()
pc4_spatial = ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = reducedDim(sce_A2, "PCA")[,"PC4"])) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "PC4")+scale_color_viridis()+theme_classic()
pc5_spatial = ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = reducedDim(sce_A2, "PCA")[,"PC5"])) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "PC5")+scale_color_viridis()+theme_classic()
pc1_spatial + pc2_spatial 
pc3_spatial + pc4_spatial + pc5_spatial

#Run mcmc
df_A2 = data.frame(reducedDim(sce_A2, "PCA"))
colnames(df_A2) = paste0("Y", 1:5)
df_A2$x = sce_A2$X1
df_A2$y = sce_A2$Y1
df_A2$j = 1:nrow(df_A2)
num_neighbors = sapply(sapply(1:nrow(df_A2), function(x){df_A2[(abs(df_A2[,"x"] -df_A2[x,"x"]) + abs(df_A2[,"y"] - df_A2[x,"y"])) == 2,"j"]}), length)
df_A2 = df_A2[num_neighbors >0, ]

df_A2_gamma2_q3 = run_mcmc_multi(df = df_A2, gamma = 2, q = 3, d = 5,nrep = 1000)
df_A2_gamma6_q3 = run_mcmc_multi(df = df_A2, gamma = 6, q = 3, d = 5,nrep = 1000)
df_A2_gamma2_q5 = run_mcmc_multi(df = df_A2, gamma = 2, q = 5, d = 5,nrep = 1000)
df_A2_gamma4_q5 = run_mcmc_multi(df = df_A2, gamma = 4, q = 5, d = 5,nrep = 1000)

#3 states, gamma = 2
df_A2$z_gamma2_q3 = apply(df_A2_gamma2_q3$z[-(1:100),], 2, Mode)
df_A2$z_gamma2_q3_alpha = pmax(colMeans(df_A2_gamma2_q3$z[-(1:100),]==1),
                         colMeans(df_A2_gamma2_q3$z[-(1:100),]==2),
                         colMeans(df_A2_gamma2_q3$z[-(1:100),]==3))
ggplot(df_A2, aes(x, y)) +
  geom_point(aes(color = factor(z_gamma2_q3), alpha = z_gamma2_q3_alpha), size = 4) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+
  theme_classic()
ggplot(as.data.frame(reducedDim(sce_A2, "TSNE")[num_neighbors>0,]), aes(x = V1, y = V2, col = factor(df_A2$z_gamma2_q3)))+
  geom_point() +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()

#3 states, gamma = 6
df_A2$z_gamma6_q3 = apply(df_A2_gamma6_q3$z[-(1:600),], 2, Mode)
df_A2$z_gamma6_q3_alpha = pmax(colMeans(df_A2_gamma6_q3$z[-(1:600),]==1),
                               colMeans(df_A2_gamma6_q3$z[-(1:600),]==2),
                               colMeans(df_A2_gamma6_q3$z[-(1:600),]==3))
ggplot(df_A2, aes(x, y)) +
  geom_point(aes(color = factor(z_gamma6_q3), alpha = z_gamma6_q3_alpha), size = 4) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+
  theme_classic()
ggplot(as.data.frame(reducedDim(sce_A2, "TSNE")[num_neighbors>0,]), aes(x = V1, y = V2, col = factor(df_A2$z_gamma6_q3)))+
  geom_point() +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()

#5 states, gamma = 2
df_A2$z_gamma2_q5 = apply(df_A2_gamma2_q5$z[-(1:500),], 2, Mode)
df_A2$z_gamma2_q5_alpha = pmax(colMeans(df_A2_gamma2_q5$z[-(1:500),]==1),
                               colMeans(df_A2_gamma2_q5$z[-(1:500),]==2),
                               colMeans(df_A2_gamma2_q5$z[-(1:500),]==3),
                               colMeans(df_A2_gamma2_q5$z[-(1:500),]==4),
                               colMeans(df_A2_gamma2_q5$z[-(1:500),]==5))
ggplot(df_A2, aes(x, y)) +
  geom_point(aes(color = factor(z_gamma2_q5), alpha = z_gamma2_q5_alpha), size = 4) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
  theme_classic()
ggplot(as.data.frame(reducedDim(sce_A2, "TSNE")[num_neighbors>0,]), aes(x = V1, y = V2, col = factor(df_A2$z_gamma2_q5)))+
  geom_point() +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()

#5 states, gamma = 4
df_A2$z_gamma4_q5 = apply(df_A2_gamma4_q5$z[-(1:500),], 2, Mode)
df_A2$z_gamma4_q5_alpha = pmax(colMeans(df_A2_gamma4_q5$z[-(1:500),]==1),
                               colMeans(df_A2_gamma4_q5$z[-(1:500),]==2),
                               colMeans(df_A2_gamma4_q5$z[-(1:500),]==3),
                               colMeans(df_A2_gamma4_q5$z[-(1:500),]==4),
                               colMeans(df_A2_gamma4_q5$z[-(1:500),]==5))
ggplot(df_A2, aes(x, y)) +
  geom_point(aes(color = factor(z_gamma4_q5), alpha = z_gamma4_q5_alpha), size = 4) +
  labs(color = "State", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
  theme_classic()
ggplot(as.data.frame(reducedDim(sce_A2, "TSNE")[num_neighbors>0,]), aes(x = V1, y = V2, col = factor(df_A2$z_gamma4_q5)))+
  geom_point() +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()

plot(df_A2_gamma4_q5$mu[,5], type = "l", ylab = "mu5")
plot(sapply(df_A2_gamma4_q5$lambda,"[[",7), type = "l", ylab = "lambda2.2")

plot(df_A2_gamma6_q3$mu[,2], type = "l", ylab = "mu5")

