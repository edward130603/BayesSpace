library(tidyverse)
library(DropletUtils)
library(scater)
library(ggbeeswarm)
library(scran)
library(viridis)
library(patchwork)
library(SingleR)


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
ggplot(as.data.frame(colData(sce_A)), aes(x = X1, y = Y1, col = discard)) + 
  geom_point(size = 4) + theme_classic()

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

tsne1 = ggplot(as.data.frame(reducedDim(sce_A, "TSNE")), aes(x = V1, y = V2, col = sce_A$cluster))+
  geom_point(size = 2) +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()
  
spatial1 = ggplot(as.data.frame(colData(sce_A)), aes(x = X1, y = Y1, col = cluster)) + 
  geom_point(size = 2) +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()

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

df_A_gamma2_q5 = run_mcmc_multi(df = df_A, gamma = 2, q = 5, d = 5,nrep = 1000)
saveRDS(df_A_gamma2_q5, "mcmc_d5_q5_A.RDS")
df_A2_gamma6_q3 = run_mcmc_multi(df = df_A2, gamma = 6, q = 3, d = 5,nrep = 1000)
df_A2_gamma2_q5 = run_mcmc_multi(df = df_A2, gamma = 2, q = 5, d = 5,nrep = 1000)
df_A2_gamma4_q5 = run_mcmc_multi(df = df_A2, gamma = 4, q = 5, d = 5,nrep = 1000)

#5 states, gamma = 2
df_A$z_gamma2_q5 = apply(df_A_gamma2_q5$z[-(1:500),], 2, Mode)
df_A$z_gamma2_q5_alpha = pmax(colMeans(df_A_gamma2_q5$z[-(1:500),]==1),
                               colMeans(df_A_gamma2_q5$z[-(1:500),]==2),
                               colMeans(df_A_gamma2_q5$z[-(1:500),]==3),
                               colMeans(df_A_gamma2_q5$z[-(1:500),]==4),
                               colMeans(df_A_gamma2_q5$z[-(1:500),]==5))
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


ggplot(df_A, aes(x, y)) +
  geom_point(aes(color = logcounts(sce_A)[grep("NCAM1$", rownames(sce_A)),]), size = 2) +
  labs(color = "CD62L", x = NULL, y = NULL) +
  scale_color_viridis()+
  theme_classic()

ggplot(as.data.frame(reducedDim(sce_A, "TSNE")[num_neighbors>0,]), aes(x = V1, y = V2, col = logcounts(sce_A)[grep("KRT20", rownames(sce_A)),]))+
  geom_point() +
  scale_color_viridis() +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()

library(SingleR)
BPE = BlueprintEncodeData()
pred = SingleR(test = sce_A, ref = BPE, labels = BPE$label.main)
pred$pruned.labels[is.na(pred$pruned.labels)] = "Unknown"
ggplot(df_A, aes(x, y)) +
  geom_point(aes(color = pred$pruned.labels), size = 4) +
  labs(color = "State", x = NULL, y = NULL) +
  theme_classic()

pred.cluster = SingleR(test = sce_A, ref = BPE, labels = BPE$label.main, method = "cluster", clusters = df_A$z_gamma2_q5)
pred.cluster$pruned.labels[is.na(pred.cluster$pruned.labels)] = "Unknown"
ggplot(df_A, aes(x, y)) +
  geom_point(aes(color = pred$pruned.labels), size = 4) +
  labs(color = "State", x = NULL, y = NULL) +
  theme_classic()


##Deconvolution
deconv1 = run_mcmc_deconv(df = df_A, gamma = 2, q = 5, d = 5,nrep = 1000, prev = df_A_gamma2_q5)
save(deconv1, file = "data/deconv1.RDS")

df2$z_deconv1 = apply(deconv1$z[900:1000,], 2, Mode)
df2$z_deconv1_alpha = pmax(colMeans(deconv1$z[900:1000,]==1),
                           colMeans(deconv1$z[900:1000,]==2),
                           colMeans(deconv1$z[900:1000,]==3),
                           colMeans(deconv1$z[900:1000,]==4),
                           colMeans(deconv1$z[900:1000,]==5))

# cairo_pdf("output/my_plot2.pdf", family = "Lucida Sans Unicode", 10, 5)
# ggplot(data = df2, aes(x = x, y = y))+
#   geom_point(aes(color = factor(z_deconv1), shape = factor(rep(1:4, each = 1500)),
#                  alpha = z_deconv1_alpha), size = 2)+
#   scale_shape_manual(values = c("\u25E3", "\u25E2", "\u25E4", "\u25E5")) +
#   scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
#   labs(color = "State")+theme_classic()+
#   guides(shape = F, size = F, alpha = F)
# dev.off()

ggplot(data = df2, aes(x = x, y = y))+
  geom_text(aes(color = factor(z_deconv1), label = rep(c("\u25E3", "\u25E2", "\u25E4", "\u25E5"), each = 1500),
                 alpha = z_deconv1_alpha), size = 4, family = "Arial Unicode MS")+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.2,1,0.1), range = c(0,1))+
  labs(color = "State")+theme_classic()+
  guides(shape = F, size = F, alpha = F)


plot(deconv1$mu[,1], type = "l", ylab = "mu")
plot(deconv1$mu[,6], type = "l", ylab = "mu")
plot(sapply(deconv1$lambda,"[[",1), type = "l", ylab = "lambda")
plot(sapply(deconv1$lambda,"[[",2), type = "l", ylab = "lambda")
