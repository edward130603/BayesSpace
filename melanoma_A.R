library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(scran)
library(viridis)
library(patchwork)

data1 = read_delim("data-raw/ST-Melanoma-Datasets_1/ST_mel1_rep2_counts.tsv", delim = "\t")
#data2 = read.delim("data-raw/PDAC/GSM3036911.tsv", stringsAsFactors=FALSE)
rowdata1 = data1$gene
positions = colnames(data1)[-1]
data1 = as.matrix(data1[,-1])
positions = data.frame(positions) %>%
  separate(positions, c("x","y"), "x") %>%
  mutate_all(as.numeric)
rownames(data1) = rowdata1


ggplot(positions, aes(y, x, col = data1["ANXA2 ENSG00000182718",]))+
  geom_point(size = 5)+coord_fixed()

#Make sce
sce_A = SingleCellExperiment(assays = list(counts = data1), colData = positions)

#Preprocess
df_qc = perCellQCMetrics(sce_A)
qc = quickPerCellQC(df_qc)
ggplot(positions, aes(y, x, col = qc$discard))+
  geom_point(size = 5)+coord_fixed()
sce_A$discard = qc$discard

#Normalization
set.seed(100)
clusters = quickCluster(sce_A)
sce_A = computeSumFactors(sce_A, clusters = clusters)
sce_A = logNormCounts(sce_A)
plot(librarySizeFactors(sce_A), sizeFactors(sce_A), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", log="xy", 
     col = clusters)

set.seed(101)
dec <- modelGeneVarByPoisson(sce_A)
top <- getTopHVGs(dec, prop=0.1)
plot(dec$mean, dec$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression",
     col = ifelse(sce_A$discard, "red", "black"))
curfit <- metadata(dec)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)

#Dim reduction and clustering
set.seed(102)
sce_A <- denoisePCA(sce_A, technical=dec, subset.row=top)
PCs = getDenoisedPCs(sce_A, technical=dec, subset.row=top, min.rank = 5, max.rank = 50)
sce_A <- runTSNE(sce_A, dimred="PCA")
sce_A <- runUMAP(sce_A, dimred="PCA")

snn.gr <- buildSNNGraph(sce_A, use.dimred="PCA", k=10)
sce_A$cluster <- factor(igraph::cluster_walktrap(snn.gr)$membership)
table(sce_A$cluster)

ggplot(as.data.frame(colData(sce_A)), aes(x = y, y = x, fill = logcounts(sce_A)[grep("CD8A", rownames(sce_A)),])) + 
  geom_point(pch = 22, size = 7)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed()
ggplot(as.data.frame(colData(sce_A)), aes(x = y, y = x, fill = reducedDim(sce_A, "UMAP")[,2])) + 
  geom_point(pch = 22, size = 7)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed()
ggplot(as.data.frame(colData(sce_A)), aes(x = y, y = x, fill = sizeFactors(sce_A))) + 
  geom_point(pch = 22, size = 7)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed()

tsne1 + spatial1 +
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "A")

#Run mcmc
km6 = kmeans(reducedDim(sce_A, "PCA"), centers = 6)$cluster
km4 = kmeans(reducedDim(sce_A, "PCA"), centers = 4)$cluster
km3 = kmeans(reducedDim(sce_A, "PCA"), centers = 3)$cluster
km7 = kmeans(reducedDim(sce_A, "PCA"), centers = 7)$cluster
km8 = kmeans(reducedDim(sce_A, "PCA"), centers = 8)$cluster
km5 = kmeans(reducedDim(sce_A, "PCA"), centers = 5)$cluster
km9 = kmeans(reducedDim(sce_A, "PCA"), centers = 9)$cluster
km10 = kmeans(reducedDim(sce_A, "PCA"), centers = 10)$cluster
clust4 = cluster(Y= PCs$components[,1:9], positions = as.matrix(positions), q = 4, init = km4, nrep = 10000, gamma = 1.5, dist = 1)
clust3 = cluster(Y= PCs$components[,1:9], positions = as.matrix(positions), q = 3, init = km3, nrep = 10000, gamma = 1.5, dist = 1)
clust4_g3 = cluster(Y= PCs$components[,1:10], positions = as.matrix(positions), q = 4, init = km4, nrep = 10000, gamma = 3, dist = 1)
clust5 = cluster(Y= PCs$components[,1:10], positions = as.matrix(positions), q = 5, init = km5, nrep = 10000, gamma = 2, dist = 1)
clust6 = cluster(Y= PCs$components[,1:10], positions = as.matrix(positions), q = 6, init = km6, nrep = 10000, gamma = 2, dist = 1)
clust8 = cluster(Y= PCs$components[,1:10], positions = as.matrix(positions), q = 8, init = km8, nrep = 10000, gamma = 3, dist = 1)
clust9 = cluster(Y= PCs$components[,1:10], positions = as.matrix(positions), q = 9, init = km9, nrep = 10000, gamma = 3, dist = 1)

plot(clust4$plogLik, type = "l")
clust4col = apply(clust4$z[9000:10000,], 2, Mode)
clust3col = apply(clust3$z[9000:10000,], 2, Mode)
clust3_alpha = pmax(colMeans(clust3$z[9000:10000,]==1),
                    colMeans(clust3$z[9000:10000,]==2),
                    colMeans(clust3$z[9000:10000,]==3))
clust4_alpha = pmax(colMeans(clust4$z[9000:10000,]==1),
                     colMeans(clust4$z[9000:10000,]==2),
                     colMeans(clust4$z[9000:10000,]==3),
                     colMeans(clust4$z[9000:10000,]==4))
clust4g3col = apply(clust4_g3$z[9000:10000,], 2, Mode)
clust5col = apply(clust5$z[9000:10000,], 2, Mode)
clust6col = apply(clust6$z[9000:10000,], 2, Mode)
clust8col = apply(clust8$z[9000:10000,], 2, Mode)
clust9col = apply(clust9$z[9000:10000,], 2, Mode)
p00 = ggplot(as.data.frame(colData(sce_A)), aes(x = y, y = x, fill = factor(km3))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, color = "Cluster") +
  #scale_fill_viridis_d(option = "A")+
  scale_fill_manual(values = c("purple", "red", "yellow", "blue"))+
  guides(fill = F)+
  theme_void()+ coord_fixed()


ggplot(as.data.frame(colData(sce_A)), aes(x = x, y = y, fill = factor(clust9col))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, color = "Cluster") +
  scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()
p0 = ggplot(as.data.frame(colData(sce_A)), aes(x = y, y = x, fill = factor(clust3col), alpha = clust3_alpha)) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  guides(alpha = F, fill = F)+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+ 
  scale_fill_manual(values = c("purple", "red", "yellow", "blue"))+
  #scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()

#Deconvolution
df = data.frame(PCs$components[,1:9])
colnames(df) = paste0("Y",1:9)
df$x = positions$x
df$y = positions$y
df$j = 1:nrow(df)
saveRDS(clust3, "data-raw/clust3_melanomaA.RDS")
saveRDS(df, "data-raw/df_melanomaA.RDS")
deconv1 = deconvolve(Y= PCs$components[,1:9], positions = positions, q = 3, init = rep(clust3col, 9), nrep = 10000, gamma = 1, dist = 0.34)
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
