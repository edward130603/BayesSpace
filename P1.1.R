library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(scran)
library(viridis)
library(patchwork)

data1 = read.delim("data-raw/prostate-twelve-2/prostate-twelve/P1.2.tsv", stringsAsFactors=FALSE)
data1 = as.matrix(data1)
positions = as.tibble((str_split(colnames(data1), "[xX]", simplify = T)))[,-1]
colnames(positions) = c("x", "y")
positions = mutate_all(positions, as.numeric)

ggplot(positions, aes(x, -y, col = data1["ACTB ENSG00000075624",]))+
  geom_point(size = 5)+coord_fixed()

#Make sce
sce_1.2 = SingleCellExperiment(assays = list(counts = data1), colData = positions)

#Preprocess
df_qc = perCellQCMetrics(sce_1.2)
qc = quickPerCellQC(df_qc)
ggplot(positions, aes(x, -y, col = qc$discard))+
  geom_point(size = 5)+coord_fixed()
sce_1.2$discard = qc$discard

#Normalization
set.seed(100)
clusters = quickCluster(sce_1.2)
sce_1.2 = computeSumFactors(sce_1.2, clusters = clusters)
sce_1.2 = logNormCounts(sce_1.2)
plot(librarySizeFactors(sce_1.2), sizeFactors(sce_1.2), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", log="xy", 
     col = ifelse(clust7col==1, "red", "black"))

set.seed(101)
dec <- modelGeneVarByPoisson(sce_1.2)
top <- getTopHVGs(dec, prop=0.1)
plot(dec$mean, dec$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression",
     col = ifelse(sce_1.2$discard, "red", "black"))
curfit <- metadata(dec)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)

#Dim reduction and clustering
set.seed(102)
sce_1.2 <- denoisePCA(sce_1.2, technical=dec, subset.row=top)
PCs = getDenoisedPCs(sce_1.2, technical=dec, subset.row=top, min.rank = 5, max.rank = 50)
sce_1.2 <- runTSNE(sce_1.2, dimred="PCA")
sce_1.2 <- runUMAP(sce_1.2, dimred="PCA")

snn.gr <- buildSNNGraph(sce_1.2, use.dimred="PCA", k=40)
sce_1.2$cluster <- factor(igraph::cluster_walktrap(snn.gr)$membership)
table(sce_1.2$cluster)

tsne1 = ggplot(as.data.frame(reducedDim(sce_1.2, "TSNE")), aes(x = V1, y = V2, col = PCs$components[,1]))+
  geom_point(size = 2) +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()
logcounts(sce_1.2)[grep("CHGA", rownames(sce_1.2)),num_neighbors>0]
ggplot(as.data.frame(colData(sce_1.2)), aes(x = x, y = -y, fill = logcounts(sce_1.2)[grep("TGFB2", rownames(sce_1.2)),])) + 
  geom_point(pch = 22, size = 7)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed()
ggplot(as.data.frame(colData(sce_1.2)), aes(x = x, y = -y, fill = reducedDim(sce_1.2, "UMAP")[,2])) + 
  geom_point(pch = 22, size = 7)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed()

tsne1 + spatial1 +
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "A")
  
#Run mcmc
km6 = kmeans(reducedDim(sce_1.2, "PCA"), centers = 6)$cluster
km4 = kmeans(reducedDim(sce_1.2, "PCA"), centers = 4)$cluster
km7 = kmeans(reducedDim(sce_1.2, "PCA"), centers = 7)$cluster
km8 = kmeans(reducedDim(sce_1.2, "PCA"), centers = 8)$cluster
km5 = kmeans(reducedDim(sce_1.2, "PCA"), centers = 5)$cluster
km9 = kmeans(reducedDim(sce_1.2, "PCA"), centers = 9)$cluster
km10 = kmeans(reducedDim(sce_1.2, "PCA"), centers = 10)$cluster
clust7 = cluster(Y= PCs$components[,1:15], positions = as.matrix(positions), q = 7, init = km7, nrep = 100000, gamma = 2)
plot(clust7$plogLik, type = "l")
clust7col = apply(clust7$z[90000:100000,], 2, Mode)
ggplot(as.data.frame(colData(sce_1.2)), aes(x = x, y = -y, fill = factor(km7))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, color = "Cluster") +
  scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()
 ggplot(as.data.frame(colData(sce_1.2)), aes(x = x, y = -y, fill = factor(clust7col))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()

ggplot(as.data.frame(colData(sce_1.2)), aes(x = x, y = -y, fill = PCs$components[,3])) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = NULL) +
  #scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed()

sce_sub = sce_1.2[, clust1col %in% c(3,5,6,7)]
de = findMarkers(sce_sub, groups = clust1col[clust1col %in% c(3,5,6,7)])
de2 = SoupX::quickMarkers(counts(sce_sub),clust1col[clust1col %in% c(3,5,6,7)] == 7,N=20)

clust7_umap = cluster(Y= reducedDim(sce_1.2, "TSNE"), positions = as.matrix(positions), q = 7, init = km7, nrep = 15000, gamma = 2)
plot(clust7_umap$plogLik, type = "l")
clust7col_umap = apply(clust7_umap$z[10000:15000,], 2, Mode)
ggplot(as.data.frame(colData(sce_1.2)), aes(x = x, y = -y, fill = factor(clust7col_umap))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()


ptm = proc.time()
cluster(Y= reducedDim(sce_1.2, "TSNE"), positions = as.matrix(positions), q = 7, init = km7, nrep = 500, gamma = 2)
proc.time-ptm()

#Deconvolution

deconv1 = run_mcmc_squaredeconv(df, nrep = 10000, q = 7, d = 15, seed = 100, prev = clust7)
deconv1 = readRDS("C:/Users/Edward Zhao/Downloads/deconv1_small.RDS")
plot(deconv3$plogLik, type = "l")
plot(deconv3$mu[,1+0*15], type = "l")
plot(deconv3$mu[,1+2*15], type = "l")
plot(sapply(deconv1$Y,"[[", 2000), type = "l")
deconv1_col = apply(deconv1$z[25000:30000,], 2, Mode)
deconv1_alpha = pmax(colMeans(deconv1$z[25000:30000,]==1),
                     colMeans(deconv1$z[25000:30000,]==2),
                     colMeans(deconv1$z[25000:30000,]==3),
                     colMeans(deconv1$z[25000:30000,]==4),
                     colMeans(deconv1$z[25000:30000,]==5),
                     colMeans(deconv1$z[25000:30000,]==6),
                     colMeans(deconv1$z[25000:30000,]==7))
Ymeans = Reduce('+', deconv1$Y[250:292])/length(252:292)
ggplot(as.data.frame(df2), aes(x = x, y = -y, 
                               alpha = deconv1_alpha,
                               fill = factor(deconv1_col))) + 
  geom_point(size = 2.1, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+ 
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #scale_fill_viridis_d(option = "A")+
  guides(alpha = F) +
  theme_void()+ coord_fixed()

ggplot(as.data.frame(df2), aes(x = x, y = -y, 
                               fill = Ymeans[,1])) + 
  geom_point(size = 2.1, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+ 
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #scale_fill_viridis_d(option = "A")+
  guides(alpha = F) +
  theme_void()+ coord_fixed()
ggplot(as.data.frame(df2[rep(clust7col %in% c(3,5,6,7),9),]), aes(x = x, y = -y, fill = factor(deconv1_col, levels = c(3,1,2,4)))) + 
  geom_point(size = 2.1, pch = 22)+
  #geom_point(data = as.data.frame(colData(sce_1.2)), aes(x = x, y = -y), size = 5.5, pch = 22, inherit.aes = F, stroke = 1)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()
saveRDS(clust7, "clust7.RDS")
deconv1 = readRDS("data-raw/deconv1_gamma6_small.RDS")
deconv1 = readRDS("data-raw/deconv1_gamma4.RDS")
deconv5 = readRDS("data-raw/deconv2_gamma2.RDS")

deconv4 = deconvolve(Y= PCs$components[,1:15], positions = positions, q = 7, init = rep(clust7col, 9), nrep = 5000, gamma = 1, dist = 0.34)
plot(deconv4$mu[,1+15*4], type = "l")
plot(deconv4$mu[,1+15*5], type = "l")
plot(deconv4$plogLik, type = "l")
plot(sapply(deconv4$lambda,"[", 3,3), type = "l")

deconv4_col = apply(deconv4$z[4000:5000,], 2, Mode)
ggplot(as.data.frame(df2), aes(x = x, y = -y, fill = factor(deconv1_col))) + 
  geom_point(size = 2.1, pch = 22)+
  #geom_point(data = as.data.frame(colData(sce_1.2)), aes(x = x, y = -y), size = 5.5, pch = 22, inherit.aes = F, stroke = 1)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()
deconv4_alpha = pmax(colMeans(deconv1$z[40000:50000,]==1),
                     colMeans(deconv1$z[40000:50000,]==2),
                     colMeans(deconv1$z[40000:50000,]==3),
                     colMeans(deconv1$z[40000:50000,]==4))


plot(deconv5$mu[,1+15*0], type = "l")
plot(deconv5$mu[,1+15*5], type = "l")
plot(deconv5$Ychange, type = "l")
plot(sapply(deconv5$lambda,"[", 3,3), type = "l")

deconv5_col = apply(deconv5$z[45000:50000,], 2, Mode)
test1 = apply(deconv5$z[seq(100,nrow(deconv5$z), 100),], 2, Mode)
test2 = apply(deconv1$z[seq(100,nrow(deconv1$z), 100),], 2, Mode)
test3 = apply(clust7$z[95000:100000,], 2, Mode)
test4 = apply(clust7$z[55000:60000,], 2, Mode)
plot(deconv5$mu[seq(1000,50000,1000),1], type = "l")
ggplot(as.data.frame(df2), aes(x = x, y = -y, fill = factor(deconv1$z[1,]))) + 
  geom_point(size = 2.1, pch = 22)+
  #geom_point(data = as.data.frame(colData(sce_1.2)), aes(x = x, y = -y), size = 5.5, pch = 22, inherit.aes = F, stroke = 1)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()
ggplot(as.data.frame(df2), aes(x = x, y = -y, fill = factor(test2))) + 
  geom_point(size = 2.1, pch = 22)+
  #geom_point(data = as.data.frame(colData(sce_1.2)), aes(x = x, y = -y), size = 5.5, pch = 22, inherit.aes = F, stroke = 1)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()
