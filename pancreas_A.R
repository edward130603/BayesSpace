library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(scran)
library(viridis)
library(patchwork)

data1 = read.delim("data-raw/PDAC/GSM3036911_PDAC-A-ST1-filtered.txt", stringsAsFactors=FALSE)
#data2 = read.delim("data-raw/PDAC/GSM3036911.tsv", stringsAsFactors=FALSE)
positions = data1$X
data1 = t(as.matrix(data1[,-1]))
positions = data.frame(positions) %>%
  separate(positions, c("x","y"), "x") %>%
  mutate_all(as.numeric)


ggplot(positions, aes(x, y, col = data1["CRISP3",]))+
  geom_point(size = 5)+coord_fixed()

#Make sce
sce_A = SingleCellExperiment(assays = list(counts = data1), colData = positions)

#Preprocess
df_qc = perCellQCMetrics(sce_A)
qc = quickPerCellQC(df_qc)
ggplot(positions, aes(x, -y, col = qc$discard))+
  geom_point(size = 5)+coord_fixed()
sce_A$discard = qc$discard  #not actually discarded.. but just marked as low quality

#Normalization
set.seed(100)

sce_A = logNormCounts(sce_A)

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


ggplot(as.data.frame(colData(sce_A)), aes(x = x, y = -y, fill = logcounts(sce_A)[grep("S100A4", rownames(sce_A)),])) + 
  geom_point(pch = 22, size = 7)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed()
ggplot(as.data.frame(colData(sce_A)), aes(x = x, y = y, fill = reducedDim(sce_A, "PCA")[,1])) + 
  geom_point(pch = 22, size = 7)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed()
ggplot(as.data.frame(colData(sce_A)), aes(x = x, y = y, fill = sizeFactors(sce_A))) + 
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
km7 = kmeans(reducedDim(sce_A, "PCA"), centers = 7)$cluster
km8 = kmeans(reducedDim(sce_A, "PCA"), centers = 8)$cluster
km5 = kmeans(reducedDim(sce_A, "PCA"), centers = 5)$cluster
km9 = kmeans(reducedDim(sce_A, "PCA"), centers = 9)$cluster
km10 = kmeans(reducedDim(sce_A, "PCA"), centers = 10)$cluster
clust4 = cluster(Y= PCs$components[,1:10], positions = as.matrix(positions), q = 4, init = km4, nrep = 10000, gamma = 1.5, dist = 1)
clust4_g3 = cluster(Y= PCs$components[,1:10], positions = as.matrix(positions), q = 4, init = km4, nrep = 10000, gamma = 3, dist = 1)
clust5 = cluster(Y= PCs$components[,1:10], positions = as.matrix(positions), q = 5, init = km5, nrep = 10000, gamma = 2, dist = 1)
clust6 = cluster(Y= PCs$components[,1:10], positions = as.matrix(positions), q = 6, init = km6, nrep = 10000, gamma = 2, dist = 1)
clust8 = cluster(Y= PCs$components[,1:10], positions = as.matrix(positions), q = 8, init = km8, nrep = 10000, gamma = 3, dist = 1)
clust9 = cluster(Y= PCs$components[,1:10], positions = as.matrix(positions), q = 9, init = km9, nrep = 10000, gamma = 3, dist = 1)

plot(clust4$plogLik, type = "l")
clust4col = apply(clust4$z[9000:10000,], 2, Mode)
clust4_alpha = pmax(colMeans(clust4$z[9000:10000,]==1),
                     colMeans(clust4$z[9000:10000,]==2),
                     colMeans(clust4$z[9000:10000,]==3),
                     colMeans(clust4$z[9000:10000,]==4))
clust4g3col = apply(clust4_g3$z[9000:10000,], 2, Mode)
clust5col = apply(clust5$z[9000:10000,], 2, Mode)
clust6col = apply(clust6$z[9000:10000,], 2, Mode)
clust8col = apply(clust8$z[9000:10000,], 2, Mode)
clust9col = apply(clust9$z[9000:10000,], 2, Mode)
ggplot(as.data.frame(colData(sce_A)), aes(x = x, y = y, fill = factor(km4))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, color = "Cluster") +
  scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()
ggplot(as.data.frame(colData(sce_A)), aes(x = x, y = y, fill = factor(clust9col))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, color = "Cluster") +
  scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()
ggplot(as.data.frame(colData(sce_A)), aes(x = x, y = y, fill = factor(clust4col), alpha = clust4_alpha)) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  guides(alpha = F)+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+ 
  scale_fill_manual(values = c("red", "green", "yellow", "blue"))+
  #scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()

#Deconvolution
df = data.frame(PCs$components[,1:10])
colnames(df) = paste0("Y",1:10)
df$x = positions$x
df$y = positions$y
df$j = 1:nrow(df)
saveRDS(clust4, "data-raw/clust4_pancreasA.RDS")
saveRDS(df, "data-raw/df_pancreasA.RDS")
deconv1 = deconvolve(Y= PCs$components[,1:10], positions = positions, q = 4, init = rep(clust4col, 9), nrep = 5000, gamma = 1, dist = 0.34)
plot(deconv1$mu[,1+15*0], type = "l")
plot(deconv1$plogLik, type = "l")
plot(sapply(deconv1$lambda,"[", 3,3), type = "l")
positions2 = positions[rep(seq_len(nrow(positions)), 9), ] #rbind 9 times
shift = rbind(expand.grid(c(1/3, -1/3, 0), c(1/3,-1/3,0)))
shift_long = shift[rep(seq_len(9), each = nrow(positions)), ]
positions2[,"x"] = positions2[,"x"] + shift_long[,"Var1"]
positions2[,"y"] = positions2[,"y"] + shift_long[,"Var2"]

deconv1_col = apply(deconv1$z[4000:5000,], 2, Mode)
deconv1_alpha = pmax(colMeans(deconv1$z[4000:5000,]==1),
                     colMeans(deconv1$z[4000:5000,]==2),
                     colMeans(deconv1$z[4000:5000,]==3),
                     colMeans(deconv1$z[4000:5000,]==4))
ggplot(positions2, aes(x = x, y = y, fill = factor(deconv1_col), alpha = deconv1_alpha+ 0.25)) + 
  geom_point(size = 2.1, pch = 22)+
  #geom_point(data = as.data.frame(colData(sce_1.2)), aes(x = x, y = -y), size = 5.5, pch = 22, inherit.aes = F, stroke = 1)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("red", "green", "yellow", "blue"))+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+ 
  guides(alpha= F) +
  #scale_fill_viridis_d(option = "A")+
  theme_void()+ coord_fixed()
deconv4_alpha = pmax(colMeans(deconv1$z[40000:50000,]==1),
                     colMeans(deconv1$z[40000:50000,]==2),
                     colMeans(deconv1$z[40000:50000,]==3),
                     colMeans(deconv1$z[40000:50000,]==4))
#deconv1 = run_mcmc_squaredeconv(df, nrep = 100, q = 4, d = 10, seed = 100, prev = clust4)
deconv2 = readRDS("data-raw/deconv1.RDS")
plot(deconv2$mu[,1+0*15], type = "l")
plot(deconv2$mu[,1+2*15], type = "l")
plot(sapply(deconv2$Y,"[[", 2000), type = "l")
deconv2_col = apply(deconv2$z[50000:60000,], 2, Mode)
deconv2_alpha = pmax(colMeans(deconv2$z[seq(100,60000,100),]==1),
                     colMeans(deconv2$z[seq(100,60000,100),]==2),
                     colMeans(deconv2$z[seq(100,60000,100),]==3),
                     colMeans(deconv2$z[seq(100,60000,100),]==4))
Ymeans = Reduce('+', deconv1$Y[250:292])/length(252:292)
ggplot(as.data.frame(positions2), aes(x = x, y = y, 
                               #alpha = deconv2_alpha,
                               fill = factor(deconv2_col))) + 
  geom_point(size = 2.1, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.3,1,0.1), range = c(0,1))+ 
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #scale_fill_viridis_d(option = "A")+
  scale_fill_manual(values = c("red", "green", "yellow", "blue"))+
  guides(alpha = F) +
  theme_void()+ coord_fixed()
