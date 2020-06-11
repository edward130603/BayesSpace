library(tidyverse)
library(DropletUtils)
library(scater)
library(scran)
library(viridis)
library(patchwork)

#Load sce
#rep1
sce_A = read10xCounts("data-raw/breast/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5", type = "HDF5")
pos = read.csv("data-raw/breast/spatial1/tissue_positions_list.csv", header=FALSE)
#rep2
sce_A = read10xCounts("data-raw/breast/V1_Breast_Cancer_Block_A_Section_2_filtered_feature_bc_matrix.h5", type = "HDF5")
pos = read.csv("data-raw/breast/spatial2/tissue_positions_list.csv", header=FALSE)
colnames(pos) = c("Barcode", "tissue", "Y1", "X1", "Y2", "X2")
xdist = coef(lm(data = pos, X2~X1))[2]
ydist = coef(lm(data = pos, Y2~Y1))[2]
dist = xdist + ydist + 5
positions = pos[,c("X2","Y2")]
colnames(positions) = c("x", "y")
rownames(positions) = pos$Barcode
positions = positions[sce_A$Barcode,]
sce_A$x = positions$x
sce_A$y = positions$y

#Preprocess
df_qc = perCellQCMetrics(sce_A)
qc = quickPerCellQC(df_qc)
ggplot(positions, aes(x, -y, col = qc$discard))+
  geom_point(size = 3)+coord_fixed()
sce_A$discard = qc$discard #not actually discarded

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

#Run mcmc
km = sapply(1:30,function(k){kmeans(PCs$components, centers = k)$cluster})
library(pbapply)
clust_list = pblapply(2:15, function(k){cluster(Y= PCs$components, positions = as.matrix(positions), q = k, init = km[,k], nrep = 10000, gamma = 1.5, dist = dist)})
clust_list2 = pblapply(16:30, function(k){cluster(Y= PCs$components, positions = as.matrix(positions), q = k, init = km[,k], nrep = 10000, gamma = 1.5, dist = dist)})
saveRDS(clust_list, "data-raw/clust_list.RDS")
saveRDS(clust_list2, "data-raw/clust_list2.RDS")
clust_list = readRDS("data-raw/clust_list.RDS")

plogLiks = sapply(clust_list, "[[", 4)
plogLiks2 = sapply(clust_list2, "[[", 4)
plogLiks = cbind(plogLiks, plogLiks2)
plot(x = 2:30, y = colMeans(plogLiks[9000:10000,]), type = "b", ylab = "Pseudo-log-likelihood")
plot(x = 2:30, y = plogLiks[10000,] -0.5* 11*(2:30)*log(nrow(positions)), type = "b", ylab = "Pseudo-BIC")

#FMM 
library(mclust)
clust_bic = Mclust(PCs$components, G = 1:20, modelNames = c("VVV"))
clust_icl = mclustICL(PCs$components, G = 1:20, modelNames = c("VVV"))

ggplot(positions, aes(x = x, y = -y))+
  geom_text(aes(color =  logcounts(sce_A)["ERBB2",]),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  coord_fixed() + theme_void() + guides(col = F) + scale_color_viridis(option = "A")
