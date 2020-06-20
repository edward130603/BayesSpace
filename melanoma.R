###Load packages, script
library(tidyverse)
library(SingleCellExperiment)
library(patchwork)
library(scater)
library(mclust)
library(scran)
library(teigen)
library(xgboost)
source("script.R")

###Load SCEs, generate input data structures
melanoma1.1 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel1_rep1.RDS")
melanoma1.2 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel1_rep2.RDS")

intersect_genes = intersect(rownames(melanoma1.1), rownames(melanoma1.2))
melanoma1.1 = melanoma1.1[intersect_genes,]
melanoma1.2 = melanoma1.2[intersect_genes,]
set.seed(100)
pca_1.2 = BiocSingular::runPCA(t(logcounts(melanoma1.2)), rank = 50)
reducedDim(melanoma1.2, "PCA") = pca_1.2$x
reducedDim(melanoma1.1, "PCA") =  scale(t(logcounts(melanoma1.1)), scale = F)%*% pca_1.2$rotation

positions1.1 = cbind(x = - melanoma1.1$col, y = - melanoma1.1$row)
positions1.2 = cbind(x = melanoma1.2$col, y = melanoma1.2$row)

Y1.1 = reducedDim(melanoma1.1)[,1:10]
Y1.2 = reducedDim(melanoma1.2)[,1:10]

plot1.1 = ggplot(data.frame(positions1.1), aes(x = x, y = y, fill = librarySizeFactors(melanoma1.1))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Library Size") +
  scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed()
plot1.2 = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = librarySizeFactors(melanoma1.2))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Library Size") +
  scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed()

###Mclust
set.seed(100)
mclust1.1 = Mclust(Y1.1, 4, "EEE")$classification
mclust1.2 = Mclust(Y1.2, 4, "EEE")$classification

plot_mclust1.1 = ggplot(data.frame(positions1.1), aes(x = x, y = y, fill = factor(mclust1.1))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("yellow", "red", "purple", "blue"))+
  guides(fill = F)+
  theme_void()+ coord_fixed()
plot_mclust1.2 = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = factor(mclust1.2))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("purple", "yellow", "red", "blue"))+
  guides(fill = F)+
  theme_void()+ coord_fixed()
plot_mclust1.1 + plot_mclust1.2

###Spatial clustering
clust1.1 = cluster(Y= Y1.1, positions = positions1.1, q = 4, model = "t", init = mclust1.1, nrep = 10000, gamma = 2, dist = 1)
clust1.2 = cluster(Y= Y1.2, positions = positions1.2, q = 4, model = "t", init = mclust1.2, nrep = 10000, gamma = 2, dist = 1)
# saveRDS(clust1.1, "data-raw/clust1.1_q4_figure4b.RDS")
# saveRDS(clust1.2, "data-raw/clust1.2_q4_figure4b.RDS")
# clust1.1 = readRDS("data-raw/clust1.1_q4_figure4b.RDS")
# clust1.2 = readRDS("data-raw/clust1.2_q4_figure4b.RDS")
clust1.1_alpha = pmax(colMeans(clust1.1$z[seq(1000,10000,1),]==1),
                      colMeans(clust1.1$z[seq(1000,10000,1),]==2),
                      colMeans(clust1.1$z[seq(1000,10000,1),]==3),
                      colMeans(clust1.1$z[seq(1000,10000,1),]==4))
clust1.2_alpha = pmax(colMeans(clust1.2$z[seq(1000,10000,1),]==1),
                      colMeans(clust1.2$z[seq(1000,10000,1),]==2),
                      colMeans(clust1.2$z[seq(1000,10000,1),]==3),
                      colMeans(clust1.2$z[seq(1000,10000,1),]==4))

plot_clust1.1 = ggplot(data.frame(positions1.1), aes(x = x, y = y, alpha = clust1.1_alpha, fill = factor(clust1.1$labels))) +
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("yellow", "red", "purple", "blue"))+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+
  guides(alpha = F, fill = F) +
  theme_void()+ coord_fixed()
plot_clust1.2 = ggplot(data.frame(positions1.2), aes(x = x, y = y, alpha = clust1.2_alpha, fill = factor(clust1.2$labels))) +
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("purple", "yellow", "red", "blue"))+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+
  guides(alpha = F, fill = F) +
  theme_void()+ coord_fixed()

plot_clust1.1+plot_clust1.2

###Enhance
deconv1.1 = deconvolve(Y= Y1.1, positions = positions1.1, q = 4, init = clust1.1$labels,  model = "t",
                      nrep = 200000, gamma = 2, xdist = 1, ydist = 1, platform = "ST", c = 0.003, jitter_scale = jitter_scale)
deconv1.1col = apply(deconv1.2$z[seq(10000,200000,10),], 2, Mode)
deconv1.1_alpha = pmax(colMeans(deconv1.1$z[seq(10000,200000,10),]==1),
                       colMeans(deconv1.1$z[seq(10000,200000,10),]==2),
                       colMeans(deconv1.1$z[seq(10000,200000,10),]==3),
                       colMeans(deconv1.1$z[seq(10000,200000,10),]==4))
deconv1.1 = list(obj = deconv1.1, cols = deconv1.1col, alpha = deconv1.1_alpha)
deconv1.2 = deconvolve(Y= Y1.2, positions = positions1.2, q = 4, init = clust1.2$labels,  model = "t",
                      nrep = 200000, gamma = 2, xdist = 1, ydist = 1, platform = "ST", c = 0.003, jitter_scale = 8)
deconv1.2col = apply(deconv1.2$z[seq(10000,200000,10),], 2, Mode)
deconv1.2_alpha = pmax(colMeans(deconv1.2$z[seq(10000,200000,10),]==1),
                       colMeans(deconv1.2$z[seq(10000,200000,10),]==2),
                       colMeans(deconv1.2$z[seq(10000,200000,10),]==3),
                       colMeans(deconv1.2$z[seq(10000,200000,10),]==4))
deconv1.2 = list(obj = deconv1.2, cols = deconv1.2col, alpha = deconv1.2_alpha)
# deconv1.1 = readRDS("data-raw/deconv1.1_6-17_c0.003_small.RDS") 
# deconv1.2 = readRDS("data-raw/deconv1.2_6-17_c0.003_small.RDS")

deconv1.1Y = Reduce(`+`, deconv1.1$obj$Y[-(1:100)])/length(c(1:length(deconv1.1$obj$Y))[-(1:100)])
deconv1.2Y = Reduce(`+`, deconv1.2$obj$Y[-(1:100)])/length(c(1:length(deconv1.2$obj$Y))[-(1:100)])

plot_deconv1.1 = ggplot(as.data.frame(deconv1.1$obj$positions),
                        aes(x = x, y = y, alpha = deconv1.2_alpha, fill = factor(deconv1.1col))) +
  geom_point(size = 2, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  guides(alpha = F, fill = F) +
  scale_fill_manual(values = c("red", "blue", "purple", "yellow"))+
  theme_void()+ coord_fixed()
plot_deconv1.2 = ggplot(as.data.frame(deconv1.2$obj$positions),
       aes(x = x, y = y, alpha = deconv1.2_alpha, fill = factor(deconv1.2col))) +
  geom_point(size = 2, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  guides(alpha = F, fill = F) +
  scale_fill_manual(values = c("red", "blue", "purple", "yellow"))+
  theme_void()+ coord_fixed()

###Predict gene expressions
{
#Tune xgboost
xgboost_hyper = matrix(nrow = length(HVG), ncol = 500)
rownames(xgboost_hyper) = HVG
lm_r2 = numeric(length(HVG))
names(lm_r2) = HVG
for(gene in HVG){
  Y1.2train = xgb.DMatrix(data = Y1.2, label = logcounts(melanoma1.2)[gene,])
  Y1.2test = xgb.DMatrix(data = Y1.1, label = logcounts(melanoma1.1)[gene,])
  ntest = nrow(Y1.2test)
  watchlist = list(train = Y1.2train, test = Y1.2test)
  train = xgb.train(data = Y1.2train, max_depth = 2, watchlist = watchlist,
                    eta = 0.03, nthread = 2, nrounds = 500, objective = "reg:squarederror", verbose = F)
  sstot = sum((mean(logcounts(melanoma1.1)[gene,]) - logcounts(melanoma1.1)[gene,])^2)
  xgboost_hyper[gene,] = 1 - (train$evaluation_log$test_rmse)^2 * ntest/sstot
  
  reg = lm(logcounts(melanoma1.2)[gene,] ~ ., data = data.frame(Y1.2))
  pred_lm = predict(reg, data.frame(Y1.1))
  lm_ssresid = sum((pred_lm - logcounts(melanoma1.1)[gene,])^2)
  lm_r2[gene] = 1 - lm_ssresid/sstot
  print(which(HVG == gene))
}
#Predict xgboost
dec1 <- modelGeneVarByPoisson(melanoma1.1)
top1 <- getTopHVGs(dec1, prop=0.1)
dec2 <- modelGeneVarByPoisson(melanoma1.2)
top2 <- getTopHVGs(dec2, prop=0.1)
markers = c("CD2", "CD3D", "CD3E", "CD3G", "CD19", "CD79A", "CD79B", "BLK",
            "CD163", "CD14", "CSF1R", "PECAM1", "VWF", "CDH5", "FAP", "THY1", "DCN",
            "MIA", "TYR", "PMEL", "MLANA", "PRAME")
HVG = union(union(top1, top2), markers)

xgboost_expression = matrix(ncol = nrow(deconv1.2Y), nrow = length(HVG))
rownames(xgboost_expression) = HVG
for (gene in HVG){
  train = xgboost(data = Y1.2, label = logcounts(melanoma1.2)[gene,], max_depth = 2,
                  eta = 0.03, nthread = 2, nrounds = 100, objective = "reg:squarederror", verbose = F)
  xgboost_expression[gene, ] = predict(train, deconv1.2Y)
}
xgboost_expression = xgboost_expression[xgboost_hyper2[,100]> 0, ]
saveRDS(xgboost_expression, "data-raw/xgboost_1.2.RDS")

colnames(xgboost_hyper ) = 1:500
reshape::melt(xgboost_hyper) %>%
  as_tibble() %>%
  filter(X2 %in% seq(1, 500, 10)) %>%
  ggplot(aes(x = X2, y = value, group = X2))+
  geom_boxplot(width = 5)

set.seed(1)
traindata = sample(1:293, 200)

xgboost_hyper2 = matrix(nrow = length(HVG), ncol = 500)
rownames(xgboost_hyper2) = HVG
lm_r2_v2 = numeric(length(HVG))
names(lm_r2_v2) = HVG

for(gene in HVG){
  Y1.2train = xgb.DMatrix(data = Y1.2[traindata,], label = logcounts(melanoma1.2)[gene,traindata])
  Y1.2test = xgb.DMatrix(data = Y1.2[-traindata,], label = logcounts(melanoma1.2)[gene,-traindata])
  ntest = nrow(Y1.2[-traindata,])
  watchlist = list(train = Y1.2train, test = Y1.2test)
  train = xgb.train(data = Y1.2train, max_depth = 2, watchlist = watchlist,
                    eta = 0.03, nthread = 2, nrounds = 500, objective = "reg:squarederror", verbose = F)
  sstot = sum((mean(logcounts(melanoma1.2)[gene,-traindata]) - logcounts(melanoma1.2)[gene,-traindata])^2)
  xgboost_hyper2[gene,] = 1 - (train$evaluation_log$test_rmse)^2 * ntest/sstot
  
  reg = lm(logcounts(melanoma1.2)[gene,traindata] ~ ., data = data.frame(Y1.2[traindata,]))
  pred_lm = predict(reg, data.frame(Y1.2[-traindata,]))
  lm_ssresid = sum((pred_lm - logcounts(melanoma1.2)[gene,-traindata])^2)
  lm_r2_v2[gene] = 1 - lm_ssresid/sstot
  print(which(HVG == gene))
}

}

###Clustering Simulation
clust1.1mu = matrix(colMeans(clust1.1$mu[-(1:1000),]), byrow = T, ncol = 10)
clust1.1lambda = Reduce(`+`, clust1.1$lambda[seq(1000,10000, 1)])/length(seq(1000,10000, 1))

clust1.2mu = matrix(colMeans(clust1.2$mu[-(1:1000),]), byrow = T, ncol = 10)
clust1.2lambda = Reduce(`+`, clust1.2$lambda[seq(1000,10000, 1)])/length(seq(1000,10000, 1))
simY = matrix(nrow = nrow(positions1.1), ncol = 10)
##next chunk run on cluster
{
set.seed(100)
for(i in 1:nrow(simY)){
  simY[i,] = mvnfast::rmvt(1, mu = clust1.2mu[clust1.2$labels[i],], sigma = solve(clust1.2lambda), df = 4)
  save()
}
ggplot(data.frame(positions1.2), aes(x, y, fill = simY[,1])) +
  geom_point(size = 7, psich = 22)+
  labs(fill = "PC1", alpha = "Proportion", x = NULL, y = NULL) +
  theme_void() + coord_fixed()

set.seed(2)
for(i in 1:nrow(positions1.2)){
  simY[i,] = mvnfast::rmvt(1, mu = clust1.2mu[clust1.2$labels[i],], sigma = solve(clust1.2lambda), df = 4)
}
saveRDS(simY, paste0("clust1.2sim/clust1.2_sim_", sid, ".RDS"))
ARIs = numeric(6)
names(ARIs) = c("km", "louvain", "mclust", "spatial_clust", "spatial_clust_t", "spatial_clust_t_true")

#K means
set.seed(101)
km_sim = kmeans(simY, centers = 4)$cluster
ARIs["km"] = mclust::adjustedRandIndex(clust1.2$labels, km_sim)

#louvain
set.seed(100)
reducedDim(melanoma1.2, "PCA_sim") = simY
g.jaccard <- buildSNNGraph(melanoma1.2, use.dimred="PCA_sim", type="jaccard")
clust.louvain <- igraph::cluster_louvain(g.jaccard)$membership
ARIs["louvain"] = mclust::adjustedRandIndex(clust1.2$labels, clust.louvain)

#Mclust
set.seed(102)
mclust_sim = Mclust(simY, G = 4, modelNames = c("EEE"))$classification
ARIs["mclust"] = mclust::adjustedRandIndex(clust1.2$labels, mclust_sim)

#teigen
# set.seed(103)
# teigen_sim = teigen(simY, Gs = 4, models = "CCCC")$classification
# if (is.null(teigen_sim)){
#   teigen_sim = teigen(simY, Gs = 4, models = "CCCC", init = "hard")$classification
# }
# ARIs["teigen"] = mclust::adjustedRandIndex(clust1.2$labels, teigen_sim)

#Spatial
clust_sim = cluster(Y= simY, positions = positions1.2, q = 4, init = mclust_sim, 
                    nrep = 10000, gamma = 2, dist = 1, model = "normal")
ARIs["spatial_clust"] = mclust::adjustedRandIndex(clust1.2$labels, clust_sim$labels)

clust_sim_t = cluster(Y= simY, positions = positions1.2, q = 4, init = mclust_sim, 
                      nrep = 10000, gamma = 2, dist = 1, model = "t")
ARIs["spatial_clust_t"] = mclust::adjustedRandIndex(clust1.2$labels, clust_sim_t$labels)

clust_sim_t_true = cluster(Y= simY, positions = positions1.2, q = 4, init = clust1.2$labels,
                           nrep = 10000, gamma = 2, dist = 1, model = "t")
ARIs["spatial_clust_t_true"] = mclust::adjustedRandIndex(clust1.2$labels, clust_sim_t_true$labels)

ari_files = list.files("data-raw/melanoma_clust_sim/")
ari_out = matrix(nrow = 30, ncol = 7)
colnames(ari_out) = names(readRDS(paste0("data-raw/melanoma_clust_sim/", ari_files[1])))
rownames(ari_out) = 1:30
for (i in 1:length(ari_files)){
  ari_out[i,] = readRDS(paste0("data-raw/melanoma_clust_sim/", ari_files[i]))  
}
ari_out = ari_out[,-4]
ari_out_long = reshape::melt(ari_out)
ari_out_long$Sample = "2"
ari_plot = ggplot(ari_out_long, aes(x = factor(X2, levels = colnames(ari_out)), y = value, group = X2, col = factor(X1)))+ #figure 4d
  geom_boxplot()+ggbeeswarm::geom_quasirandom()+
  theme_light()+guides(col = F)+
  scale_x_discrete(labels = c("k-means", "Louvain", "mclust", "SC (normal)", "SC (t)", 
                              "SC (t, truth init.)"))+
  labs(x = NULL, y = "ARI")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ari_files1.1 = list.files("data-raw/melanoma_clust_sim_1.1/")
ari_out1.1 = matrix(nrow = 30, ncol = 6)
colnames(ari_out1.1) = names(readRDS(paste0("data-raw/melanoma_clust_sim_1.1/", ari_files1.1[1])))
rownames(ari_out1.1) = 1:30
for (i in 1:length(ari_files1.1)){
  ari_out1.1[i,] = readRDS(paste0("data-raw/melanoma_clust_sim_1.1/", ari_files1.1[i]))  
}
ari_out_long1.1 = reshape::melt(ari_out1.1)
ari_out_long1.1$Sample = "1"
rbind(ari_out_long, ari_out_long1.1) %>%
  ggplot(aes(col = factor(X2, levels = colnames(ari_out)), y = value, x = Sample))+ #figure 4d
  geom_boxplot(aes(group = interaction(X2, Sample)), color = "black", position = position_dodge(0.7), width = 0.6, outlier.shape = NA)+
  #geom_point(position = position_dodge(0.7))+
  ggbeeswarm::geom_quasirandom(dodge.width = 0.7)+
  theme_light()+
  scale_color_discrete(labels = c("k-means", "Louvain", "mclust", "SC (normal)", "SC (t)", 
                              "SC (t, truth init.)"))+
  scale_y_continuous(breaks = seq(0.4, 1, 0.1), minor_breaks = seq(0.4,1,0.05), limits = c(0.4,1))+
  labs(x = "Replicate", y = "ARI", col = "Method")#+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
}

###Deconvolution Simulation
deconv1.1mu = matrix(colMeans(deconv1.1$obj$mu[100:2000,]), byrow = T, ncol = 10)
deconv1.1lambda = Reduce(`+`, deconv1.1$obj$lambda[100:2000])/length(100:2000)
deconv1.1Y = Reduce(`+`, deconv1.1$obj$Y[-(1:100)])/length(c(1:length(deconv1.1$obj$Y))[-(1:100)])
deconv1.1positions = deconv1.1$obj$positions
deconv1.1col = deconv1.1$cols
deconv1.1alpha = deconv1.1$alpha

deconv1.2mu = matrix(colMeans(deconv1.2$obj$mu[100:2000,]), byrow = T, ncol = 10)
deconv1.2lambda = Reduce(`+`, deconv1.2$obj$lambda[100:2000])/length(100:2000)
deconv1.2Y = Reduce(`+`, deconv1.2$obj$Y[-(1:100)])/length(c(1:length(deconv1.2$obj$Y))[-(1:100)])
deconv1.2positions = deconv1.2$obj$positions
deconv1.2col = deconv1.2$cols
deconv1.2alpha = deconv1.2$alpha

##next chunk  run on cluster 
{ 
simY1.1 = matrix(nrow = nrow(deconv1.1positions), ncol = 10)
set.seed(1)
for(i in 1:nrow(deconv1.1positions)){
  simY1.1[i,] = mvnfast::rmvt(1, mu = deconv1.1mu[deconv1.1col[i],], sigma = solve(deconv1.1lambda), df = 4)
}

data_sim1.1 = data.frame(simY1.1)
data_sim1.1$j = rep(1:nrow(Y1.1), 9)
data_sim1.1 %>%
  group_by(j) %>%
  summarise_all(mean)->
  data_sim_mean1.1

ARIs = numeric(9)
names(ARIs) = c("km", "louvain", "mclust", "spatial_clust_t", "deconv_t_c0.003", "deconv_t_c0.01", 
                "spot_rmse", "deconv_t_c0.003_rmse","deconv_t_c0.01_rmse")

set.seed(101)
km_sim = kmeans(data_sim_mean1.1[,-1], centers = 4)$cluster
ARIs["km"] = mclust::adjustedRandIndex(deconv1.1col, rep(km_sim,9))
set.seed(103)
#louvain
set.seed(100)
sim1.1 = melanoma1.1[,rep(seq_len(nrow(Y1.1)),9)]
reducedDim(sim1.1, "PCA_sim") = simY1.1
g.jaccard <- buildSNNGraph(sim1.1, use.dimred="PCA_sim", type="jaccard")
clust.louvain <- igraph::cluster_louvain(g.jaccard)$membership
ARIs["louvain"] = mclust::adjustedRandIndex(deconv1.1col, clust.louvain)

data_sim_mean1.1long = as.matrix(data_sim_mean1.1[rep(1:nrow(data_sim_mean1.1),9),-1])
ARIs["spot_rmse"] = mean((data_sim_mean1.1long - simY1.1)^2)


set.seed(102)
mclust_sim = Mclust(data_sim_mean1.1[,-1], G = 4, modelNames = c("EEE"))$classification
ARIs["mclust"] = mclust::adjustedRandIndex(deconv1.1col, rep(mclust_sim,9))

clust_sim = cluster(Y= data_sim_mean1.1[,-1], positions = as.matrix(positions1.1), q = 4, init = mclust_sim, 
                    nrep = 10000, gamma = 2, dist = 1, model = "t")
ARIs["spatial_clust_t"] = mclust::adjustedRandIndex(deconv1.1col, rep(clust_sim$labels,9))
clust_col = clust_sim$labels
rm(clust_sim)
deconv_sim_t = deconvolve(Y= data_sim_mean1.1[,-1], positions = positions1.1, q = 4, init = clust_col, model = "t",
                          nrep = 200000, gamma = 2, xdist = 1, ydist = 1, platform = "ST", c = 0.003, jitter_scale = 1)
deconv_simcol_t = apply(deconv_sim_t$z[seq(10000,200000,10),], 2, Mode)
ARIs["deconv_t_c0.003"] = mclust::adjustedRandIndex(deconv1.1col, deconv_simcol_t)
summary(deconv_sim_t$Ychange)
deconv_run_Y1.1 = Reduce(`+`, deconv_sim_t$Y[-(1:100)])/length(c(1:length(deconv_sim_t$Y))[-(1:100)])
ARIs["deconv_t_c0.003_rmse"] = mean((deconv_run_Y1.1 - simY1.1)^2)

rm(deconv_sim_t)

deconv_sim_t = deconvolve(Y= data_sim_mean1.1[,-1], positions = positions1.1, q = 4, init = clust_col, model = "t",
                          nrep = 200000, gamma = 2, xdist = 1, ydist = 1, platform = "ST", c = 0.01, jitter_scale = 1)
deconv_simcol_t = apply(deconv_sim_t$z[seq(10000,200000,10),], 2, Mode)
ARIs["deconv_t_c0.01"] = mclust::adjustedRandIndex(deconv1.1col, deconv_simcol_t)
summary(deconv_sim_t$Ychange)

deconv_run_Y1.1 = Reduce(`+`, deconv_sim_t$Y[-(1:100)])/length(c(1:length(deconv_sim_t$Y))[-(1:100)])
ARIs["deconv_t_c0.01_rmse"] = mean((deconv_run_Y1.1 - simY1.1)^2)
rm(deconv_sim_t)

filename = paste0("ari_deconv_sim1.1_", sid, ".RDS")
saveRDS(ARIs, filename)


##Sample 1.2
simY1.2 = matrix(nrow = nrow(deconv1.2positions), ncol = 10)
set.seed(1)
for(i in 1:nrow(deconv1.2positions)){
  simY1.2[i,] = mvnfast::rmvt(1, mu = deconv1.2mu[deconv1.2col[i],], sigma = solve(deconv1.2lambda), df = 4)
}

data_sim1.2 = data.frame(simY1.2)
data_sim1.2$j = rep(1:nrow(Y1.2), 9)
data_sim1.2 %>%
  group_by(j) %>%
  summarise_all(mean)->
  data_sim_mean1.2

ARIs = numeric(9)
names(ARIs) = c("km", "louvain", "mclust", "spatial_clust_t", "deconv_t_c0.003", "deconv_t_c0.01", 
                "spot_rmse", "deconv_t_c0.003_rmse","deconv_t_c0.01_rmse")

set.seed(101)
km_sim = kmeans(data_sim_mean1.2[,-1], centers = 4)$cluster
ARIs["km"] = mclust::adjustedRandIndex(deconv1.2col, rep(km_sim,9))
set.seed(103)
#louvain
set.seed(100)
sim1.2 = melanoma1.2[,rep(seq_len(nrow(Y1.2)),9)]
reducedDim(sim1.2, "PCA_sim") = simY1.2
g.jaccard <- buildSNNGraph(sim1.2, use.dimred="PCA_sim", type="jaccard")
clust.louvain <- igraph::cluster_louvain(g.jaccard)$membership
ARIs["louvain"] = mclust::adjustedRandIndex(deconv1.2col, clust.louvain)

data_sim_mean1.2long = as.matrix(data_sim_mean1.2[rep(1:nrow(data_sim_mean1.2),9),-1])
ARIs["spot_rmse"] = mean((data_sim_mean1.2long - simY1.2)^2)


set.seed(102)
mclust_sim = Mclust(data_sim_mean1.2[,-1], G = 4, modelNames = c("EEE"))$classification
ARIs["mclust"] = mclust::adjustedRandIndex(deconv1.2col, rep(mclust_sim,9))

clust_sim = cluster(Y= data_sim_mean1.2[,-1], positions = as.matrix(positions1.2), q = 4, init = mclust_sim, 
                    nrep = 10000, gamma = 2, dist = 1, model = "t")
ARIs["spatial_clust_t"] = mclust::adjustedRandIndex(deconv1.2col, rep(clust_sim$labels,9))
clust_col = clust_sim$labels
rm(clust_sim)
deconv_sim_t = deconvolve(Y= data_sim_mean1.2[,-1], positions = positions1.2, q = 4, init = clust_col, model = "t",
                          nrep = 200000, gamma = 2, xdist = 1, ydist = 1, platform = "ST", c = 0.003, jitter_scale = 4)
deconv_simcol_t = apply(deconv_sim_t$z[seq(10000,200000,10),], 2, Mode)
ARIs["deconv_t_c0.003"] = mclust::adjustedRandIndex(deconv1.2col, deconv_simcol_t)
summary(deconv_sim_t$Ychange)
deconv_run_Y1.2 = Reduce(`+`, deconv_sim_t$Y[-(1:100)])/length(c(1:length(deconv_sim_t$Y))[-(1:100)])
ARIs["deconv_t_c0.003_rmse"] = mean((deconv_run_Y1.2 - simY1.2)^2)

rm(deconv_sim_t)

deconv_sim_t = deconvolve(Y= data_sim_mean1.2[,-1], positions = positions1.2, q = 4, init = clust_col, model = "t",
                          nrep = 200000, gamma = 2, xdist = 1, ydist = 1, platform = "ST", c = 0.01, jitter_scale = 4)
deconv_simcol_t = apply(deconv_sim_t$z[seq(10000,200000,10),], 2, Mode)
ARIs["deconv_t_c0.01"] = mclust::adjustedRandIndex(deconv1.2col, deconv_simcol_t)
summary(deconv_sim_t$Ychange)

deconv_run_Y1.2 = Reduce(`+`, deconv_sim_t$Y[-(1:100)])/length(c(1:length(deconv_sim_t$Y))[-(1:100)])
ARIs["deconv_t_c0.01_rmse"] = mean((deconv_run_Y1.2 - simY1.2)^2)
rm(deconv_sim_t)

filename = paste0("ari_deconv_sim1.2_", sid, ".RDS")
saveRDS(ARIs, filename)
}

ari_files = list.files("data-raw/deconv_sim_1.1/")
ari_out = matrix(nrow = 30, ncol = 9)
colnames(ari_out) = names(readRDS(paste0("data-raw/deconv_sim_1.1/", ari_files[1])))
rownames(ari_out) = 1:30
for (i in 1:length(ari_files)){
  ari_out[i,] = readRDS(paste0("data-raw/deconv_sim_1.1/", ari_files[i]))  
}
ari_out_long = reshape::melt(ari_out[,c(1,3:6)]) %>% 
  filter(X2 %in% c("deconv_t_c0.003", "deconv_t_c0.01"))
modes = numeric(nrow(Y1.1))
for(i in 1:nrow(Y1.1)){
  modes[i] =  Mode(deconv1.1col[i + (0:8)*nrow(Y1.1)])
}

yint1 = adjustedRandIndex(rep(modes,9), deconv1.1col)

ari_plot_1.1 = ggplot(ari_out_long, aes(x = factor(X2), y = value, group = factor(X2), col = factor(X1)))+ #figure 4d
  geom_boxplot(outlier.size = 0, outlier.stroke = 0)+ggbeeswarm::geom_quasirandom()+
  geom_hline(color = "red", aes(yintercept = yint1, linetype = "Optimal spot-level clustering"), show.legend = T) + 
  theme_light()+guides(col = F)+
  scale_x_discrete(labels = c("Deconvolution (c = 0.003)",
                              "Deconvolution (c = 0.01)"))+
  scale_linetype_manual(values = 2) + 
  labs(x = NULL, y = "ARI", linetype = NULL, title = "Replicate 1.1")

mse_out_long = reshape::melt(ari_out[,7:9])
mse_plot_1.1 = ggplot(mse_out_long, aes(x = factor(X2, levels = colnames(ari_out)[7:9]), y = value, group = X2, col = factor(X1)))+ #figure 4d
  geom_boxplot(outlier.size = NA)+ggbeeswarm::geom_quasirandom()+
  theme_light()+guides(col = F)+
  scale_x_discrete(labels = c("Spot mean", "Deconvolution (c = 0.003)",
                              "Deconvolution (c = 0.01)"))+
  labs(x = NULL, y = "MSE")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

ari_files = list.files("data-raw/deconv_sim_1.2/")
ari_out = matrix(nrow = 30, ncol = 9)
colnames(ari_out) = names(readRDS(paste0("data-raw/deconv_sim_1.2/", ari_files[1])))
rownames(ari_out) = 1:30
for (i in 1:length(ari_files)){
  ari_out[i,] = readRDS(paste0("data-raw/deconv_sim_1.2/", ari_files[i]))  
}
ari_out_long = reshape::melt(ari_out[,c(1,3:6)]) %>% 
  filter(X2 %in% c("deconv_t_c0.003", "deconv_t_c0.01"))
modes = numeric(nrow(Y1.2))
for(i in 1:nrow(Y1.2)){
  modes[i] =  Mode(deconv1.2col[i + (0:8)*nrow(Y1.2)])
}

yint2 = adjustedRandIndex(rep(modes,9), deconv1.2col)

ari_plot_1.2 = ggplot(ari_out_long, aes(x = factor(X2), y = value, group = factor(X2), col = factor(X1)))+ #figure 4d
  geom_boxplot(outlier.size = 0, outlier.stroke = 0)+ggbeeswarm::geom_quasirandom()+
  geom_hline(color = "red", aes(yintercept = yint2, linetype = "Optimal spot-level clustering"), show.legend = T) + 
  theme_light()+guides(col = F)+
  scale_x_discrete(labels = c("Deconvolution (c = 0.003)",
                              "Deconvolution (c = 0.01)"))+
  scale_linetype_manual(values = 2) + 
  labs(x = NULL, y = "ARI", linetype = NULL, title = "Replicate 1.2")

ari_plot_1.1+ari_plot_1.2+plot_layout(guides = "collect") & theme(legend.position = 'bottom')
mse_out_long = reshape::melt(ari_out[,7:9])
mse_plot_1.2 = ggplot(mse_out_long, aes(x = factor(X2, levels = colnames(ari_out)[7:9]), y = value, group = X2, col = factor(X1)))+ #figure 4d
  geom_boxplot(outlier.size = NA)+ggbeeswarm::geom_quasirandom()+
  theme_light()+guides(col = F)+
  scale_x_discrete(labels = c("Spot mean", "Deconvolution (c = 0.003)",
                              "Deconvolution (c = 0.01)"))+
  labs(x = NULL, y = "MSE")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

#Integration
spotlight = read.csv("data-raw/cell_type_proportions.no_unclassified.csv")
colnames(spotlight)
plot_melanoma = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = spotlight$Melanoma)) +
  geom_point(size = 5.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Proportion", title = "Melanoma") +
  theme_void() + 
  coord_fixed()
plot_macro = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = spotlight$Macrophage)) +
  geom_point(size = 5.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Proportion", title = "Macrophage") +
  theme_void() + 
  coord_fixed()
plot_b = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = spotlight$B.cell)) +
  geom_point(size = 5.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Proportion", title = "B cell") +
  theme_void() + 
  coord_fixed()
plot_t = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = spotlight$T.cell)) +
  geom_point(size = 5.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Proportion", title = "T cell") +
  theme_void() + 
  coord_fixed()
plot_caf = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = spotlight$CAF)) +
  geom_point(size = 5.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Proportion", title = "CAF") +
  theme_void() + 
  coord_fixed()
plot_endothelial = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = spotlight$Endothelial)) +
  geom_point(size = 5.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Proportion", title = "Endothelial") +
  theme_void() + 
  coord_fixed()
(plot_melanoma | plot_caf|plot_endothelial)/(plot_b |plot_t |plot_macro) + plot_layout(guides = "collect") & 
  scale_fill_gradient(low="#ffffff", high="#e50000", breaks=c(0, 0.25, 0.5, 0.75, 1), 
                      oob = scales::squish, labels = c("0", "0.25", ">0.5", ">0.75", "1"),
                      limits = c(0, 0.5)) 


pred_proportion = matrix(nrow = 7, ncol = nrow(deconv1.2Y))
rownames(pred_proportion) = colnames(spotlight)[2:8]
colnames(deconv1.2Y) = paste0("PC", 1:10)
for (i in 2:8){
  reg = lm(spotlight[,i]~ ., data = data.frame(Y1.2))
  pred_proportion[i-1,] = predict(reg, newdata = data.frame(deconv1.2Y))
}
pred_proportion2 = pred_proportion
pred_proportion2[pred_proportion2<0] = 0
pred_proportion2[pred_proportion2>1] = 1
pred_proportion_normalized = t(t(pred_proportion2)/colSums(pred_proportion2))

plot_melanoma = ggplot(data.frame(deconv1.2$obj$positions), aes(x = x, y = y, fill = pred_proportion_normalized["Melanoma",])) +
  geom_point(size = 1.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Proportion", title = "Melanoma") +
  theme_void() + 
  coord_fixed()
plot_macro = ggplot(data.frame(deconv1.2$obj$positions), aes(x = x, y = y, fill = pred_proportion_normalized["Macrophage",])) +
  geom_point(size = 1.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Proportion", title = "Macrophage") +
  theme_void() + 
  coord_fixed()
plot_b = ggplot(data.frame(deconv1.2$obj$positions), aes(x = x, y = y, fill = pred_proportion_normalized["B.cell",])) +
  geom_point(size = 1.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Proportion", title = "B cell") +
  theme_void() + 
  coord_fixed()
plot_t = ggplot(data.frame(deconv1.2$obj$positions), aes(x = x, y = y, fill = pred_proportion_normalized["T.cell",])) +
  geom_point(size = 1.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Proportion", title = "T cell") +
  theme_void() + 
  coord_fixed()
plot_caf = ggplot(data.frame(deconv1.2$obj$positions), aes(x = x, y = y, fill = pred_proportion_normalized["CAF",])) +
  geom_point(size = 1.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Proportion", title = "CAF") +
  theme_void() + 
  coord_fixed()
plot_endothelial = ggplot(data.frame(deconv1.2$obj$positions), aes(x = x, y = y, fill = pred_proportion_normalized["Endothelial",])) +
  geom_point(size = 1.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Proportion", title = "Endothelial") +
  theme_void() + 
  coord_fixed()
plot_deconv = (plot_melanoma | plot_caf|plot_endothelial)/(plot_b |plot_t |plot_macro) + plot_layout(guides = "collect") & 
  scale_fill_gradient(low="#ffffff", high="#e50000", breaks=c(0, 0.25, 0.5, 0.75, 1), 
                      oob = scales::squish, labels = c("0", "0.25", "0.5", "0.75", "1"),
                      limits = c(0,0.75)) 
ggplot(data.frame(deconv1.2$obj$positions), aes(x = x, y = y, fill = rownames(pred_proportion_normalized)[apply(pred_proportion_normalized,2, which.max)],
                                                alpha = apply(pred_proportion_normalized,2, max))) +
  geom_point(size = 2, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Proportion", title = "Endothelial") +
  guides(alpha = F)+
  scale_fill_manual(values = c("yellow", "red", "purple", "blue", "green", "orange"))+
  scale_alpha_continuous(range = c(0,1))+
  theme_void() + 
  coord_fixed()
 