library(tidyverse)
library(SingleCellExperiment)
library(patchwork)
library(scater)
library(mclust)
library(scran)
library(teigen)
library(xgboost)
source("script.R")

melanoma1.1 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel1_rep1.RDS")
melanoma1.2 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel1_rep2.RDS")

# dec1 <- modelGeneVarByPoisson(melanoma1.1)
# top1 <- getTopHVGs(dec1, prop=0.1)
# dec2 <- modelGeneVarByPoisson(melanoma1.2)
# top2 <- getTopHVGs(dec2, prop=0.1)

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


# markers = c("CD2", "CD3D", "CD3E", "CD3G", "CD19", "CD79A", "CD79B", "BLK",
#             "CD163", "CD14", "CSF1R", "PECAM1", "VWF", "CDH5", "FAP", "THY1", "DCN",
#             "MIA", "TYR", "PMEL", "MLANA", "PRAME")
# HVG = union(union(top1, top2), markers)

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
# Cluster
clust1.1 = cluster(Y= Y1.1, positions = positions1.1, q = 4, model = "t", init = mclust1.1, nrep = 10000, gamma = 2, dist = 1)
clust1.2 = cluster(Y= Y1.2, positions = positions1.2, q = 4, model = "t", init = mclust1.2, nrep = 10000, gamma = 2, dist = 1)
saveRDS(clust1.1, "data-raw/clust1.1_q4_figure4b.RDS")
saveRDS(clust1.2, "data-raw/clust1.2_q4_figure4b.RDS")
clust1.1 = readRDS("data-raw/clust1.1_q4_figure4b.RDS")
clust1.2 = readRDS("data-raw/clust1.2_q4_figure4b.RDS")
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

#Deconvolution
deconv1.1 = readRDS("data-raw/deconv1.1_6-17_c0.003_small.RDS") 
deconv1.2 = readRDS("data-raw/deconv1.2_6-17_c0.003_small.RDS")
deconv1.1 = deconvolve(Y= Y1.1, positions = positions1.1, q = 4, init = clust1.1$labels,  model = "t",
                       nrep = 5000, gamma = 2, xdist = 1, ydist = 1, platform = "ST", c = 0.003, jitter_scale = 1.6)
deconv1.2 = deconvolve(Y= Y1.2, positions = positions1.2, q = 4, init = clust1.2$labels,  model = "t",
                       nrep = 500, gamma = 2, xdist = 1, ydist = 1, platform = "ST", c = 0.003, jitter_scale = 8)
deconv1.2col = apply(deconv1.2$obj$z[seq(500,2000,1),], 2, Mode)
deconv1.2_alpha = pmax(colMeans(deconv1.2$obj$z[seq(500,2000,1),]==1),
                       colMeans(deconv1.2$obj$z[seq(500,2000,1),]==2),
                       colMeans(deconv1.2$obj$z[seq(500,2000,1),]==3),
                       colMeans(deconv1.2$obj$z[seq(500,2000,1),]==4))
deconv1.2Y = Reduce(`+`, deconv1.2$obj$Y[-(1:100)])/length(c(1:length(deconv1.2$obj$Y))[-(1:100)])
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
plot(colQuantiles(xgboost_hyper, probs = 0.05))
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

plot_deconv1.2 = ggplot(as.data.frame(deconv1.2$obj$positions),
       aes(x = x, y = y, alpha = deconv1.2_alpha, fill = factor(deconv1.2col))) +
  geom_point(size = 2, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  guides(alpha = F, fill = F) +
  scale_fill_manual(values = c("red", "blue", "purple", "yellow"))+
  theme_void()+ coord_fixed()

placeholder = ggplot(as.data.frame(deconv1.2$obj$positions),
                     aes(x = y, y = x, alpha = deconv1.2$alpha, fill = factor(deconv1.2$cols)))+
  geom_blank()+theme_void()
(placeholder | plot_clust1.2 |plot_deconv1.2) + plot_annotation(tag_levels = "A")

#Clustering Simulation
{
clust1.1mu = matrix(colMeans(clust1.1$mu[-(1:1000),]), byrow = T, ncol = 10)
clust1.1lambda = Reduce(`+`, clust1.1$lambda[seq(1000,10000, 1)])/length(seq(1000,10000, 1))

clust1.2mu = matrix(colMeans(clust1.2$mu[-(1:1000),]), byrow = T, ncol = 10)
clust1.2lambda = Reduce(`+`, clust1.2$lambda[seq(1000,10000, 1)])/length(seq(1000,10000, 1))
simY = matrix(nrow = nrow(positions1.1), ncol = 10)

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

# Deconvolution Simulation
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
 
simY = matrix(nrow = nrow(deconv1.1positions), ncol = 10)
set.seed(100)
for(i in 1:nrow(deconv1.1positions)){
  simY[i,] = mvnfast::rmvt(1, mu = deconv1.1mu[deconv1.1col[i],], sigma = solve(deconv1.1lambda), df = 4)
}

data_sim = data.frame(simY)
data_sim$j = rep(1:nrow(Y1.1), 9)
data_sim %>%
  group_by(j) %>%
  summarise_all(mean)->
  data_sim_mean

agg1 = ggplot(data.frame(positions1.2), aes(x, y, fill = data_sim_mean$X1)) +
  geom_point(size = 7, pch = 22)+
  labs(fill = "PC1", alpha = "Proportion", x = NULL, y = NULL) +
  theme_void() + coord_fixed()
agg2 = ggplot(data.frame(positions1.2), aes(x, y, fill = data_sim_mean$X2)) +
  geom_point(size = 7, pch = 22)+
  labs(fill = "PC2", alpha = "Proportion", x = NULL, y = NULL) +
  theme_void() + coord_fixed()
agg3 = ggplot(data.frame(positions1.2), aes(x, y, fill = data_sim_mean$X3)) +
  geom_point(size = 7, pch = 22)+
  labs(fill = "PC3", alpha = "Proportion", x = NULL, y = NULL) +
  theme_void() + coord_fixed()
agg1+agg2+agg3

set.seed(101)
km_sim = kmeans(data_sim_mean[,-1], centers = 3)$cluster
mclust::adjustedRandIndex(deconv1.2col, rep(km_sim,9))
set.seed(102)
mclust_sim = Mclust(data_sim_mean[,-1], G = 3, modelNames = c("EEE"))$classification
mclust::adjustedRandIndex(deconv1.2col, rep(mclust_sim,9))

clust_sim = cluster(Y= data_sim_mean[,-1], positions = as.matrix(positions1.2), q = 3, init = mclust_sim, 
                    nrep = 10000, gamma = 2, dist = 1)
mclust::adjustedRandIndex(deconv1.2col, rep(clust_sim$labels,9))

km_sim_out = ggplot(data.frame(positions1.2), aes(x, y, fill = factor(km_sim))) +
  geom_point(size = 7, pch = 22)+
  guides(fill = F)+ 
  scale_color_viridis_d(option = "A") + 
  labs(fill = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  scale_fill_manual(values = c("purple", "red", "yellow", "grey"))+
  theme_void() + coord_fixed()
clust_sim_out = ggplot(data.frame(positions1.2), aes(x, y, fill = factor(clust_sim$labels))) +
  geom_point(size = 7, pch = 22)+
  guides(fill = F)+ 
  scale_color_viridis_d(option = "A") + 
  labs(fill = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  scale_fill_manual(values = c("yellow", "red", "purple", "grey"))+
  theme_void() + coord_fixed()

deconv_sim = deconvolve(Y= data_sim_mean[,-1], positions = positions1.2, q = 3, init = clust_sim$labels, 
                        nrep = 1000, gamma = 2, xdist = 1, ydist = 1, platform = "ST", c = 0.003)
deconv_simcol = apply(deconv_sim$z[seq(10000,100000,10),], 2, Mode)
deconv_sim_alpha = pmax(colMeans(deconv_sim$z[seq(10000,100000,10),]==1),
                        colMeans(deconv_sim$z[seq(10000,100000,10),]==2),
                        colMeans(deconv_sim$z[seq(10000,100000,10),]==3))
deconv_sim = readRDS("data-raw/deconvsim_g2_c0.01.RDS")
deconv_sim_col_c0.01 = deconv_sim$cols
deconv_sim_alpha_c0.01 = deconv_sim$alpha
rm(deconv_sim)
deconv_sim = readRDS("data-raw/deconvsim_g2_c0.003.RDS")
deconv_sim_col_c0.003 = deconv_sim$cols
deconv_sim_alpha_c0.003 = deconv_sim$alpha
rm(deconv_sim)
deconv_sim = readRDS("data-raw/deconvsim_g2_c0.001.RDS")
deconv_sim_col_c0.001 = deconv_sim$cols
deconv_sim_alpha_c0.001 = deconv_sim$alpha
rm(deconv_sim)

mclust::adjustedRandIndex(deconv1.2col, deconv_sim_col_c0.01)
mclust::adjustedRandIndex(deconv1.2col, deconv_sim_col_c0.003)
mclust::adjustedRandIndex(deconv1.2col, deconv_sim_col_c0.001)

plot_deconv_sim_0.01 = ggplot(as.data.frame(deconv_sim$obj$positions),
                        aes(x = x, y = y, alpha = deconv_sim_alpha_c0.01, fill = factor(deconv_sim_col_c0.01))) +
  geom_point(size = 2.1, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  guides(alpha = F, fill = F) +
  scale_fill_manual(values = c("yellow", "red", "purple", "grey"))+
  theme_void()+ coord_fixed()


plot_deconv_sim_0.003 = ggplot(as.data.frame(deconv_sim$obj$positions),
                              aes(x = x, y = y, alpha = deconv_sim_alpha_c0.003, fill = factor(deconv_sim_col_c0.003))) +
  geom_point(size = 2.1, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  guides(alpha = F, fill = F) +
  scale_fill_manual(values = c("yellow", "red", "purple", "grey"))+
  theme_void()+ coord_fixed()


plot_deconv_sim_0.001 = ggplot(as.data.frame(deconv_sim$obj$positions),
                              aes(x = x, y = y, alpha = deconv_sim_alpha_c0.001, fill = factor(deconv_sim_col_c0.001))) +
  geom_point(size = 2.1, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  guides(alpha = F, fill = F) +
  scale_fill_manual(values = c("yellow", "red", "purple", "grey"))+
  theme_void()+ coord_fixed()
plot_deconv_sim_0.01 + plot_deconv_sim_0.003+ plot_deconv_sim_0.001

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

(placeholder | plot_clust1.2 )/(plot_deconv1.2 | ari_plot) + plot_annotation(tag_levels = "A")



pblank = ggplot(ari_out_long)+geom_blank()+theme_void()
(pblank+pblank+pblank)/ari_plot+plot_annotation(tag_levels = "A")
#Deconvolved PCs => lognormcounts
dec <- modelGeneVarByPoisson(melanoma1.2)
top <- getTopHVGs(dec, prop=0.1)
logNormCounts(melanoma1.1)["IGLL5",]
colnames(deconv1.2Y) = paste0("PC",1:10)
deconv_data = data.frame(deconv1.2Y)
deconv_expression = matrix(nrow = length(top), ncol = nrow(deconv_data))
rownames(deconv_expression) = top
Y1.2data = data.frame(Y1.2)
rsquared = numeric(length(top))
names(rsquared) = top
for (gene in top){
  train = lm(logcounts(melanoma1.2)[gene,]~. , data = Y1.2data)
  rsquared[gene] = summary(train)$r.squared
  deconv_expression[gene,] = predict(train, newdata = deconv_data)
}
deconv_expression2 = predictExpression(melanoma1.2, newdata = deconv1.2Y, components = 10)

saveRDS(deconv_expression, "data-raw/deconv_expression.RDS")
deconv_expression = readRDS("data-raw/deconv_expression.RDS")

cd6_plot1 = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = logcounts(melanoma1.2)["CD6",])) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Expression") +
  scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed() +
  guides(fill=F)
cd6_plot2 = ggplot(as.data.frame(deconv1.2$obj$positions),
       aes(x = x, y = y, fill = deconv_expression$expression["CD6",])) +
  geom_point(size = 2, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_viridis(option="A")+
  guides(alpha = F, fill = F) +
  theme_void()+ coord_fixed()
(cd2_plot1|cd2_plot2)/
  (cd6_plot1|cd6_plot2)/
  (cd14_plot1|cd14_plot2)/
  (cd19_plot1|cd19_plot2)/
  (cd20_plot1|cd20_plot2)/
  (cd45_plot1|cd45_plot2)
  
  
  

#logliks for cluster number 
logliks = readRDS("data-raw/logliks_2.1.RDS")
logliks
library(tidyverse)
reshape::melt(logliks) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = X1, y = value, col = X2))+
  geom_point()+geom_line()+
  labs(x = "Number of clusters", y = "Pseudo-log-likelihood")

logliks2 = matrix(nrow = 6, ncol = 2)
colnames(logliks2) = c("normal_var", "t_var")

for (q in 2:6){
  set.seed(101)
  mclust_q = Mclust(Y1.2, G = q, modelNames = "EEE")$classification
  logliks2[q,1] = mean(cluster(Y= Y1.2, positions = positions1.2, q = q, init = mclust_q, nrep = 1000, gamma = 2, dist = 1, 
                               precision = "variable")$plogLik[100:1000])
  logliks2[q,2] = mean(cluster(Y= Y1.2, positions = positions1.2, q = q, init = mclust_q, nrep = 1000, gamma = 2, dist = 1, alpha = 10, beta = 10,
                               precision = "variable", model = "t")$plogLik[100:1000])
  print(logliks2[q,])
}

plot_pll = reshape::melt(logliks) %>% rbind(reshape::melt(logliks2)) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = X1, y = value, col = X2))+
  geom_point()+geom_line()+
  scale_x_continuous(breaks = 1:15, labels = 1:15) +
  theme_light()+
  labs(x = "Number of clusters", y = "Pseudo-log-likelihood", col = "Model")
logliksv2 = logliks
for (i in 2:15){
  logliksv2[i, ] = logliks[i,] - 1/2 * i*10 * log(nrow(Y1.2))
}
logliks2v2 = logliks2
for (i in 2:6){
  logliks2v2[i, ] = logliks2[i,] - 1/2 * i*(10+55) * log(nrow(Y1.2))
}

plot_bic = reshape::melt(logliksv2) %>% rbind(reshape::melt(logliks2v2)) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = X1, y = value, col = X2))+
  geom_point()+geom_line()+
  scale_x_continuous(breaks = 1:15, labels = 1:15) +
  theme_light()+
  labs(x = "Number of clusters", y = "BIC", col = "Model")
plot_pll+plot_bic+plot_layout(guides = "collect")

#4 clusters
set.seed(101)
km1.2 = kmeans(Y1.2, centers = 5)$cluster
set.seed(100)
mclust1.1 = Mclust(Y1.1, G = 4, modelNames = "EEE")$classification
mclust1.2 = Mclust(Y1.2, G = 5, modelNames = "EEE")$classification
mclust2.1 = Mclust(Y2.1, G = 6, modelNames = "EEE")$classification
mclust2.2 = Mclust(Y2.2, G = 5, modelNames = "EEE")$classification

set.seed(100)
g.jaccard1.1 <- buildSNNGraph(melanoma1.1, use.dimred="PCA", type="jaccard")
clust.walktrap1.1 <- igraph::cut_at(igraph::cluster_walktrap(g.jaccard1.1), n = 4)
g.jaccard1.2 <- buildSNNGraph(melanoma1.2, use.dimred="PCA", type="jaccard")
clust.walktrap1.2 <- igraph::cut_at(igraph::cluster_walktrap(g.jaccard1.2), n = 4)
ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = factor(km1.2))) +
  geom_point(size = 4.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #scale_fill_manual(values = c("purple", "red", "yellow", "black", "grey"))+
  guides(fill = F) +
  theme_void()+ coord_fixed()

clust1.1 = cluster(Y= Y1.1, positions = positions1.1, q = 4, init = clust.walktrap1.1, nrep = 10000, gamma = 2, dist = 1, model = "t")
clust1.2 = cluster(Y= Y1.2, positions = positions1.2, q = 4, init = clust.walktrap1.2, nrep = 10000, gamma = 2, dist = 1, model = "t")
clust2.1 = cluster(Y= Y2.1, positions = positions2.1, q = 6, init = mclust2.1, nrep = 10000, gamma = 2, dist = 1, model = "t")
clust2.2 = cluster(Y= Y2.2, positions = positions2.2, q = 5, init = mclust2.2, nrep = 10000, gamma = 2, dist = 1, model = "t")
clust1.2v2 = cluster(Y= Y1.2, positions = positions1.2, q = 4, init = km1.2_4, 
                     nrep = 1000, gamma = 2, dist = 1, model = "t")
clust1.2q5 = cluster(Y= Y1.2, positions = positions1.2, q = 5, init = mclust1.2, 
                     nrep = 1000, gamma = 2, dist = 1, model = "t")
clust5_1.1 = ggplot(data.frame(positions1.1), aes(x = x, y = y, fill = factor(clust1.1$labels))) +
  geom_point(pch =15, color = "grey", size = 10)+
  geom_point(size = 4.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("red", "violet", "yellow", "purple"))+
  guides(fill = F) +
  theme_void()+ coord_fixed()
clust5_1.2 = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = factor(clust1.2$labels))) +
  geom_point(pch =15, color = "grey", size = 10)+
  geom_point(size = 4.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("purple", "red", "violet", "yellow"))+
  guides(fill = F) +
  theme_void()+ coord_fixed()
clust5_2.1 = ggplot(data.frame(positions2.1), aes(x = x, y = y, fill = factor(clust2.1$labels))) +
  geom_point(size = 4.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #scale_fill_manual(values = c("purple", "red", "yellow", "black"))+
  guides(fill = F) +
  theme_void()+ coord_fixed()
clust5_2.2 = ggplot(data.frame(positions2.2), aes(x = x, y = y, fill = factor(clust2.2$labels))) +
  geom_point(size = 4.5, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #scale_fill_manual(values = c("purple", "red", "yellow", "black"))+
  guides(fill = F) +
  theme_void()+ coord_fixed()
(clust5_1.1 | clust5_1.2) / (clust5_2.1 | clust5_2.2) + plot_annotation(tag_levels = "A")

deconv1.2 = deconvolve(Y= Y1.2, positions = positions1.2, q = 4, init = km1.2_4, nrep = 20000, gamma = 2, xdist = 1, ydist = 1, platform = "ST", c = 0.003)
deconv1.2col = apply(deconv1.2$z[seq(1000,20000,10),] ,2, Mode)
deconv1.2_alpha = pmax(colMeans(deconv1.2$z[seq(1000,20000,10),]==1),
                       colMeans(deconv1.2$z[seq(1000,20000,10),]==2),
                       colMeans(deconv1.2$z[seq(1000,20000,10),]==3),
                       colMeans(deconv1.2$z[seq(1000,20000,10),]==4))
ggplot(as.data.frame(deconv1.2$positions),
       aes(x = x, y = y, alpha = deconv1.2_alpha, fill = factor(deconv1.2col))) +
  geom_point(size = 2, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  guides(alpha = F, fill = F) +
  scale_fill_manual(values = c("purple", "red", "yellow", "black"))+
  theme_void()+ coord_fixed()
table(deconv1.2col)
saveRDS(list(clust1.1 = clust1.1$labels, 
     clust1.2 = clust1.2$labels,
     clust2.1 = clust2.1$labels,
     clust2.2 = clust2.2$labels,
     deconv1.2_3 = deconv1.2$cols,
     deconv1.2_4 = deconv1.2col), "data-raw/labels.RDS")

#xgboost
xgboost_expression = matrix(ncol = nrow(deconv1.2Y), nrow = length(top))
rownames(xgboost_expression) = top
for (gene in top){
  train = xgboost(data = Y1.2, label = logcounts(melanoma1.2)[gene,], max_depth = 2,
                  eta = 0.03, nthread = 2, nrounds = 200, objective = "reg:squarederror", verbose = F)
  xgboost_expression[gene, ] = predict(train, deconv1.2Y)
}
set.seed(1)
traindata = sample(1:293, 200)

xgboost_hyper = matrix(nrow = length(top), ncol = 500)
rownames(xgboost_hyper) = top
for(gene in top){
  Y1.2train = xgb.DMatrix(data = Y1.2[traindata,], label = logcounts(melanoma1.2)[gene,traindata])
  Y1.2test = xgb.DMatrix(data = Y1.2[-traindata,], label = logcounts(melanoma1.2)[gene,-traindata])
  watchlist = list(train = Y1.2train, test = Y1.2test)
  train = xgb.train(data = Y1.2train, max_depth = 2, watchlist = watchlist,
                    eta = 0.03, nthread = 2, nrounds = 500, objective = "reg:squarederror", verbose = F)
  xgboost_hyper[gene,] = train$evaluation_log$test_rmse
}
colnames(xgboost_hyper ) = 1:500
plot(colQuantiles(xgboost_hyper, probs = 0.05))
reshape::melt(xgboost_hyper) %>%
  as_tibble() %>%
  filter(X2 %in% seq(1, 500, 10)) %>%
  ggplot(aes(x = X2, y = value, group = X2))+
  geom_boxplot(width = 5)

Y1.2train = xgb.DMatrix(data = Y1.2[traindata,], label = logcounts(melanoma1.2)["CD3G",traindata])
Y1.2test = xgb.DMatrix(data = Y1.2[-traindata,], label = logcounts(melanoma1.2)["CD3G",-traindata])
watchlist = list(train = Y1.2train, test = Y1.2test)
train = xgb.train(data = Y1.2train, max_depth = 2, watchlist = watchlist,
        eta = 0.03, nthread = 2, nrounds = 200, objective = "reg:squarederror")
test = predict(train, deconv1.2Y)
test2 = predictExpression(melanoma1.2, deconv1.2Y, genes = "CD14")
plot(x = test2$expression[1,], y = test, xlab = "Linear regression", ylab = "xgboost")
