library(tidyverse)
library(SingleCellExperiment)
library(viridis)
library(patchwork)
# melanoma1.1 = readRDS("ST_mel1_rep1.RDS")
# melanoma1.2 = readRDS("ST_mel1_rep2.RDS")
# melanoma2.1 = readRDS("ST_mel2_rep1.RDS")
# melanoma2.2 = readRDS("ST_mel2_rep2.RDS")
melanoma1.1 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel1_rep1.RDS")
melanoma1.2 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel1_rep2.RDS")
melanoma2.1 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel2_rep1.RDS")
melanoma2.2 = readRDS("data-raw/ST-Melanoma-Datasets_1/ST_mel2_rep2.RDS")
library(scater)
library(mclust)
source("script.R")
set.seed(100)
melanoma1.1 = runPCA(melanoma1.1)
melanoma1.2 = runPCA(melanoma1.2)
melanoma2.1 = runPCA(melanoma2.1)
melanoma2.2 = runPCA(melanoma2.2)
positions1.1 = cbind(x = - melanoma1.1$col, y = - melanoma1.1$row)
positions1.2 = cbind(x = melanoma1.2$col, y = melanoma1.2$row)
positions2.1 = cbind(y = melanoma2.1$col, x = - melanoma2.1$row)
positions2.2 = cbind(y = - melanoma2.2$col, x = melanoma2.2$row)
Y1.1 = reducedDim(melanoma1.1)[,1:10]
Y1.2 = reducedDim(melanoma1.2)[,1:10]
Y2.1 = reducedDim(melanoma2.1)[,1:10]
Y2.2 = reducedDim(melanoma2.2)[,1:10]

plot1.2 = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = librarySizeFactors(melanoma1.2))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Library Size") +
  scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed()
plot2.1 = ggplot(data.frame(positions2.1), aes(x = x, y = y, fill = librarySizeFactors(melanoma2.1))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Library Size") +
  scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed()

set.seed(100)
km1.1 = kmeans(Y1.1, 4)$cluster
km1.2 = kmeans(Y1.2, 4)$cluster
km2.1 = kmeans(Y2.1, 4)$cluster
km2.2 = kmeans(Y2.2, 4)$cluster
plot_km1.2 = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = factor(km1.2))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("purple", "yellow", "red", "grey"))+
  theme_void()+ coord_fixed()
plot_km2.1 = ggplot(data.frame(positions2.1), aes(x = x, y = y, fill = factor(km2.1))) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #scale_fill_manual(values = c("orange", "black", "grey", "red", "green"))+
  theme_void()+ coord_fixed()

# Cluster
clust1.2 = cluster(Y= Y1.2, positions = positions1.2, q = 4, model = "t", init = km1.2, nrep = 10000, gamma = 2, dist = 1)
saveRDS(clust1.2, "data-raw/clust1.2_q4_figure4b.RDS")
clust1.2 = readRDS("data-raw/clust1.2_q4_figure4b.RDS")
clust1.2_alpha = pmax(colMeans(clust1.2$z[seq(1000,10000,1),]==1),
                      colMeans(clust1.2$z[seq(1000,10000,1),]==2),
                      colMeans(clust1.2$z[seq(1000,10000,1),]==3),
                      colMeans(clust1.2$z[seq(1000,10000,1),]==4))
# clust1.2 = readRDS("data-raw/clust1.2_figure3b.RDS") #figure 4b

clust2.1 = cluster(Y= Y2.1, positions = positions2.1, q = 3, init = km2.1, nrep = 10000, gamma = 2, dist = 1)

plot_clust1.2 = ggplot(data.frame(positions1.2), aes(x = x, y = y, alpha = clust1.2_alpha, fill = factor(clust1.2$labels))) +
  geom_point(size = 6, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("red", "blue", "purple", "yellow"))+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+
  guides(alpha = F, fill = F) +
  theme_void()+ coord_fixed()
plot_clust2.1 = ggplot(data.frame(positions2.1), aes(x = x, y = y, fill = factor(clust2.1$labels))) +
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_manual(values = c("purple", "yellow", "black", "grey"))+
  theme_void()+ coord_fixed()

#Deconvolution
#deconv1.2 = readRDS("data-raw/deconv1.2_g2_c0.003_withY.RDS") #FIGURE 4C
deconv1.2 = readRDS("data-raw/deconv1.2_normal_q4_g2_c0.003.RDS") #FIGURE 4C
deconv1.2 = readRDS("data-raw/deconv1.2_t_q4_g2_c0.003.RDS") #FIGURE 4C
deconv1.2 = readRDS("data-raw/deconv1.2_t_q4_g2_c0.003_200k_small.RDS") #FIGURE 4C
# deconv1.2 = deconvolve(Y= Y1.2, positions = positions1.2, q = 4, init = km1.2,  model = "t", 
#                        nrep = 100, gamma = 2, xdist = 1, ydist = 1, platform = "ST", c = 0.003)
deconv1.2col = apply(deconv1.2$obj$z[seq(500,2000,1),], 2, Mode)
# deconv1.2_alpha = pmax(colMeans(deconv1.2$z[seq(10000,100000,10),]==1),
#                        colMeans(deconv1.2$z[seq(10000,100000,10),]==2),
#                        colMeans(deconv1.2$z[seq(10000,100000,10),]==3),
#                        colMeans(deconv1.2$z[seq(10000,100000,10),]==4))

# saveRDS(list(obj = deconv1.2, cols = deconv1.2col, alpha = deconv1.2_alpha), "deconv1.2_g1_c0.01.RDS")
plot_deconv1.2 = ggplot(as.data.frame(deconv1.2$obj$positions),
       aes(x = x, y = y, alpha = deconv1.2$alpha, fill = factor(deconv1.2col))) +
  geom_point(size = 2, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  guides(alpha = F, fill = F) +
  scale_fill_manual(values = c("red", "blue", "purple", "yellow"))+
  theme_void()+ coord_fixed()


deconv2.1 = readRDS("data-raw/deconv2.1_g2_c0.003.RDS")
# deconv2.1 = deconvolve(Y= Y2.1, positions = positions2.1, q = 3, init = km2.1, nrep = 100000, gamma = 1, xdist = 1, ydist = 1, platform = "ST", c = 0.01)
# deconv2.1col = apply(deconv2.1$z[seq(10000,100000,10),], 2, Mode)
# deconv2.1_alpha = pmax(colMeans(deconv2.1$z[seq(10000,100000,10),]==1),
#                        colMeans(deconv2.1$z[seq(10000,100000,10),]==2),
#                        colMeans(deconv2.1$z[seq(10000,100000,10),]==3))
# saveRDS(list(obj = deconv2.1, cols = deconv2.1col, alpha = deconv2.1_alpha), "deconv2.1_g1_c0.01.RDS")

plot_deconv2.1 = ggplot(as.data.frame(deconv2.1$obj$positions),
       aes(x = y, y = x, alpha = deconv2.1$alpha, fill = factor(deconv2.1$cols))) +
  geom_point(size = 1, pch = 22)+
  scale_alpha_continuous(limits = c(0,1), range = c(0,1))+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  guides(alpha = F) +
  scale_fill_manual(values = c("purple", "yellow", "purple", "grey"))+
  theme_void()+ coord_fixed()

placeholder = ggplot(as.data.frame(deconv1.2$obj$positions),
                     aes(x = y, y = x, alpha = deconv1.2$alpha, fill = factor(deconv1.2$cols)))+
  geom_blank()+theme_void()

(placeholder | plot_clust1.2 |plot_deconv1.2) + plot_annotation(tag_levels = "A")

#Clustering Simulation
clust1.2mu = matrix(colMeans(clust1.2$mu[-(1:1000),]), byrow = T, ncol = 10)
clust1.2lambda = Reduce(`+`, clust1.2$lambda[seq(1000,10000, 1)])/length(seq(1000,10000, 1))
simY = matrix(nrow = nrow(positions1.2), ncol = 10)
set.seed(100)
for(i in 1:nrow(positions1.2)){
  simY[i,] = mvnfast::rmvt(1, mu = clust1.2mu[clust1.2$labels[i],], sigma = solve(clust1.2lambda), df = 4)
}
ggplot(data.frame(positions1.2), aes(x, y, fill = simY[,1])) +
  geom_point(size = 7, pch = 22)+
  labs(fill = "PC1", alpha = "Proportion", x = NULL, y = NULL) +
  theme_void() + coord_fixed()

set.seed(101)
km_sim = kmeans(simY, centers = 4)$cluster
mclust::adjustedRandIndex(clust1.2$labels, km_sim)
set.seed(102)
mclust_sim = Mclust(simY, G = 4, modelNames = c("EEE"))$classification
mclust::adjustedRandIndex(clust1.2$labels, mclust_sim)

clust_sim = cluster(Y= simY, positions = positions1.2, q = 4, init = mclust_sim, 
                    nrep = 10000, gamma = 2, dist = 1, model = "normal")
mclust::adjustedRandIndex(clust1.2$labels, clust_sim$labels)

clust_sim_t = cluster(Y= simY, positions = positions1.2, q = 4, init = mclust_sim, 
                    nrep = 10000, gamma = 2, dist = 1, model = "t")
mclust::adjustedRandIndex(clust1.2$labels, clust_sim_t$labels)

ari_files = list.files("data-raw/melanoma_clust_sim/")
ari_out = matrix(nrow = 30, ncol = 5)
colnames(ari_out) = names(readRDS(paste0("data-raw/melanoma_clust_sim/", ari_files[1])))
rownames(ari_out) = 1:30
for (i in 1:length(ari_files)){
  ari_out[i,] = readRDS(paste0("data-raw/melanoma_clust_sim/", ari_files[i]))  
}
ari_out_long = reshape::melt(ari_out)
ari_plot = ggplot(ari_out_long, aes(x = factor(X2, levels = colnames(ari_out)), y = value, group = X2, col = factor(X1)))+ #figure 4d
  geom_boxplot()+geom_point()+
  theme_light()+guides(col = F)+
  scale_x_discrete(labels = c("k-means", "mclust", "Spatial clustering (normal)", "Spatial clustering (t)", "SC (t), truth init"))+
  labs(x = NULL, y = "ARI")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

# Deconvolution Simulation
deconv1.2mu = matrix(colMeans(deconv1.2$obj$mu[50000:200000,]), byrow = T, ncol = 10)
deconv1.2lambda = Reduce(`+`, deconv1.2$obj$lambda[50000:200000])/length(50000:200000)
deconv1.2Y = Reduce(`+`, deconv1.2$obj$Y[-(1:500)])/length(c(1:length(deconv1.2$obj$Y))[-(1:500)])
deconv_positions = deconv1.2$obj$positions
deconv1.2col = apply(deconv1.2$obj$z[seq(500,2000,1),], 2, Mode)
 
simY = matrix(nrow = nrow(deconv_positions), ncol = 10)
set.seed(100)
for(i in 1:nrow(deconv_positions)){
  simY[i,] = mvnfast::rmvt(1, mu = deconv1.2mu[deconv1.2col[i],], sigma = solve(deconv1.2lambda), df = 4)
}

data_sim = data.frame(simY)
data_sim$j = rep(1:nrow(Y1.2), 9)
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

ari_files = list.files("data-raw/melanoma_sim_ari/")
ari_out = matrix(nrow = 10, ncol = 6)
colnames(ari_out) = names(readRDS(paste0("data-raw/melanoma_sim_ari/", ari_files[i])))
rownames(ari_out) = 1:10
for (i in 1:length(ari_files)){
  ari_out[i,] = readRDS(paste0("data-raw/melanoma_sim_ari/", ari_files[i]))  
}
ari_out_long = reshape::melt(ari_out) %>% filter(! X2 %in% c("deconv_c0.001", "deconv_c0.01")) 
ari_plot = ggplot(ari_out_long, aes(x = factor(X2, levels = colnames(ari_out)), y = value, group = X2, col = factor(X1)))+ #figure 4d
  geom_boxplot()+geom_point()+theme_light()+guides(col = F)+
  scale_x_discrete(labels = c("k-means", "mclust", "Spatial clustering", "Spatial deconvolution"))+
  labs(x = NULL, y = "ARI")#+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
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
saveRDS(deconv_expression, "data-raw/deconv_expression.RDS")
deconv_expression = readRDS("data-raw/deconv_expression.RDS")

predictExpression = function(sce, newdata, dimred = "PCA", genes = rownames(sce), components = ncol(newdata)){
  actual_data = data.frame(reducedDim(sce, type = dimred))[,1:components]
  newdata = as.data.frame(newdata)
  if (ncol(actual_data) != ncol(newdata)){
    stop("number of components do not match")
  }
  if (! all(colnames(newdata) == colnames(actual_data))){
    warning("colnames of reducedDim and newdata do not match. Setting newdata colnames to match reducedDim")
    colnames(newdata) = colnames(actual_data)
  }
  rsquared = numeric(length(genes))
  names(rsquared) = genes
  deconv_expression = matrix(nrow = length(genes), ncol = nrow(newdata))
  rownames(deconv_expression) = genes
  colnames(deconv_expression) = rownames(newdata)
  for (gene in genes){
    train = lm(logcounts(sce)[gene,]~. , data = actual_data)
    rsquared[gene] = summary(train)$r.squared
    deconv_expression[gene,] = predict(train, newdata = newdata)
  }
  list(expression = deconv_expression, r2 = rsquared)

}
cd14_plot1 = ggplot(data.frame(positions1.2), aes(x = x, y = y, fill = logcounts(melanoma1.2)["CD3D",])) + 
  geom_point(size = 7, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Expression") +
  scale_fill_viridis(option = "A")+
  theme_void()+ coord_fixed() +
  guides(fill=F)
cd14_plot2 = ggplot(as.data.frame(deconv1.2$obj$positions),
       aes(x = x, y = y, fill = deconv_expression["CD3D",])) +
  geom_point(size = 2, pch = 22)+
  labs(x = NULL, y = NULL, fill = "Cluster") +
  scale_fill_viridis(option="A")+
  guides(alpha = F, fill = F) +
  theme_void()+ coord_fixed()
(cd14_plot1|cd14_plot2)

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
