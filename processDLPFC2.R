library(dplyr)
library(ggplot2)
library(reshape)
library(spatialLIBD)
library(scran)
library(patchwork)

#10 PCs
ARIs = list.files("data-raw/DLPFC3/PC10")
out_matrix = matrix(nrow = length(ARIs), ncol = 9)
colnames(out_matrix) = c("Walktrap (markers)", "Walktrap (HVG)", "k-means", "mclust (EEE)", "mclust (VVV)", 
                         "SC (normal/equal)",  "SC (t/equal)",  "SC (normal/variable)",  "SC (t/variable)" )
for (i in 1:length(ARIs)){
  out_matrix[i,] = readRDS(paste0("data-raw/DLPFC3/PC10/", ARIs[i]))
}
rownames(out_matrix) = matrix(unlist(strsplit(ARIs, "[_.]")),ncol = 3, byrow = T)[,2]

sce <- fetch_data(type = 'sce')

ARI_louvain = numeric(length = 12)
names(ARI_louvain) = rownames(out_matrix)
ARI_stLearn_pca = numeric(length=12)
names(ARI_stLearn_pca) = rownames(out_matrix)
ARI_stLearn_umap = numeric(length=12)
names(ARI_stLearn_umap) = rownames(out_matrix)
ARI_Giotto_3 = numeric(length = 12)
names(ARI_Giotto_3) = rownames(out_matrix)
ARI_Giotto_9 = numeric(length = 12)
names(ARI_Giotto_9) = rownames(out_matrix)

for (name in rownames(out_matrix)){
  sce_sample = sce[,sce$sample_name == name]
  q_sample = length(unique(na.omit(sce_sample$layer_guess)))
  g.jaccard <- buildSNNGraph(sce_sample, use.dimred="PCA", type="jaccard")
  clust.louvain <- igraph::cluster_louvain(g.jaccard)$membership
  ARI_louvain[name] = mclust::adjustedRandIndex(sce_sample$layer_guess, clust.louvain)
  
  stLearn_clusters = read.csv(paste0("data-raw/stLearn_HVGs/clusters/", name, ".csv"))
  ARI_stLearn_pca[name] = mclust::adjustedRandIndex(sce_sample$layer_guess, stLearn_clusters$pca_kmeans)
  ARI_stLearn_umap[name] = mclust::adjustedRandIndex(sce_sample$layer_guess, stLearn_clusters$umap_kmeans)
  
  giotto_i = read.csv(paste0("data-raw/Giotto_HMRF_rerun/", name, ".csv"))
  ARI_Giotto_3[name] = mclust::adjustedRandIndex(sce_sample$layer_guess, giotto_i[,2])
  ARI_Giotto_9[name] = mclust::adjustedRandIndex(sce_sample$layer_guess, giotto_i[,3])
}

out_matrix = cbind(out_matrix, Louvain = ARI_louvain)
out_matrix = out_matrix[,c(1:3, 10, 4:9)]
data_out = melt(out_matrix)
colnames(data_out) = c("sample", "method", "ARI")
data_out$method = factor(data_out$method, 
                         levels = colnames(out_matrix))

data_out = melt(out_matrix)
colnames(data_out) = c("sample", "method", "ARI")
data_out$method = factor(data_out$method, 
                         levels = colnames(out_matrix))
p10 = ggplot(data_out, aes(x = method, y = ARI))+
  geom_boxplot() +
  geom_point(aes(color = factor(sample))) +
  labs(x = "Method", color = "Sample")+
  theme_light()+
  coord_cartesian(ylim = c(0,0.6)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

ARIs = list.files("data-raw/DLPFC3/")[1:12]
out_matrix = matrix(nrow = length(ARIs), ncol = 9)
colnames(out_matrix) = c("Walktrap (markers)", "Walktrap (HVG)", "k-means", "mclust", "mclust (VVV)", 
                         "SC (normal)",  "SC (t)",  "SC (normal/variable)",  "SC (t/variable)" )
for (i in 1:length(ARIs)){
  out_matrix[i,] = readRDS(paste0("data-raw/DLPFC3/", ARIs[i]))
}
rownames(out_matrix) = matrix(unlist(strsplit(ARIs, "[_.]")),ncol = 3, byrow = T)[,2]

out_matrix = cbind(out_matrix, Louvain = ARI_louvain, "stLearn" = ARI_stLearn_pca, "stLearn (UMAP)" = ARI_stLearn_umap,
                   "Giotto" = ARI_Giotto_9)
out_matrix = out_matrix[,c(1:3, 10, 4, 11, 13, 6:7)]
data_out = melt(out_matrix)
colnames(data_out) = c("sample", "method", "ARI")
data_out$method = factor(data_out$method, 
                         levels = colnames(out_matrix))
p15 = ggplot(data_out, aes(x = method, y = ARI))+
  geom_boxplot() +
  geom_point(aes(color = factor(sample))) +
  labs(x = "Method", color = "Sample")+
  theme_light()+
  coord_cartesian(ylim = c(0,0.6)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

p10+p15+plot_layout(guides = "collect")

#Giotto


giotto = read.csv("data-raw/HMRF_domains.betas.csv")
sce_sample = sce[,sce$sample_name == "151673"]
giottos = numeric(22)
names(giottos) = colnames(giotto)[2:23]
names(giottos)[1:4] = paste0("HMRF_PCA_Delaunay_k7_b.", c(3,9,18,27))
for(i in 1:22){
  giottos[i] = mclust::adjustedRandIndex(sce_sample$layer_guess, giotto[,1+i])
}
data.frame(name = names(giottos), ARI = giottos) %>%
  tidyr::separate(name, into = c("HMRF", "input", "network", "k", "gamma"), sep = "_")  %>% 
  #mutate(gamma = as.numeric(gamma)) %>%
  #filter(gamma == "b.1") %>%
  ggplot(aes(x = input, y = ARI, col = factor(gamma, levels = c("b.0.01", "b.0.11", "b.0.33", "b.0.67", "b.1", "b.3", "b.9", "b.18", "b.27")), 
             size = factor(gamma)))+
  geom_point()+labs(col = "gamma", size = "gamma")
