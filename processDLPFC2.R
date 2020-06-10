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
for (name in rownames(out_matrix)){
  sce_sample = sce[,sce$sample_name == name]
  q_sample = length(unique(na.omit(sce_sample$layer_guess)))
  g.jaccard <- buildSNNGraph(sce_sample, use.dimred="PCA", type="jaccard")
  clust.louvain <- igraph::cluster_louvain(g.jaccard)$membership
  ARI_louvain[name] = mclust::adjustedRandIndex(sce_sample$layer_guess, clust.louvain)
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
colnames(out_matrix) = c("Walktrap (markers)", "Walktrap (HVG)", "k-means", "mclust (EEE)", "mclust (VVV)", 
                         "SC (normal/equal)",  "SC (t/equal)",  "SC (normal/variable)",  "SC (t/variable)" )
for (i in 1:length(ARIs)){
  out_matrix[i,] = readRDS(paste0("data-raw/DLPFC3/", ARIs[i]))
}
rownames(out_matrix) = matrix(unlist(strsplit(ARIs, "[_.]")),ncol = 3, byrow = T)[,2]

out_matrix = cbind(out_matrix, Louvain = ARI_louvain)
out_matrix = out_matrix[,c(1:3, 10, 4:9)]
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

