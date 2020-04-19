library(spatialLIBD)
library(scran)
library(ggbeeswarm)
library(tidyverse)

#Load sce
sce = readRDS("data/sce.RDS")
samples = unique(sce$sample_name)

#Preprocess sce
sce_list = vector(mode = "list", length = length(samples))
for (i in 1:length(samples)){
  name = samples[i]
  df = data.frame(reducedDim(sce, "PCA"))[sce$sample_name == name, 1:5]
  df$x = sce$imagecol[sce$sample_name == name]
  df$y = sce$imagerow[sce$sample_name == name]
  df$j = 1:nrow(df)
  num_neighbors = sapply(sapply(1:nrow(df), function(x){df[(abs(df[,"x"] -df[x,"x"]) + abs(df[,"y"] - df[x,"y"])) <= 9 &
                                                                           (abs(df[,"x"] -df[x,"x"]) + abs(df[,"y"] - df[x,"y"])) > 0,"j"]}), length)
  print(sum(num_neighbors==0))
  sce_list[i] = sce[,sce$sample_name == name][,num_neighbors>0]
}
rm(sce)
sce = do.call("cbind", sce_list)

data1 = data.frame(name = sce$sample_name,
                   x = sce$imagecol,
                   y = sce$imagerow,
                   truth = sce$layer_guess_reordered,
                   unsupervised = sce$HVG_PCA_spatial,
                   markers = sce$markers_PCA_spatial,
                   semisupervised = sce$pseudobulk_PCA_spatial)
rm(sce)

#Load and process mcmc
file_list = list.files(path = "data-raw/DLPFC2/")
z = vector(mode = "list", length = length(file_list))
avgCor = rep(NA, length = length(file_list))
avgCor3 = rep(NA, length = length(file_list))
logLik = matrix(nrow = length(file_list), ncol = 500)
out = data.frame(matrix(nrow = 300, ncol = 6))
colnames(out) = c("i", "sample", "input", "dim", "clusters", "ARI")
for (i in 1:length(file_list)){
  temp_data = readRDS(paste0("data-raw/DLPFC2/",file_list[i]))
  z[[i]] = apply(temp_data$z[4000:5000,], 2, Mode)
  temp_cor = cov2cor(Reduce(`+`, temp_data$lambda[4000:5000]))
  avgCor[i] = mean(abs(temp_cor[upper.tri(temp_cor)]))
  #avgCor3[i] = mean(abs(temp_cor[upper.tri(temp_cor)][1:min(3, length(temp_cor[upper.tri(temp_cor)]))]))
  logLik[i,] = temp_data$logLik
  
  temp_row = unlist(strsplit(file_list[i], c("_|clust|q|\\.|RDS")))
  out[i,] = c(i, temp_row[2], sub("[^[:alpha:]]+", "", temp_row[3]),
              sub("[^[:digit:]]+", "", temp_row[3]), temp_row[5],
              mclust::adjustedRandIndex(data1$truth[data1$name == temp_row[2]], z[[i]]))
  print(i)
}
out$dim = paste0("Dimensions: ", out$dim)

#Process reference clustering
out_append = as_tibble(matrix(nrow = length(samples), ncol = 7))
colnames(out_append) = c("i", "sample", "input", "dim", "Unsupervised", "Markers", "Semi-supervised")
for (i in 1:length(samples)){
  out_append[i,] = c(NA, sample = as.character(samples[i]),
                     input = "Maynard et al., 2020", dim = "",
                     Unsupervised = mclust::adjustedRandIndex(data1$truth[data1$name == samples[i]], 
                                                              data1$unsupervised[data1$name == samples[i]]),
                     Markers = mclust::adjustedRandIndex(data1$truth[data1$name == samples[i]], 
                                                              data1$markers[data1$name == samples[i]]),
                     `Semi-supervised` = mclust::adjustedRandIndex(data1$truth[data1$name == samples[i]], 
                                                              data1$semisupervised[data1$name == samples[i]]))
}
out_append %>%
  gather("clusters", "ARI", Unsupervised:`Semi-supervised`) %>%
  rbind(out) ->
  out2

out_snn_k10 = data.frame(i = NA, sample = samples, input = "SNN", dim = "", clusters = NA, ARI = NA)
for (i in 1:nrow(out_snn_k10)){
  snn.gr = buildSNNGraph(sce_list[[i]], use.dimred="PCA", k=10)
  snn = factor(igraph::cluster_walktrap(snn.gr)$membership)
  out_snn_k10[i, "clusters"] = length(unique(snn))
  out_snn_k10[i, "ARI"] = mclust::adjustedRandIndex(snn, sce_list[[i]]$layer_guess)
}

out_snn_k30 = data.frame(i = NA, sample = samples, input = "SNN", dim = "", clusters = NA, ARI = NA)
for (i in 1:nrow(out_snn_k30)){
  snn.gr = buildSNNGraph(sce_list[[i]], use.dimred="PCA", k=30)
  snn = factor(igraph::cluster_walktrap(snn.gr)$membership)
  out_snn_k30[i, "clusters"] = length(unique(snn))
  out_snn_k30[i, "ARI"] = mclust::adjustedRandIndex(snn, sce_list[[i]]$layer_guess)
}

out_snn_k50 = data.frame(i = NA, sample = samples, input = "SNN", dim = "", clusters = NA, ARI = NA)
for (i in 1:nrow(out_snn_k50)){
  snn.gr = buildSNNGraph(sce_list[[i]], use.dimred="PCA", k=50)
  snn = factor(igraph::cluster_walktrap(snn.gr)$membership)
  out_snn_k50[i, "clusters"] = length(unique(snn))
  out_snn_k50[i, "ARI"] = mclust::adjustedRandIndex(snn, sce_list[[i]]$layer_guess)
}

out_snn_k10$clusters = "k = 10"
out_snn_k30$clusters = "k = 30"
out_snn_k50$clusters = "k = 50"

out_kmeans = data.frame(i = NA, sample = rep(samples,each = 5), input = "K-means", dim = "", clusters = 4:8, ARI = NA)
for (i in 1:nrow(out_kmeans)){
  kmeans = kmeans(reducedDim(sce_list[[floor((i-1)/5 + 1)]]), centers = out_kmeans[i, "clusters"])$cluster
  out_kmeans$ARI[i] = mclust::adjustedRandIndex(kmeans, sce_list[[floor((i-1)/5 + 1)]]$layer_guess)
}

out3 = rbind(out2, out_kmeans, out_snn_k10, out_snn_k30, out_snn_k50)

#Plot
ggplot(out2, aes(x = factor(clusters, levels = c(4:8, "Markers", "Unsupervised", "Semi-supervised")), y = as.numeric(ARI), color = sample)) +
  geom_boxplot(aes(group = clusters)) +
  geom_beeswarm(priority = "descending") +
  facet_wrap(nrow = 1, ~factor(input, levels = c("Maynard et al., 2020", "TSNE", "PCA"))+
               factor(dim, levels = c("", "Dimensions: 2", "Dimensions: 3", "Dimensions: 5", "Dimensions: 10", "Dimensions: 15")), scales = "free_x") +
  labs(x = "Clusters", y = "ARI", color = "Sample")

ggplot(out2[out2$input != "TSNE",], aes(x = factor(clusters, levels = c(4:8, "Markers", "Unsupervised", "Semi-supervised")), y = as.numeric(ARI), color = (sample %in% 151669:151672))) +
  geom_boxplot(aes(group = clusters)) +
  geom_beeswarm(priority = "descending") +
  facet_wrap(nrow = 1, ~factor(input, levels = c("Maynard et al., 2020", "TSNE", "PCA"))+
               factor(dim, levels = c("", "Dimensions: 2", "Dimensions: 3", "Dimensions: 5", "Dimensions: 10", "Dimensions: 15")), scales = "free_x") +
  labs(x = "Clusters", y = "ARI", color = "Sample")+
  scale_y_continuous(breaks = seq(0,0.7,0.1))

ggplot(out2[out2$input != "PCA",], aes(x = factor(clusters, levels = c(4:8, "Markers", "Unsupervised", "Semi-supervised")), y = as.numeric(ARI), color = sample)) +
  geom_boxplot(aes(group = clusters)) +
  geom_beeswarm(priority = "descending") +
  facet_wrap(nrow = 1, ~factor(input, levels = c("Maynard et al., 2020", "TSNE", "PCA"))+
               factor(dim, levels = c("", "Dimensions: 2", "Dimensions: 3", "Dimensions: 5", "Dimensions: 10", "Dimensions: 15")), scales = "free_x") +
  labs(x = "Clusters", y = "ARI", color = "Sample")+
  scale_y_continuous(breaks = seq(0,0.7,0.1))

ggplot(out2, aes(x = factor(clusters, levels = c(4:8, "Markers", "Unsupervised", "Semi-supervised")), y = as.numeric(ARI), color = sample)) +
  geom_boxplot(aes(group = clusters)) +
  geom_beeswarm(priority = "descending") +
  facet_wrap(nrow = 1, ~factor(input, levels = c("Maynard et al., 2020", "TSNE", "PCA"))+
               factor(dim, levels = c("", "Dimensions: 2", "Dimensions: 3", "Dimensions: 5", "Dimensions: 10", "Dimensions: 15")), scales = "free_x") +
  labs(x = "Clusters", y = "ARI", color = "Sample")

out3$input = factor(out3$input, levels = c("Maynard et al., 2020", "SNN", "K-means","PCA", "TSNE"))
ggplot(out3[!out3$input %in% c("PCA", "TSNE"),], aes(x = factor(clusters, levels = c(4:8, "Markers", "Unsupervised", "Semi-supervised",
                                                 "k = 10", "k = 30", "k = 50")), y = as.numeric(ARI), color = sample)) +
  geom_boxplot(aes(group = clusters)) +
  geom_beeswarm(priority = "descending") +
  facet_wrap(~factor(input)+
               factor(dim, levels = c("", "Dimensions: 2", "Dimensions: 3", "Dimensions: 5", "Dimensions: 10", "Dimensions: 15")), scales = "free_x") +
  labs(x = "Clusters", y = "ARI", color = "Sample")


out$logLik = rowMeans(logLik[,400:500])
ggplot(out, aes(x = factor(clusters), y = logLik, color = sample)) +
  geom_line(aes(group = sample)) +
  geom_point(size = 2) +
  facet_wrap(nrow = 1, ~factor(input, levels = c("TSNE", "PCA"))+
               factor(dim, levels = c("Dimensions: 2", "Dimensions: 3", "Dimensions: 5", "Dimensions: 10", "Dimensions: 15")), scales = "free") +
  labs(x = "Clusters", y = "Pseudo-log-likelihood", color = "Sample")

ggplot(out, aes(x = factor(clusters), y = as.numeric(ARI), color = sample)) +
  geom_line(aes(group = sample)) +
  geom_point(size = 2) +
  facet_wrap(nrow = 1, ~factor(input, levels = c("TSNE", "PCA"))+
               factor(dim, levels = c("Dimensions: 2", "Dimensions: 3", "Dimensions: 5", "Dimensions: 10", "Dimensions: 15")), scales = "free") +
  labs(x = "Clusters", y = "ARI", color = "Sample")

ggplot(out, aes(x = factor(clusters), y = avgCor, color = sample)) +
  geom_line(aes(group = sample)) +
  geom_point(size = 2) +
  facet_wrap(nrow = 1, ~factor(input, levels = c("TSNE", "PCA"))+
               factor(dim, levels = c("Dimensions: 2", "Dimensions: 3", "Dimensions: 5", "Dimensions: 10", "Dimensions: 15")), scales = "fixed") +
  labs(x = "Clusters", y = "Mean magnitude of correlation coefficient", color = "Sample")

ggplot(out, aes(x = factor(clusters), y = avgCor3, color = sample)) +
  geom_line(aes(group = sample)) +
  geom_point(size = 2) +
  facet_wrap(nrow = 1, ~factor(input, levels = c("TSNE", "PCA"))+
               factor(dim, levels = c("Dimensions: 2", "Dimensions: 3", "Dimensions: 5", "Dimensions: 10", "Dimensions: 15")), scales = "fixed") +
  labs(x = "Clusters", y = "Mean magnitude of correlation coefficient", color = "Sample")

ggplot(out, aes(y = as.numeric(ARI), x = logLik, color = sample)) +
  geom_smooth(aes(group = sample), method = "lm", se = F) +
  geom_point(size = 2) +
  facet_wrap(nrow = 1, ~factor(input, levels = c("TSNE", "PCA"))+
               factor(dim, levels = c("Dimensions: 2", "Dimensions: 3", "Dimensions: 5", "Dimensions: 10", "Dimensions: 15")), scales = "free") +
  labs(y = "ARI", x = "Pseudo-log-likelihood", color = "Sample")
