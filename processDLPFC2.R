library(dplyr)
library(ggplot2)
library(reshape)
ARIs = list.files("data-raw/DLPFC3/PC10")
out_matrix = matrix(nrow = length(ARIs), ncol = 9)
colnames(out_matrix) = c("SNN (Markers)", "SNN (HVG)", "k-means", "mclust (EEE)", "mclust (VVV)", 
                         "SC (normal/equal)",  "SC (t/equal)",  "SC (normal/variable)",  "SC (t/variable)" )
for (i in 1:length(ARIs)){
  out_matrix[i,] = readRDS(paste0("data-raw/DLPFC3/PC10/", ARIs[i]))
}
rownames(out_matrix) = matrix(unlist(strsplit(ARIs, "[_.]")),ncol = 3, byrow = T)[,2]
data_out = melt(out_matrix)
colnames(data_out) = c("sample", "method", "ARI")
data_out$method = factor(data_out$method, 
                         levels = colnames(out_matrix))
ggplot(data_out, aes(x = method, y = ARI))+
  geom_boxplot() +
  geom_point(aes(color = factor(sample))) +
  labs(x = "Method", color = "Sample")+
  theme_light()+
  coord_cartesian(ylim = c(0,0.6))

ARIs = list.files("data-raw/DLPFC3/")[1:12]
out_matrix = matrix(nrow = length(ARIs), ncol = 9)
colnames(out_matrix) = c("SNN (Markers)", "SNN (HVG)", "k-means", "mclust (EEE)", "mclust (VVV)", 
                         "SC (normal/equal)",  "SC (t/equal)",  "SC (normal/variable)",  "SC (t/variable)" )
for (i in 1:length(ARIs)){
  out_matrix[i,] = readRDS(paste0("data-raw/DLPFC3/", ARIs[i]))
}
rownames(out_matrix) = matrix(unlist(strsplit(ARIs, "[_.]")),ncol = 3, byrow = T)[,2]
data_out = melt(out_matrix)
colnames(data_out) = c("sample", "method", "ARI")
data_out$method = factor(data_out$method, 
                         levels = colnames(out_matrix))
ggplot(data_out, aes(x = method, y = ARI))+
  geom_boxplot() +
  geom_point(aes(color = factor(sample))) +
  labs(x = "Method", color = "Sample")+
  theme_light()+
  coord_cartesian(ylim = c(0,0.6)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
