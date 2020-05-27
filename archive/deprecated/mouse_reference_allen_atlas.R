library(scrattch.io)
options(stringsAsFactors = FALSE)
setwd("/fh/fast/gottardo_r/ezhao_working/transcriptome/")
tome        = "transcrip.tome"
exons       = read_tome_dgCMatrix(tome,"data/t_exon")    # (or data/exon)
#introns     = read_tome_dgCMatrix(tome,"data/t_intron")  # (or data/intron)
sample_name = read_tome_sample_names(tome)  
gene_name   = read_tome_gene_names(tome)

library(SingleCellExperiment)
library(tibble)
library(scater)
library(scran)
library(SingleR)
sce = SingleCellExperiment(assays = list(counts = exons))
rownames(sce) = gene_name
colnames(sce) = sample_name

atlas = read.csv("sample_annotations.csv")
gene_info = as_tibble(read.csv("mouse_VISp_2018-06-14_genes-rows.csv"))

rowData(sce) = gene_info
rowData(sce)
atlas2 = atlas[match(sample_name, atlas$sample_name),] 
colData(sce) = DataFrame(atlas2)

gene_list = readRDS("sce_genes.RDS")
sce2 = sce[rowData(sce)$gene_symbol %in% gene_list$Symbol,sce$class_label != "Exclude"]
df_qc = perCellQCMetrics(sce2)
qc = quickPerCellQC(df_qc)
sce2 = sce2[,!qc$discard]
sce2 = logNormCounts(sce2)

dec = modelGeneVarByPoisson(sce2)
sce2 = sce2[dec$bio > 0.2, ]

colData(sce2) %>%
  as_tibble() %>%
  rowid_to_column("index") %>%
  group_by(region_label) %>%
  sample_n(min(n(),200)) ->
  sample_region
sce_region = sce2[,sample_region$index]

colData(sce2) %>%
  as_tibble() %>%
  rowid_to_column("index") %>%
  group_by(subclass_label) %>%
  sample_n(min(n(),150)) ->
  sample_subclass
sce_subclass = sce2[,sample_subclass$index]

sce_A2 = readRDS("sce_A2.RDS")
#cluster, region
pred = SingleR(test = sce_A2, ref = sce_region, labels = sce_region$region_label, method = "cluster",
               clusters = sce_A2$cluster)
pred$pruned.labels[is.na(pred$pruned.labels)] = "Unknown"
merge_data = data.frame(cluster = factor(1:7), label = pred$pruned.labels)
merge_data_long = left_join(data.frame(cluster = sce_A2$cluster), merge_data)
sce_A2$singleR_region_cluster = merge_data_long$label
spatial1 = ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = singleR_region_cluster)) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "Region")+
  theme_classic()

tsne1 = ggplot(as.data.frame(reducedDim(sce_A2, "TSNE")), aes(x = V1, y = V2, col = sce_A2$singleR_region_cluster))+
  geom_point() +
  labs(x = NULL, y = NULL, color = "Region") +
  theme_classic()

tsne1 + spatial1 +
  plot_layout(guides = 'collect')
#cluster, subclass
pred2 = SingleR(test = sce_A2, ref = sce_subclass, labels = sce_subclass$subclass_label, method = "cluster",
                clusters = sce_A2$cluster)
pred2$pruned.labels[is.na(pred2$pruned.labels)] = "Unknown"
merge_data = data.frame(cluster = factor(1:7), label = pred2$pruned.labels)
merge_data_long = left_join(data.frame(cluster = sce_A2$cluster), merge_data)
sce_A2$singleR_subclass_cluster = merge_data_long$label
spatial1 = ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = singleR_subclass_cluster)) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "Subclass")+
  theme_classic()

tsne1 = ggplot(as.data.frame(reducedDim(sce_A2, "TSNE")), aes(x = V1, y = V2, col = sce_A2$singleR_subclass_cluster))+
  geom_point() +
  labs(x = NULL, y = NULL, color = "Subclass") +
  theme_classic()

tsne1 + spatial1 +
  plot_layout(guides = 'collect')
pred_single = SingleR(test = sce_A2, ref = sce_region, labels = sce_region$region_label)
pred_single$pruned.labels[is.na(pred_single$pruned.labels)] = "Unknown"

spatial1 = ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = pred_single$pruned.labels)) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "Region")+
  theme_classic()

tsne1 = ggplot(as.data.frame(reducedDim(sce_A2, "TSNE")), aes(x = V1, y = V2, col = pred_single$pruned.labels))+
  geom_point() +
  labs(x = NULL, y = NULL, color = "Region") +
  theme_classic()

tsne1 + spatial1 +
  plot_layout(guides = 'collect')

pred_single2 = SingleR(test = sce_A2, ref = sce_subclass, labels = sce_subclass$subclass_label)
pred_single2$pruned.labels[is.na(pred_single2$pruned.labels)] = "Unknown"

spatial1 = ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = pred_single2$pruned.labels)) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "Region")+
  theme_classic()

tsne1 = ggplot(as.data.frame(reducedDim(sce_A2, "TSNE")), aes(x = V1, y = V2, col = pred_single2$pruned.labels))+
  geom_point() +
  labs(x = NULL, y = NULL, color = "Region") +
  theme_classic()

tsne1 + spatial1 +
  plot_layout(guides = 'collect')
