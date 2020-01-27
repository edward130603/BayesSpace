library(DropletUtils)
library(scater)
library(ggplot2)
library(ggbeeswarm)
library(tibble)
library(AnnotationHub)
library(scran)
library(SingleR)
library(viridis)
library(patchwork)


#Make sce
sce_A2 = read10xCounts("data-raw/A2")
sce_A2$Sample = "A2"
rownames(sce_A2) = uniquifyFeatureNames(rowData(sce_A2)$ID, rowData(sce_A2)$Symbol)
colnames(sce_A2) = paste0(sce_A2$Sample, '.', sce_A2$Barcode)

#Add position
pos = read.csv("data-raw/A2/tissue_positions_list.txt", header=FALSE)
colnames(pos) = c("Barcode", "tissue", "Y1", "X1", "Y2", "X2")
colData(sce_A2) = merge(colData(sce_A2), pos, by = "Barcode")

#Preprocess
ens.mm.v97 = AnnotationHub()[["AH73905"]]
rowData(sce_A2)$SEQNAME = mapIds(ens.mm.v97, keys=rowData(sce_A2)$ID,
                                   keytype="GENEID", column="SEQNAME")

is.mito = which(rowData(sce_A2)$SEQNAME == "MT")
is.RBC = grep("^Hb", rownames(sce_A2))
df_qc = perCellQCMetrics(sce_A2, subsets = list(Mito = is.mito,
                                          RBC = is.RBC))
qc = quickPerCellQC(df_qc, percent_subsets = "subsets_Mito_percent")
sce_A2$discard = qc$discard
ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = discard)) + 
  geom_point(size = 4)

#Normalization
set.seed(100)
clusters = quickCluster(sce_A2)
sce_A2 = computeSumFactors(sce_A2, clusters = clusters)
sce_A2 = logNormCounts(sce_A2)
plot(librarySizeFactors(sce_A2), sizeFactors(sce_A2), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", log="xy", 
     col = ifelse(sce_A2$discard, "red", "black"))

set.seed(101)
dec <- modelGeneVarByPoisson(sce_A2)
top <- getTopHVGs(dec, prop=0.1)
plot(dec$mean, dec$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)

#Dim reduction and clustering
set.seed(102)
sce_A2 <- denoisePCA(sce_A2, technical=dec, subset.row=top)
sce_A2 <- runTSNE(sce_A2, dimred="PCA")

snn.gr <- buildSNNGraph(sce_A2, use.dimred="PCA", k=100)
sce_A2$cluster <- factor(igraph::cluster_walktrap(snn.gr)$membership)
table(sce_A2$cluster)

tsne1 = ggplot(as.data.frame(reducedDim(sce_A2, "TSNE")), aes(x = V1, y = V2, col = sce_A2$cluster))+
  geom_point() +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()
  
spatial1 = ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = cluster)) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "Cluster") +
  theme_classic()

tsne1 + spatial1 +
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "A")
  
ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = logcounts(sce_A2)["Hpca",])) + 
  geom_point(size = 4) + scale_color_viridis(option = "A")
pc1 = ggplot(as.data.frame(reducedDim(sce_A2, "PCA")), aes(x = PC1))+geom_histogram(bins = 20)+
  facet_wrap(~sce_A2$cluster, scales = "free") +
  labs(y = NULL) +
  theme_classic()

#SingleR annotation
mouse_ref = MouseRNAseqData()
pred_main = SingleR(test = sce_A2, ref = mouse_ref, labels = mouse_ref$label.main)
pred_main$pruned.labels[is.na(pred_main$pruned.labels)] = "Unknown"
spatial2 = ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = pred_main$pruned.labels)) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "SingleR main label")+
  scale_color_manual(values = c("Astrocytes" = "#de2d26", "Epithelial cells" = "#756bb1",
                                "Fibroblasts" = "#e6550d", "Neurons" = "#3182bd",
                                "Oligodendrocytes" = "#31a354", "Unknown" = "#000000")) +
  theme_classic()

tsne2 = ggplot(as.data.frame(reducedDim(sce_A2, "TSNE")), aes(x = V1, y = V2, col = pred_main$pruned.labels))+
  geom_point() +
  labs(x = NULL, y = NULL, color = "SingleR main label") +
  scale_color_manual(values = c("Astrocytes" = "#de2d26", "Epithelial cells" = "#756bb1",
                                "Fibroblasts" = "#e6550d", "Neurons" = "#3182bd",
                                "Oligodendrocytes" = "#31a354", "Unknown" = "#000000")) +
  theme_classic()


pred_fine = SingleR(test = sce_A2, ref = mouse_ref, labels = mouse_ref$label.fine)
pred_fine$pruned.labels[is.na(pred_fine$pruned.labels)] = "Unknown"
spatial3 = ggplot(as.data.frame(colData(sce_A2)), aes(x = X1, y = Y1, col = pred_fine$pruned.labels)) + 
  geom_point() +
  labs(x = NULL, y = NULL, color = "SingleR fine label")+
  scale_color_manual(values = c("aNSCs" = "#6baed6", "Astrocytes" = "#de2d26", "Ependymal" = "#756bb1",
                                "Neurons" = "#3182bd", "NPCs" = "#08519c",
                                "Oligodendrocytes" = "#31a354", "qNSCs" = "#bdd7e7","Unknown" = "#000000")) +
  theme_classic()

tsne3 = ggplot(as.data.frame(reducedDim(sce_A2, "TSNE")), aes(x = V1, y = V2, col = pred_fine$pruned.labels))+
  geom_point() +
  labs(x = NULL, y = NULL, color = "SingleR fine label") +
  scale_color_manual(values = c("aNSCs" = "#6baed6", "Astrocytes" = "#de2d26", "Ependymal" = "#756bb1",
                                "Neurons" = "#3182bd", "NPCs" = "#08519c",
                                "Oligodendrocytes" = "#31a354", "qNSCs" = "#bdd7e7","Unknown" = "#000000")) +
  theme_classic()

(tsne2 + spatial2 + plot_layout(guides = 'collect')) /
  (tsne3 + spatial3 + plot_layout(guides = 'collect'))+
  plot_annotation(tag_levels = "A")


#Next steps
#annotation with allen mouse brain atlas reference
#https://github.com/AllenInstitute/scrattch.io
#https://portal.brain-map.org/atlases-and-data/rnaseq

#Principal components
pc2 = ggplot(as.data.frame(reducedDim(sce_A2, "PCA")), aes(x = PC2))+geom_histogram(bins = 20)+
  facet_wrap(~sce_A2$cluster, scales = "free") +
  labs(y = NULL) +
  theme_classic()
pc1 / pc2 + plot_annotation(tag_levels = "A")
