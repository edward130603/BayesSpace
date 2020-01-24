library(DropletUtils)
library(scater)
library(ggplot2)
library(ggbeeswarm)
library(tibble)
library(AnnotationHub)
library(scran)

#Make sce
sce_A2 = read10xCounts("data-raw/A2")
sce_A2$Sample = "A2"
rownames(sce_A2) = uniquifyFeatureNames(rowData(sce_A2)$ID, rowData(sce_A2)$Symbol)
colnames(sce_A2) = paste0(sce_A2$Sample, '.', sce_A2$Barcode)

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

#Normalization
clusters = quickCluster(sce_A2)
sce_A2 = computeSumFactors(sce_A2, clusters = clusters)
sce_A2 = logNormCounts(sce_A2)
plot(librarySizeFactors(sce_A2), sizeFactors(sce_A2), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", log="xy", 
     col = ifelse(sce_A2$discard, "red", "black"))

dec <- modelGeneVarByPoisson(sce_A2)
top <- getTopHVGs(dec, prop=0.1)
plot(dec$mean, dec$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)

#Dim reduction and clustering
sce_A2 <- denoisePCA(sce_A2, technical=dec, subset.row=top)
sce_A2 <- runTSNE(sce_A2, dimred="PCA")

snn.gr <- buildSNNGraph(sce_A2, use.dimred="PCA", k=25)
sce_A2$cluster <- factor(igraph::cluster_walktrap(snn.gr)$membership)
table(sce_A2$cluster)
plotTSNE(sce_A2, colour_by="cluster")

