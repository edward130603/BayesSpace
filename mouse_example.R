library(DropletUtils)
library(scater)
library(annotables)
library(ggplot2)
library(ggbeeswarm)
library(tibble)
library(AnnotationHub)

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
qc.lib = isOutlier(df_qc$sum, log=TRUE, type="lower")
attr(qc.lib, "thresholds")
qc.nexprs = isOutlier(df_qc$detected, log=TRUE, type="lower")
attr(qc.nexprs, "thresholds")
qc.mito = isOutlier(df_qc$subsets_Mito_percent, type="higher")
attr(qc.mito, "thresholds")
discard = qc.lib | qc.nexprs | qc.mito
