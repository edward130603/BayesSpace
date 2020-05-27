##Note: this example runs on the old non-Rcpp code

.libPaths(.libPaths()[2])
library(scater)
library(scran)
library(mvnfast)
source("script.R")

slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')
print(slurm_arrayid)
sid = as.numeric(slurm_arrayid)

sce = readRDS("sce.RDS")
array = rbind(expand.grid(unique(sce$sample_name),"PCA", c(5,10,15)),
              expand.grid(unique(sce$sample_name),"TSNE", c(2,3)))#,
#expand.grid(unique(sce$sample_name),"UMAP", c(2,3)))
colnames(array) = c("sample_name", "input", "dim")
name = array$sample_name[sid]
input = array$input[sid]
n = array$dim[sid]

sce = sce[,sce$sample_name == name]

#Dimension reduction
set.seed(101)
dec <- modelGeneVarByPoisson(sce)
top <- getTopHVGs(dec, prop=0.1)

print("Running dimension reduction")
set.seed(102)
sce <- denoisePCA(sce, technical=dec, subset.row=top, min.rank = 15, max.rank = 50)
PCs = getDenoisedPCs(sce, technical=dec, subset.row=top, min.rank = 15, max.rank = 50)
sce <- runTSNE(sce, dimred="PCA", ncomponents = 2, name = "TSNE2")
sce <- runTSNE(sce, dimred="PCA", ncomponents = 3, name = "TSNE3")
#sce <- runUMAP(sce, dimred="PCA", ncomponents = 2, name = "UMAP2")
#sce <- runUMAP(sce, dimred="PCA", ncomponents = 3, name = "UMAP3")

#Setup clustering
df_PCA = data.frame(reducedDim(sce, "PCA"))[,1:15]
df_TSNE2 = data.frame(reducedDim(sce, "TSNE2"))
df_TSNE3 = data.frame(reducedDim(sce, "TSNE3"))
#df_UMAP2 = data.frame(reducedDim(sce, "UMAP2"))
#df_UMAP3 = data.frame(reducedDim(sce, "UMAP3"))
df = cbind(df_PCA, df_TSNE2, df_TSNE3)#, df_UMAP2, df_UMAP3)
colnames(df) = c(paste0("Y", 1:15, "_PCA"),
                 paste0("Y", 1:2 , "_TSNE2"),
                 paste0("Y", 1:3 , "_TSNE3"))#,
#paste0("Y", 1:2 , "_UMAP2"),
#paste0("Y", 1:3 , "_UMAP3"))
df$kmeans4 = kmeans(reducedDim(sce, "PCA"), centers = 4)$cluster
df$kmeans5 = kmeans(reducedDim(sce, "PCA"), centers = 5)$cluster
df$kmeans6 = kmeans(reducedDim(sce, "PCA"), centers = 6)$cluster
df$kmeans7 = kmeans(reducedDim(sce, "PCA"), centers = 7)$cluster
df$kmeans8 = kmeans(reducedDim(sce, "PCA"), centers = 8)$cluster
df$x = sce$imagecol
df$y = sce$imagerow
df$j = 1:nrow(df)
num_neighbors = sapply(sapply(1:nrow(df), function(x){df[(abs(df[,"x"] -df[x,"x"]) + abs(df[,"y"] - df[x,"y"])) <= 9 &
                                                           (abs(df[,"x"] -df[x,"x"]) + abs(df[,"y"] - df[x,"y"])) > 0,"j"]}), length)
df = df[num_neighbors>0,]
addCols = colnames(df)[(ncol(df)-7):ncol(df)]
rm(sce)
#Run clustering
print("Running clustering")
q = 4:8
for (i in 1:5){
  if (input == "PCA"){
    clust = run_mcmc_multi(df = df[,c(paste0("Y", 1:n, "_PCA"), addCols)],
                           gamma = 4,
                           q = q[i],
                           d = n,
                           nrep = 5000,
                           init = df[,paste0("kmeans", q[i])])
    saveRDS(clust, paste0("clust", name, "_PCA", n, "_q", q[i], ".RDS"))
  } else {
    clust = run_mcmc_multi(df = df[,c(paste0("Y", 1:n, "_", input, n), addCols)],
                           gamma = 4,
                           q = q[i],
                           d = n,
                           nrep = 5000,
                           init = df[,paste0("kmeans", q[i])])
    saveRDS(clust, paste0("clust", name, "_", input, n, "_q", q[i], ".RDS"))
  }
  print(paste0("Clustering complete for input= ", input, ", dim = ", n))
}
