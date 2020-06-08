.libPaths(.libPaths()[2])
library(scater)
library(scran)
library(mclust)
source("script.R")

slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')
print(slurm_arrayid)
sid = as.numeric(slurm_arrayid)

sce = readRDS("sce.RDS")
array = unique(sce$sample_name)
name = array[sid]

sce = sce[,sce$sample_name == name]
q = length(unique(na.omit(sce$layer_guess)))

#Dimension reduction
set.seed(101)
dec <- modelGeneVarByPoisson(sce)
top <- getTopHVGs(dec, prop=0.1)

message("Running dimension reduction")
set.seed(102)
sce <- denoisePCA(sce, technical=dec, subset.row=top, min.rank = 15, max.rank = 50)
PCs = getDenoisedPCs(sce, technical=dec, subset.row=top, min.rank = 15, max.rank = 50)

#Setup clustering
Y = PCs$components[,1:15]
set.seed(103)
km = kmeans(Y, centers = q)$cluster
set.seed(104)
mclust_eee = Mclust(Y, G = q, modelNames = "EEE")$classification
set.seed(105)
mclust_vvv = Mclust(Y, G = q, modelNames = "VVV")$classification


positions = cbind(sce$imagecol, sce$imagerow) #positions
colnames(positions) = c("x", "y") 
xdist = coef(lm(sce$imagecol~sce$col))[2] #x distance between neighbors
ydist = coef(lm(sce$imagerow~sce$row))[2] #y distance between neighbors
dist = xdist + ydist + 0.2
truth = sce$layer_guess
snn = sce$HVG_PCA_spatial
markers = sce$markers_PCA_spatial
rm(sce)
#Normal model
clust1 = cluster(Y= Y, positions = positions, 
                 model = "normal", precision = "equal",
                 q = q, init = mclust_eee, nrep = 50000, gamma = 3, dist = dist)
clust1col = clust1$labels
rm(clust1)
clust2 = cluster(Y= Y, positions = positions, 
                 model = "t", precision = "equal",
                 q = q, init = mclust_eee, nrep = 50000, gamma = 3, dist = dist)
clust2col = clust2$labels
rm(clust2)
clust3 = cluster(Y= Y, positions = positions, 
                 model = "normal", precision = "variable",
                 q = q, init = mclust_eee, nrep = 50000, gamma = 3, dist = dist)
clust3col = clust3$labels
rm(clust3)
clust4 = cluster(Y= Y, positions = positions, 
                 model = "t", precision = "variable",
                 q = q, init = mclust_eee, nrep = 50000, gamma = 3, dist = dist)
clust4col = clust4$labels
rm(clust4)
out = c(markers = adjustedRandIndex(truth, markers),
        snn = adjustedRandIndex(truth, snn),
        km = adjustedRandIndex(truth, km),
        mclust_eee = adjustedRandIndex(truth, mclust_eee),
        mclust_vvv = adjustedRandIndex(truth, mclust_vvv),
        spatial_ne = adjustedRandIndex(truth, clust1col),
        spatial_te = adjustedRandIndex(truth, clust2col),
        spatial_nv = adjustedRandIndex(truth, clust3col),
        spatial_tv = adjustedRandIndex(truth, clust4col))
filename = paste0("ari_", name, ".RDS")
saveRDS(out, filename)