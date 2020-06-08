#library(devtools)
#install_github("LieberInstitute/spatialLIBD")
library(tidyverse)
library(spatialLIBD)
library(scater)
library(scran)
library(extrafont)
library(patchwork)
Rcpp::sourceCpp('script.cpp')

sce <- fetch_data(type = 'sce')
sce_image_clus(
  sce = sce,
  clustervar = 'layer_guess_reordered',
  sampleid = '151673',
  colors = libd_layer_colors,
  ... = ' DLPFC Human Brain Layers\nMade with github.com/LieberInstitute/spatialLIBD'
)

sce = sce[,sce$sample_name == "151673"]

#Dim reduction and clustering
set.seed(101)
dec <- modelGeneVarByPoisson(sce)
top <- getTopHVGs(dec, prop=0.1)
plot(dec$mean, dec$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)

set.seed(102)
sce <- denoisePCA(sce, technical=dec, subset.row=top)
PCs = getDenoisedPCs(sce, technical=dec, subset.row=top, min.rank = 5, max.rank = 50)
sce <- runTSNE(sce, dimred="PCA")
#sce2 = sce[,sce$col<25 & sce$row < 20]


#Spatial clustering

#choose input
km = sapply(1:10,function(k){kmeans(PCs$components, centers = k)$cluster})
positions = cbind(sce$imagecol, sce$imagerow) #save positions as df
colnames(positions) = c("x", "y") 
xdist = coef(lm(sce$imagecol~sce$col))[2] #x distance between neighbors
ydist = coef(lm(sce$imagerow~sce$row))[2] #y distance between neighbors
dist = xdist + ydist + 0.2
#Normal model
clust1 = cluster(Y= PCs$components, positions = as.matrix(positions), q = 7, init = km[,7], nrep = 10000, gamma = 1.5, dist = dist)
clust1g3 = cluster(Y= PCs$components, positions = as.matrix(positions), q = 7, init = km[,7], nrep = 5000, gamma = 3, dist = dist)
clust1col = apply(clust1$z[9000:10000,], 2, Mode)
clust1alpha = pmax(colMeans(clust1$z[9000:10000,]==1),
                   colMeans(clust1$z[9000:10000,]==2),
                   colMeans(clust1$z[9000:10000,]==3),
                   colMeans(clust1$z[9000:10000,]==4),
                   colMeans(clust1$z[9000:10000,]==5),
                   colMeans(clust1$z[9000:10000,]==6),
                   colMeans(clust1$z[9000:10000,]==7))
clust1_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(clust1col, levels = c(7,1,4,3,6,2,5))), alpha = clust1alpha,
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F, col = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()

truth_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(sce$layer_guess_reordered)),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()

km_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(km[,7], c(7,1,4,3,6,2,5))),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F, col = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()

maynardclust_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(sce$HVG_PCA_spatial, levels = c(5,1,3,2,7,8,6,4))),
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F, col = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()

truth_out + clust_out #Truth vs spatial clustering
truth_out + maynardclust_out #Truth vs best clustering implementation from Maynard

#t model
clust2 = cluster(Y= PCs$components, positions = as.matrix(positions), model = "t", q = 7, init = km[,7], nrep = 10000, gamma = 1.5, dist = dist)
clust2col = apply(clust2$z[9000:10000,], 2, Mode)
clust2alpha = pmax(colMeans(clust2$z[9000:10000,]==1),
                   colMeans(clust2$z[9000:10000,]==2),
                   colMeans(clust2$z[9000:10000,]==3),
                   colMeans(clust2$z[9000:10000,]==4),
                   colMeans(clust2$z[9000:10000,]==5),
                   colMeans(clust2$z[9000:10000,]==6),
                   colMeans(clust2$z[9000:10000,]==7))
clust2_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(clust2col, levels = c(7,1,4,3,6,2,5))), alpha = clust2alpha,
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F, col = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()

###variable normal model
clust3 = cluster(Y= PCs$components, positions = as.matrix(positions), precision = "variable", q = 7, init = km[,7], nrep = 10000, gamma = 1.5, dist = dist)
clust3col = apply(clust3$z[9000:10000,], 2, Mode)
clust3alpha = pmax(colMeans(clust3$z[9000:10000,]==1),
                   colMeans(clust3$z[9000:10000,]==2),
                   colMeans(clust3$z[9000:10000,]==3),
                   colMeans(clust3$z[9000:10000,]==4),
                   colMeans(clust3$z[9000:10000,]==5),
                   colMeans(clust3$z[9000:10000,]==6),
                   colMeans(clust3$z[9000:10000,]==7))
clust3_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(clust3col, levels = c(7,1,4,3,6,2,5))), alpha = clust3alpha,
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F, col = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()

###variable t model
clust4 = cluster(Y= PCs$components, positions = as.matrix(positions), model = "t", precision = "variable", q = 7, init = km[,7], nrep = 10000, gamma = 1.5, dist = dist)
clust4col = apply(clust4$z[9000:10000,], 2, Mode)
clust4alpha = pmax(colMeans(clust4$z[9000:10000,]==1),
                   colMeans(clust4$z[9000:10000,]==2),
                   colMeans(clust4$z[9000:10000,]==3),
                   colMeans(clust4$z[9000:10000,]==4),
                   colMeans(clust4$z[9000:10000,]==5),
                   colMeans(clust4$z[9000:10000,]==6),
                   colMeans(clust4$z[9000:10000,]==7))
clust4_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(clust4col, levels = c(7,1,4,3,6,2,5))), alpha = clust4alpha,
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F, col = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()

((truth_out | maynardclust_out | km_out) + plot_layout(guides = "collect"))/(clust1_out | clust2_out|clust3_out|clust4_out) 

#calculate ARI
mclust::adjustedRandIndex(sce$layer_guess, clust1col)
mclust::adjustedRandIndex(sce$layer_guess, clust2col)
mclust::adjustedRandIndex(sce$layer_guess, clust3col)
mclust::adjustedRandIndex(sce$layer_guess, clust4col)
mclust::adjustedRandIndex(clust3col, clust4col)
mclust::adjustedRandIndex(sce$layer_guess, sce$HVG_PCA_spatial)
mclust::adjustedRandIndex(sce$layer_guess, km[,7])
mclust::adjustedRandIndex(clust2col, clust1col)

#Spatial deconvolution
ptm = proc.time()
deconv1 = deconvolve(Y = PCs$components, positions = positions, nrep = 10000, gamma = 2, xdist = xdist, ydist = ydist, init = clust1col, q = 7)
proc.time()-ptm #10000 reps: 5082 seconds runtime
saveRDS(deconv1, "data-raw/deconv_sim1.RDS")
deconv1col = apply(deconv1$z[seq(1000,10000,10),], 2, Mode)
deconv1mu = colMeans(deconv1$mu[1000:10000,])
deconv1mu_matrix = matrix(deconv1mu, byrow = T, ncol = ncol(deconv1$lambda[[1]]))
deconv1lambda = Reduce(`+`, deconv1$lambda[1000:10000])/length(1000:10000)
labels = ggplot(data.frame(deconv1$positions), aes(x, -y)) +
  geom_text(color = "grey", label = "\u2B22", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  geom_text(aes(color = factor(deconv1col, levels = c(7,1,4,3,6,2,5))), alpha = 1,
            size = 2, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "Cluster", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A") +
  theme_void() + coord_fixed()
means1 = ggplot(data.frame(deconv1$positions), aes(x, -y)) +
  geom_text(color = "grey", label = "\u2B22", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  geom_text(aes(color = deconv1mu_matrix[deconv1col,1]), alpha = 1,
            size = 2, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "PC1", x = NULL, y = NULL) +
  theme_void() + coord_fixed()
means2 = ggplot(data.frame(deconv1$positions), aes(x, -y)) +
  geom_text(color = "grey", label = "\u2B22", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  geom_text(aes(color = deconv1mu_matrix[deconv1col,2]), alpha = 1,
            size = 2, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "PC2", x = NULL, y = NULL) +
  theme_void() + coord_fixed()
labels + means1 + means2
#simulation 1
simY = matrix(nrow = nrow(deconv1$positions), ncol = 9)
set.seed(100)
for(i in 1:nrow(deconv1$positions)){
  simY[i,] = mvtnorm::rmvnorm(1, mean = deconv1mu_matrix[deconv1col[i],], sigma = solve(deconv1lambda))
}
sim1 = ggplot(data.frame(deconv1$positions), aes(x, -y)) +
  geom_text(color = "grey", label = "\u2B22", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  geom_text(aes(color = simY[,1]), alpha = 1,
            size = 2, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "PC1", alpha = "Proportion", x = NULL, y = NULL) +
  theme_void() + coord_fixed()
sim2 = ggplot(data.frame(deconv1$positions), aes(x, -y)) +
  geom_text(color = "grey", label = "\u2B22", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  geom_text(aes(color = simY[,2]), alpha = 1,
            size = 2, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "PC2", alpha = "Proportion", x = NULL, y = NULL) +
  theme_void() + coord_fixed()
sim3 = ggplot(data.frame(deconv1$positions), aes(x, -y)) +
  geom_text(color = "grey", label = "\u2B22", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  geom_text(aes(color = simY[,3]), alpha = 1,
            size = 2, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "PC3", alpha = "Proportion", x = NULL, y = NULL) +
  theme_void() + coord_fixed()
sim1+sim2+sim3
#aggregation
data_sim = data.frame(simY)
data_sim$j = rep(1:nrow(PCs$components), 7)
data_sim %>%
  group_by(j) %>%
  summarise_all(mean)->
  data_sim_mean
agg1 = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = data_sim_mean$X1), alpha = 1,
            size = 5, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 5, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "PC1", alpha = "Proportion", x = NULL, y = NULL) +
  theme_void() + coord_fixed()
agg2 = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = data_sim_mean$X2), alpha = 1,
            size = 5, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 5, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "PC2", alpha = "Proportion", x = NULL, y = NULL) +
  theme_void() + coord_fixed()
agg3 = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = data_sim_mean$X3), alpha = 1,
            size = 5, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 5, family = "Lucida Sans Unicode", show.legend = F) +
  labs(color = "PC3", alpha = "Proportion", x = NULL, y = NULL) +
  theme_void() + coord_fixed()
agg1+agg2+agg3

# set.seed(7)
# km_sim = kmeans(data_sim_mean[,-1], centers = 7)$cluster
set.seed(7)
mclust_sim = Mclust(data_sim_mean[,-1], G = 7, modelNames = c("EEE"))
clust_sim = cluster(Y= data_sim_mean[,-1], positions = as.matrix(positions), q = 7, init = mclust_sim$classification, 
                    nrep = 10000, gamma = 2, dist = dist)
clust_sim2 = cluster(Y= data_sim_mean[,-1], positions = as.matrix(positions), q = 7, init = mclust_sim$classification, 
                    nrep = 50000, gamma = 3, dist = dist)

simcol = apply(clust_sim$z[seq(1000,10000,10),], 2, Mode)
simcol2 = apply(clust_sim2$z[seq(10000,50000,10),], 2, Mode)

mclust_sim_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color =factor(mclust_sim$classification, levels = c(2,4,5,7,1,3,6))), alpha = 1,
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  guides(col = F)+ 
  scale_color_viridis_d(option = "A") + 
  labs(color = "PC3", alpha = "Proportion", x = NULL, y = NULL) +
  theme_void() + coord_fixed()
clust_sim_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color =factor(simcol, levels = c(2,4,5,7,1,3,6))), alpha = 1,
            size = 6, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  guides(col = F)+ 
  scale_color_viridis_d(option = "A") + 
  labs(color = "PC3", alpha = "Proportion", x = NULL, y = NULL) +
  theme_void() + coord_fixed()
labels+guides(col = F) + mclust_sim_out + clust_sim_out

#re deconvolve
deconv2 = deconvolve(Y = data_sim_mean[,-1], positions = positions, nrep = 50000, gamma = 2, xdist = xdist, ydist = ydist, init = simcol, q = 7)
deconv2 = deconvolve(Y = data_sim_mean[,-1], positions = positions, nrep = 4000, gamma = 2, xdist = xdist, ydist = ydist, init = simcol, q = 7)
deconv2col = apply(deconv2$z[seq(10,4000,10),], 2, Mode)
deconv2colg1 = readRDS("data-raw/deconv2colg1.RDS")
deconv2colg3 = readRDS("data-raw/deconv2colg3.RDS")
g2 = ggplot(data.frame(deconv1$positions), aes(x, -y)) +
  geom_text(color = "grey", label = "\u2B22", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  geom_text(aes(color = factor(deconv2col, levels = c(6,1,3,4,5,7,2))), alpha = 1,
            size = 2, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "Cluster", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A") +
  guides(col = F)+
  theme_void() + coord_fixed()
g1 = ggplot(data.frame(deconv1$positions), aes(x, -y)) +
  geom_text(color = "grey", label = "\u2B22", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  geom_text(aes(color = factor(deconv2colg1, levels = c(6,1,3,4,5,7,2))), alpha = 1,
            size = 2, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "Cluster", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A") +
  guides(col = F)+
  theme_void() + coord_fixed()
g3 = ggplot(data.frame(deconv1$positions), aes(x, -y)) +
  geom_text(color = "grey", label = "\u2B22", size = 6, family = "Lucida Sans Unicode", show.legend = F) +
  geom_text(aes(color = factor(deconv2colg3, levels = c(6,1,3,4,5,7,2))), alpha = 1,
            size = 2, label = "\u2B22", family = "Lucida Sans Unicode") +
  labs(color = "Cluster", x = NULL, y = NULL) +
  scale_color_viridis_d(option = "A") +
  guides(col = F)+
  theme_void() + coord_fixed()
g1+g2+g3

mclust::adjustedRandIndex(deconv1col, rep(mclust_sim$classification,7))
mclust::adjustedRandIndex(deconv1col, rep(simcol,7))
mclust::adjustedRandIndex(deconv1col, g5col)

#Giotto
Sys.setenv(RETICULATE_PYTHON = "C:/Users/Edward Zhao/.conda/envs/Giotto/")
reticulate::py_config()
giotto_data = createGiottoObject(counts(sce),
                                 positions)
giotto_data <- normalizeGiotto(gobject = giotto_data, scalefactor = 6000, verbose = T)
spatPlot(gobject = giotto_data)
giotto_data <- calculateHVG(gobject = giotto_data, method = 'cov_loess', difference_in_cov = 0.1,
                        save_param = list(save_name = '3_HVGplot', base_height = 5, base_width = 5))
gene_metadata = fDataDT(giotto_data)
featgenes = gene_metadata[hvg == 'yes']$gene_ID
giotto_data <- Giotto::runPCA(gobject = giotto_data, genes_to_use = featgenes, scale_unit = F)
plotPCA(gobject = giotto_data)
giotto_data = createSpatialNetwork(gobject = giotto_data, maximum_distance_delaunay = 9)
spatPlot(gobject=giotto_data, show_network = T, network_color = "blue",
         spatial_network_name = 'Delaunay_network')
rank_spatialgenes = binSpect(giotto_data, bin_method = 'rank')
results_folder = "C:/Users/Edward Zhao/Google Drive/Gottardo/spatial"
hmrf_folder = paste0(results_folder,'/','HMRF/')

test = doHMRF(gobject = giotto_data,
              dim_reduction_name = "pca", dim_reduction_to_use = "pca",
              k = 7, betas = c(7,2,3),
              output_folder = paste0(hmrf_folder, '/', 'PCA/PCA_k9'))
viewHMRFresults2D(giotto_data,
                  test,
                  k = 7,
                  betas_to_view = 7)
VC_test = addHMRF(gobject = giotto_data,
                  HMRFoutput = test,
                  k = 7, betas_to_add = c(7),
                  hmrf_name = 'HMRF_7')
