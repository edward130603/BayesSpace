library(SingleCellExperiment)

test_that("clustering matches", {
  sce <- readRDS(system.file("testdata/maynard_151673_subset.rds", package="BayesSpace"))
  
  positions = cbind(sce$imagecol, sce$imagerow)
  colnames(positions) = c("x", "y") 
  
  set.seed(941)
  out <- cluster(Y = metadata(sce)$PCs, 
                 positions = as.matrix(positions), 
                 q = 7, 
                 z0 = colData(sce)$init_km7, 
                 nrep = 1000, 
                 gamma = 1.5, 
                 neighborhood.radius = metadata(sce)$dist,
                 model="normal",
                 precision="equal")
  
  labels <- apply(out$z[900:1000,], 2, Mode)
  
  expect_true(all(labels == colData(sce)$labels_normal_equal))
})
