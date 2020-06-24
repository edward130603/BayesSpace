library(SingleCellExperiment)

test_that("clustering matches", {
  sce <- readRDS(system.file("testdata/maynard_151673_subset2.rds", package="BayesSpace"))
  
  positions = cbind(sce$imagecol, sce$imagerow)
  colnames(positions) = c("x", "y") 
  
  set.seed(149)
  out <- cluster(Y = reducedDim(sce, "PCA"), 
                 positions = as.matrix(positions), 
                 q = 7, 
                 z0 = sce$km_init, 
                 nrep = 1000, 
                 gamma = 1.5, 
                 neighborhood.radius = metadata(sce)$dist,
                 model="normal",
                 precision="equal")
  
  labels <- apply(out$z[900:1000,], 2, Mode)
  
  expect_true("truth" %in% names(colData(sce)))
  expect_true(all(labels == sce$truth))
})

