library(SingleCellExperiment)

test_that("clustering matches", {
  # sce <- readRDS("../../data/maynard_151673.rds")
  sce <- readRDS("../../data/maynard_151673_subset.rds")
  
  positions = cbind(sce$imagecol, sce$imagerow)
  colnames(positions) = c("x", "y") 
  
  set.seed(941)
  out <- cluster(Y = metadata(sce)$PCs, 
                 positions = as.matrix(positions), 
                 q = 7, 
                 init = colData(sce)$init_km7, 
                 nrep = 1000, 
                 gamma = 1.5, 
                 dist = metadata(sce)$dist,
                 model="normal",
                 precision="equal")
  
  labels <- apply(out$z[900:1000,], 2, Mode)
  
  expect_true(all(labels == colData(sce)$labels_normal_equal))
})
