library(SingleCellExperiment)

test_that("SCE deconvolution matches", {
  sce <- readRDS(system.file("testdata/maynard_151673_subset2.rds", package="BayesSpace"))
  truth <- read.csv(system.file("testdata/maynard_151673_subset2.deconv_truth.csv", package="BayesSpace"))
  
  positions <- cbind(sce$imagecol, sce$imagerow) 
  colnames(positions) <- c("x", "y") 
  
  set.seed(149)
  enhanced <- spatialEnhance(sce, 7, use.dimred="PCA", init=sce$truth, positions=positions)

  expect_true(all(enhanced$spatial.cluster == truth$x))
})
