library(SingleCellExperiment)

PCs <- read.csv(system.file("testdata/maynard_151673_subset.PCs.csv.gz", package="BayesSpace"), row.names=1)
cdata <- read.csv(system.file("testdata/maynard_151673_subset.colData.csv.gz", package="BayesSpace"), row.names=1)
meta <- readRDS(system.file("testdata/maynard_151673_subset.metadata.rds", package="BayesSpace"))
positions <- cdata[, c("imagecol", "imagerow")]

test_that("Refactored clustering matches", {
  skip("Set to match results from iterations 900:1000")
  set.seed(149)
  out <- cluster(Y = PCs, 
                 positions = as.matrix(positions), 
                 q = 7,
                 init = cdata$km.init, 
                 nrep = 1000, 
                 gamma = 1.5, 
                 radius = meta$dist,
                 model="normal",
                 precision="equal")
  
  labels <- apply(out$z[900:1000,], 2, BayesSpace:::Mode)
  
  expect_true("spatial.cluster" %in% names(cdata))
  expect_true(all(labels == cdata$spatial.cluster))
})

test_that("SCE clustering matches", {
  skip("Set to match results from iterations 900:1000")
  sce <- SingleCellExperiment(assays=list(), reducedDims=list("PCA"=PCs), colData=cdata[, -6])
  
  set.seed(149)
  sce <- spatialCluster(sce, 7, use.dimred="PCA", init=sce$km.init, gamma=1.5)

  expect_true(all(sce$spatial.cluster == cdata$spatial.cluster))
})

test_that("Refactored deconvolution (SCE) matches", {
  skip("Set to match results from iterations 900:1000")
  truth <- read.csv(system.file("testdata/maynard_151673_subset.enhance_truth.csv.gz", package="BayesSpace"))
  sce <- SingleCellExperiment(assays=list(), reducedDims=list("PCA"=PCs), colData=cdata)

  set.seed(149)
  enhanced <- spatialEnhance(sce, 7, use.dimred="PCA", init=sce$spatial.cluster)

  expect_true(all(enhanced$spatial.cluster == truth$x))
})

