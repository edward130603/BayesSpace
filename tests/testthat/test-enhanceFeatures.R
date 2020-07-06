## Example spot and subspot SCEs
set.seed(149)
sce <- exampleSCE(nrow=8, ncol=12, n_genes=20)
enhanced <- exampleSCE(nrow=21, ncol=32, n_genes=20)
enhanced <- enhanceFeatures(enhanced, sce, assay.type="logcounts")

test_that("enhanced features are added to assay", {
  expect_true("logcounts" %in% assayNames(enhanced))
  expect_equal(nrow(logcounts(enhanced)), nrow(sce))
  expect_equal(ncol(logcounts(enhanced)), ncol(sce) * 7)
})

test_that("genes are predicted with linear model", {
  gene <- "gene_4"
  spot_PCs <- as.data.frame(reducedDim(sce, "PCA"))
  subspot_PCs <- as.data.frame(reducedDim(enhanced, "PCA"))
  
  fit <- lm(logcounts(sce)[gene, ] ~ ., data=spot_PCs)
  new <- predict(fit, newdata=subspot_PCs)
  
  expect_equal(new, logcounts(enhanced)[gene, ])
})

## Example feature matrix (cell type proportions over 5 types)
set.seed(149)
xs <- matrix(runif(5 * ncol(sce)), ncol=ncol(sce))
props <- sweep(xs, 2, colSums(xs), "/")
colnames(props) <- colnames(sce)

test_that("rownames are required", {
  altExp(sce, "cell_type") <- SummarizedExperiment(assays=list("cell_type"=props))
  
  expect_error(enhanced <- enhanceFeatures(enhanced, sce, altExp.type="cell_type"))
  expect_error(enhanced.features <- enhanceFeatures(enhanced, sce, feature.matrix=props))
})

## Add cell type proportions to sce as altExp
rownames(props) <- paste0("type_", seq_len(nrow(props)))
altExp(sce, "cell_type") <- SummarizedExperiment(assays=list("cell_type"=props))

test_that("enhanced features are added to altExp", {
  enhanced <- enhanceFeatures(enhanced, sce, altExp.type="cell_type")
  
  expect_true("cell_type" %in% altExpNames(enhanced))
  expect_equal(nrow(altExp(enhanced, "cell_type")), nrow(props))
  expect_equal(ncol(altExp(enhanced, "cell_type")), ncol(sce) * 7)
})

test_that("enhanced features are returned as matrix", {
  enhanced.features <- enhanceFeatures(enhanced, sce, feature.matrix=props)
  
  expect_true(is.matrix(enhanced.features))
  expect_equal(nrow(enhanced.features), nrow(props))
  expect_equal(ncol(enhanced.features), ncol(sce) * 7)
})