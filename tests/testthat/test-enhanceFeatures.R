set.seed(149)
sce <- exampleSCE(nrow=8, ncol=12)
enhanced <- exampleSCE(nrow=21, ncol=32)
enhanced.lm <- enhanceFeatures(enhanced, sce, assay.type="logcounts")

test_that("enhanced features are added", {
  expect_true("logcounts" %in% assayNames(enhanced.lm))
  expect_equal(nrow(logcounts(enhanced.lm)), nrow(sce))
  expect_equal(ncol(logcounts(enhanced.lm)), ncol(sce) * 7)
})

test_that("genes are predicted with linear model", {
  gene <- "gene_4"
  spot_PCs <- as.data.frame(reducedDim(sce, "PCA"))
  subspot_PCs <- as.data.frame(reducedDim(enhanced, "PCA"))
  
  fit <- lm(logcounts(sce)[gene, ] ~ ., data=spot_PCs)
  new <- predict(fit, newdata=subspot_PCs)
  
  expect_equal(new, logcounts(enhanced.lm)[gene, ])
})
