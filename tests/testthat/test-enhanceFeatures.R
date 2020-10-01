## ========================
## Setup
## ========================
set.seed(149)
sce <- exampleSCE(nrow=8, ncol=12, n_genes=20)
enhanced <- exampleSCE(nrow=21, ncol=32, n_genes=20)

X.ref <- reducedDim(sce, "PCA")
X.enhanced <- reducedDim(enhanced, "PCA")
Y.ref <- logcounts(sce)
gene <- "gene_4"  # example gene for expression prediction

## Example feature matrix (cell type proportions over 5 types)
xs <- matrix(runif(5 * ncol(sce)), ncol=ncol(sce))
props <- sweep(xs, 2, colSums(xs), "/")
colnames(props) <- colnames(sce)
rownames(props) <- paste0("type_", seq_len(nrow(props)))
altExp(sce, "cell_type") <- SummarizedExperiment(assays=list("cell_type"=props))

## ========================
## Test helpers
## ========================

test_that("genes are predicted with linear model", {
  fit <- lm(Y.ref[gene, ] ~ ., data=as.data.frame(X.ref))
  new <- predict(fit, newdata=as.data.frame(X.enhanced))
  Y.enhanced <- .lm_enhance(as.data.frame(X.ref), as.data.frame(X.enhanced), Y.ref, c(gene))
  
  expect_equal(new, Y.enhanced[gene, ], check.names=FALSE)
})

test_that("genes are predicted with xgboost", {
  gene <- "gene_4"
  fit <- xgboost(data=X.ref, label=Y.ref[gene, ], objective="reg:squarederror",
                 max_depth=2, eta=0.03, nrounds=100, nthread=1, verbose=FALSE)
  new <- predict(fit, newdata=X.enhanced)
  
  Y.enhanced <- .xgboost_enhance(X.ref, X.enhanced, Y.ref, c(gene), nrounds = 100)
  
  expect_equal(new, Y.enhanced[gene, ], check.names=FALSE)
  expect_true("diagnostic" %in% names(attributes(Y.enhanced)))
  expect_true("rmse" %in% names(attr(Y.enhanced, "diagnostic")))
})

test_that("genes are predicted with tuned xgboost", {
  gene <- "gene_4"
  set.seed(100)
  
  train.index <- sample(1:ncol(Y.ref), 64)
  data.train <- xgb.DMatrix(data=X.ref[train.index, ],
                            label=Y.ref[gene, train.index])
  data.test  <- xgb.DMatrix(data=X.ref[-train.index, ],
                            label=Y.ref[gene, -train.index])
  watchlist <- list(train=data.train, test=data.test)
  
  fit.train <- xgb.train(data=data.train, max_depth=2, watchlist=watchlist,
                         eta=0.03, nrounds=500, objective="reg:squarederror",
                         verbose=FALSE)
  
  nrounds <- which.min(fit.train$evaluation_log$test_rmse)
  
  fit <- xgboost(data=X.ref, label=Y.ref[gene, ], objective="reg:squarederror",
                 max_depth=2, eta=0.03, nrounds=nrounds, nthread=1, 
                 verbose=FALSE)
  new <- predict(fit, newdata=X.enhanced)
  
  set.seed(100)
  Y.enhanced <- .xgboost_enhance(X.ref, X.enhanced, Y.ref, c(gene),
                                 nrounds=0, train.n=64)
  
  expect_equal(new, Y.enhanced[gene, ], check.names=FALSE)
  expect_true("diagnostic" %in% names(attributes(Y.enhanced)))
  expect_true("rmse" %in% names(attr(Y.enhanced, "diagnostic")))
})

test_that("proportions are predicted with Dirichlet model", { 
  features <- DR_data(t(as.data.frame(props)))
  fit <- DirichReg(features ~ ., data = as.data.frame(X.ref))
  new <- t(predict(fit, newdata = as.data.frame(X.enhanced)))
  enhanced.props <- .dirichlet_enhance(as.data.frame(X.ref), as.data.frame(X.enhanced), props)
  
  expect_equal(as.vector(new), as.vector(enhanced.props), check.names=FALSE)
  expect_equal(colSums(enhanced.props), rep(1, ncol(enhanced)), check.names=FALSE)
})

## ========================
## Test enhanceFeatures
## ========================
enhanced <- enhanceFeatures(enhanced, sce, assay.type="logcounts")

test_that("rownames are required", {
  rownames(props) <- NULL
  altExp(sce, "cell_type") <- SummarizedExperiment(assays=list("cell_type"=props))
  
  expect_error(enhanced <- enhanceFeatures(enhanced, sce, altExp.type="cell_type"))
  expect_error(enhanced.features <- enhanceFeatures(enhanced, sce, feature.matrix=props))
})

test_that("enhanced features are added to assay", {
  expect_true("logcounts" %in% assayNames(enhanced))
  expect_equal(nrow(logcounts(enhanced)), nrow(sce))
  expect_equal(ncol(logcounts(enhanced)), ncol(sce) * 7)
})

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
