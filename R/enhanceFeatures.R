#' Predict log-normalized expression vectors from deconvolved PCs using linear 
#' regression.
#' 
#' @param sce SingleCellExperiment object with original PCs and logcounts
#' @param newdata Deconvolved PCs (rows are subspots, columns are PCs)
#' @param dimred Name of dimension reduction to use
#' @param genes List of genes to predict expression for
#' @param components Number of reduced dimensions (PCs) to use
#' 
#' @return Returns list with names:
#'    * \code{expression} - Deconvolved expression values. 
#'      (Rows are genes, columns are subspots)
#'    * \code{r2} - Percent of variation in original gene expression 
#'      explained by PCs
#'
#' @name enhanceFeatures
NULL

#' @importFrom assertthat assert_that
.enhance_features <- function(X.enhanced, X.ref, Y.ref, 
    features = rownames(Y.ref), model = c("lm", "xgboost")) {

    assert_that(ncol(X.enhanced) == ncol(X.ref))
    assert_that(ncol(Y.ref) == nrow(X.ref))
    model <- match.arg(model)
    
    X.ref <- as.data.frame(X.ref)
    X.enhanced <- as.data.frame(X.enhanced)
    
    if (!all(colnames(X.enhanced) == colnames(X.ref))) {
        warning("colnames of reducedDim and X.enhanced do not match.")
        warning("Setting X.enhanced colnames to match reducedDim.")
        colnames(X.enhanced) <- colnames(X.ref)
    }

    if (model == "lm") {
        .lm_enhance(X.ref, X.enhanced, Y.ref, features)
    }
}

#' @importFrom stats lm predict
.lm_enhance <- function(X.ref, X.enhanced, Y.ref, features) {
    r.squared <- numeric(length(features))
    names(r.squared) <- features
    
    Y.enhanced <- matrix(nrow=length(features), ncol=nrow(X.enhanced))
    rownames(Y.enhanced) <- features
    colnames(Y.enhanced) <- rownames(X.enhanced)
    
    for (feature in features) {
        fit <- lm(Y.ref[feature, ] ~ ., data=X.ref)
        r.squared[feature] <- summary(fit)$r.squared
        Y.enhanced[feature, ] <- predict(fit, newdata=X.enhanced)
    }
    
    attr(Y.enhanced, "r.squared") <- r.squared
    Y.enhanced
}

## TODO: store deconvolved data in sce assays or altExps or reducedDims

#' @export
#' @rdname enhanceFeatures
setGeneric("enhanceFeatures", function(enhanced, reference, ...) standardGeneric("enhanceFeatures"))

#' @export
#' @rdname enhanceFeatures
setMethod("enhanceFeatures", c("ANY", "ANY"), .enhance_features)

#' @export
#' @rdname enhanceFeatures
setMethod("enhanceFeatures", 
    c("SingleCellExperiment", "ANY"), 
    function (enhanced, reference, ..., use.dimred = "PCA")
{
        
})

#' @export
#' @rdname enhanceFeatures
setMethod("enhanceFeatures", 
    c("SingleCellExperiment", "SingleCellExperiment"), 
    function (enhanced, reference, ..., use.dimred = "PCA") {
    
}
