#' Predict log-normalized expression vectors from deconvolved PCs using linear 
#' regression.
#' 
#' @param sce.enhanced SingleCellExperiment object with enhanced PCs.
#' @param sce.ref SingleCellExperiment object with original PCs and expression..
#' @param use.dimred Name of dimension reduction to use.
#' @param assay.type Expression matrix in \code{assays(sce.ref)} to predict.
#' @param altExp.type Expression matrix in \code{altExps(sce.ref)} to predict.
#'   Overrides \code{assay.type} if specified.
#' @param feature.matrix Expression/feature matrix to predict, if not directly
#'   attached to \code{sce.ref}. Must have columns corresponding to the spots in
#'   \code{sce.ref}. Overrides \code{assay.type} and \code{altExp.type} if
#'   specified.
#' @param features List of genes/features to predict expression/values for.
#' @param model Model used to predict enhanced values.
#' 
#' @return If \code{assay.type} or \code{altExp.type} are specified, the
#'   enhanced features are stored in the corresponding slot of
#'   \code{sce.enhanced} and the modified SingleCellExperiment object is
#'   returned.
#'   
#'   If \code{feature.matrix} is specified, the enhanced features are returned
#'   directly as a matrix.
#' 
#' @details 
#' Enhanced features are computed by fitting a predictive model to a
#' low-dimensional representation of the original expression vectors. By
#' default, a linear model is fit for each gene using the top 15 principal
#' compoenents from each spot, i.e. \code{lm(gene ~ PCs)}, and the fitted model
#' is used to predict the enhanced expression for each gene from the subspots'
#' principal components.
#' 
#' @examples
#' set.seed(149)
#' sce <- exampleSCE()
#' sce <- spatialCluster(sce, 7)
#' enhanced <- spatialEnhance(sce, 7, init=sce$spatial.cluster)
#' enhanced <- enhanceFeatures(enhanced, sce, assay.type="logcounts")
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

#' @export
#' @importFrom SingleCellExperiment reducedDim altExp altExp<-
#' @importFrom SummarizedExperiment assay assay<- SummarizedExperiment
#' @rdname enhanceFeatures
enhanceFeatures <- function(sce.enhanced, sce.ref, use.dimred = "PCA",
    assay.type="logcounts", altExp.type = NULL, feature.matrix = NULL, ...) {
    
    X.enhanced <- reducedDim(sce.enhanced, use.dimred)
    X.ref <- reducedDim(sce.ref, use.dimred)
    
    if (!is.null(feature.matrix)) {
        Y.ref <- feature.matrix
    } else if (!is.null(altExp.type)) {
        Y.ref <- assay(altExp(sce.ref, altExp.type), altExp.type)
    } else {
        Y.ref <- assay(sce.ref, assay.type)
    }
    
    msg <- "Spot features must have assigned rownames."
    assert_that(!is.null(rownames(Y.ref)), msg=msg)
    Y.enhanced <- .enhance_features(X.enhanced, X.ref, Y.ref, ...)
    
    ## TODO: add option to specify destination of enhanced features.
    ## For now, return in same form as input
    if (!is.null(feature.matrix)) {
        return(Y.enhanced)
    } else if (!is.null(altExp.type)) {
        Y.enhanced <- SummarizedExperiment(assays=list(altExp.type=Y.enhanced))
        altExp(sce.enhanced, altExp.type) <- Y.enhanced
    } else {
        assay(sce.enhanced, assay.type) <- Y.enhanced
    }
    
    return(sce.enhanced)
}