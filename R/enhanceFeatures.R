#' Predict feature vectors from enhanced PCs.
#' 
#' @param sce.enhanced SingleCellExperiment object with enhanced PCs.
#' @param sce.ref SingleCellExperiment object with original PCs and expression.
#' @param feature_names List of genes/features to predict expression/values for.
#' @param model Model used to predict enhanced values.
#' @param use.dimred Name of dimension reduction to use.
#' @param assay.type Expression matrix in \code{assays(sce.ref)} to predict.
#' @param altExp.type Expression matrix in \code{altExps(sce.ref)} to predict.
#'   Overrides \code{assay.type} if specified.
#' @param feature.matrix Expression/feature matrix to predict, if not directly
#'   attached to \code{sce.ref}. Must have columns corresponding to the spots in
#'   \code{sce.ref}. Overrides \code{assay.type} and \code{altExp.type} if
#'   specified.
#' @param nrounds Nonnegative integer to set the \code{nrounds} parameter
#'   (max number of boosting iterations) for xgboost. \code{nrounds = 100}
#'   works reasonably well in most cases. If \code{nrounds} is set to 0, the
#'   parameter will be tuned using a train-test split. We recommend tuning
#'   \code{nrounds} for improved feature prediction, but note this will increase
#'   runtime.
#' @param train.n Number of spots to use in the training dataset for tuning
#'   nrounds. By default, 2/3 the total number of spots are used.
#' 
#' @return If \code{assay.type} or \code{altExp.type} are specified, the
#'   enhanced features are stored in the corresponding slot of
#'   \code{sce.enhanced} and the modified SingleCellExperiment object is
#'   returned.
#'   
#'   If \code{feature.matrix} is specified, or if a subset of features are
#'   requested, the enhanced features are returned directly as a matrix.
#' 
#' @details 
#' Enhanced features are computed by fitting a predictive model to a
#' low-dimensional representation of the original expression vectors. By
#' default, a linear model is fit for each gene using the top 15 principal
#' components from each spot, i.e. \code{lm(gene ~ PCs)}, and the fitted model
#' is used to predict the enhanced expression for each gene from the subspots'
#' principal components.
#' 
#' Diagnostic measures, such as RMSE for \code{xgboost} or R.squared for linear
#' regression, are added to the `rowData` of the enhanced experiment if the
#' features are an assay of the original experiment. Otherwise they are stored
#' as an attribute of the returned matrix/altExp.
#' 
#' Note that feature matrices will be returned and are expected to be input as
#' \eqn{p \times n} matrices of \eqn{p}-dimensional feature vectors over the
#' \eqn{n} spots.
#' 
#' @examples
#' set.seed(149)
#' sce <- exampleSCE()
#' sce <- spatialCluster(sce, 7, nrep=200, burn.in=20)
#' enhanced <- spatialEnhance(sce, 7, init=sce$spatial.cluster, nrep=200, burn.in=20)
#' enhanced <- enhanceFeatures(enhanced, sce, assay.type="logcounts")
#'
#' @name enhanceFeatures
NULL

#' @importFrom assertthat assert_that
.enhance_features <- function(X.enhanced, X.ref, Y.ref, 
    feature_names = rownames(Y.ref), model = c("xgboost", "dirichlet", "lm"), 
    nrounds, train.n) {

    assert_that(ncol(X.enhanced) == ncol(X.ref))
    assert_that(ncol(Y.ref) == nrow(X.ref))
    model <- match.arg(model)
    
    if (model %in% c("lm", "dirichlet")) {
        X.ref <- as.data.frame(X.ref)
        X.enhanced <- as.data.frame(X.enhanced)
    }
    
    if (!all(colnames(X.enhanced) == colnames(X.ref))) {
        warning("colnames of reducedDim and X.enhanced do not match.\n",
                "  Setting X.enhanced colnames to match reducedDim.")
        colnames(X.enhanced) <- colnames(X.ref)
    }

    if (model == "lm") {
        .lm_enhance(X.ref, X.enhanced, Y.ref, feature_names)
    } else if (model == "dirichlet") {
        .dirichlet_enhance(X.ref, X.enhanced, Y.ref)
    } else if (model == "xgboost") {
        .xgboost_enhance(X.ref, X.enhanced, Y.ref, feature_names, nrounds, train.n)
    }
}

#' @importFrom stats lm predict
.lm_enhance <- function(X.ref, X.enhanced, Y.ref, feature_names) {
    r.squared <- numeric(length(feature_names))
    names(r.squared) <- feature_names
    
    Y.enhanced <- matrix(nrow=length(feature_names), ncol=nrow(X.enhanced))
    rownames(Y.enhanced) <- feature_names
    colnames(Y.enhanced) <- rownames(X.enhanced)
    
    for (feature in feature_names) {
        fit <- lm(Y.ref[feature, ] ~ ., data=X.ref)
        r.squared[feature] <- summary(fit)$r.squared
        Y.enhanced[feature, ] <- predict(fit, newdata=X.enhanced)
    }
    
    diagnostic <- list("r.squared"=r.squared)
    attr(Y.enhanced, "diagnostic") <- diagnostic
    Y.enhanced
}

#' @importFrom DirichletReg DR_data DirichReg
#' @importFrom stats predict
.dirichlet_enhance <- function(X.ref, X.enhanced, Y.ref) {
    features <- DR_data(t(Y.ref))
    
    fit <- DirichReg(features ~ ., data = X.ref)
    Y.enhanced <- t(predict(fit, newdata = X.enhanced))
    
    rownames(Y.enhanced) <- rownames(Y.ref)
    colnames(Y.enhanced) <- rownames(X.enhanced)
    
    attr(Y.enhanced, "diagnostic") <- list()
    Y.enhanced
}

#' @importFrom xgboost xgboost xgb.DMatrix xgb.train
.xgboost_enhance <- function(X.ref, X.enhanced, Y.ref, feature_names, 
                             nrounds, train.n) {

    Y.enhanced <- matrix(nrow=length(feature_names), ncol=nrow(X.enhanced))
    rownames(Y.enhanced) <- feature_names
    colnames(Y.enhanced) <- rownames(X.enhanced)
    
    rmse <- numeric(length(feature_names))
    names(rmse) <- feature_names
    
    if (nrounds == 0){
        train.index <- sample(seq_len(ncol(Y.ref)), train.n)
    }
    default.nrounds <- nrounds
    
    ## TODO: think about extracting tuning function, using apply instead of loop
    for (feature in feature_names) {
        nrounds <- default.nrounds
        if (nrounds == 0){
            data.train <- xgb.DMatrix(data=X.ref[train.index, ],  
                                      label=Y.ref[feature, train.index])
            data.test  <- xgb.DMatrix(data=X.ref[-train.index, ], 
                                      label=Y.ref[feature, -train.index])
            watchlist <- list(train=data.train, test=data.test)
            
            fit.train <- xgb.train(data=data.train, max_depth=2, 
                                   watchlist=watchlist, eta=0.03, nrounds=500,
                                   objective="reg:squarederror",
                                   verbose=FALSE)
            nrounds <- which.min(fit.train$evaluation_log$test_rmse)
        }
        
        fit <- xgboost(data=X.ref, label=Y.ref[feature, ], 
            objective="reg:squarederror", max_depth=2, eta=0.03,
            nrounds=nrounds, nthread=1, verbose=FALSE)
        
        Y.enhanced[feature, ] <- predict(fit, X.enhanced)
        rmse[feature] <- fit$evaluation_log$train_rmse[nrounds]
    }
    
    diagnostic <- list("rmse"=rmse)
    attr(Y.enhanced, "diagnostic") <- diagnostic
    
    Y.enhanced
}

#' @export
#' @importFrom SingleCellExperiment reducedDim altExp altExp<-
#' @importFrom SummarizedExperiment assay assay<- SummarizedExperiment rowData<-
#' @rdname enhanceFeatures
enhanceFeatures <- function(sce.enhanced, sce.ref, feature_names = NULL,
    model=c("xgboost", "dirichlet", "lm"), use.dimred = "PCA",
    assay.type="logcounts", altExp.type = NULL, feature.matrix = NULL,
    nrounds = 0, train.n = round(ncol(sce.ref)*2/3)) {
    
    X.enhanced <- reducedDim(sce.enhanced, use.dimred)
    X.ref <- reducedDim(sce.ref, use.dimred)
    
    ## If user specified clustering with fewer PCs than in the ref dataset,
    ## there will be fewer PCs in the enhanced SCE. Only use these.
    d <- min(ncol(X.enhanced), ncol(X.ref))
    X.enhanced <- X.enhanced[, seq_len(d)]
    X.ref <- X.ref[, seq_len(d)]
    
    if (!is.null(feature.matrix)) {
        Y.ref <- feature.matrix
    } else if (!is.null(altExp.type)) {
        Y.ref <- assay(altExp(sce.ref, altExp.type), altExp.type)
    } else {
        Y.ref <- assay(sce.ref, assay.type)
    }
    
    msg <- "Spot features must have assigned rownames."
    assert_that(!is.null(rownames(Y.ref)), msg=msg)
    
    if (is.null(feature_names)) {
        feature_names <- rownames(Y.ref)
    } else {
        feature_names <- intersect(feature_names, rownames(Y.ref))
        skipped_features <- setdiff(feature_names, rownames(Y.ref))
        if (length(skipped_features) > 0) {
            message("Skipping ", length(skipped_features), " features not in sce.ref.")
        }
    }
    
    Y.enhanced <- .enhance_features(X.enhanced, X.ref, Y.ref, feature_names, 
                                    model, nrounds, train.n)

    ## Clip negative predicted expression
    Y.enhanced <- pmax(Y.enhanced, 0)
    
    ## TODO: add option to specify destination of enhanced features.
    ## For now, return in same form as input
    if (!is.null(feature.matrix)) {
        return(Y.enhanced)
    } else if (!is.null(altExp.type)) {
        Y.enhanced <- SummarizedExperiment(assays=list(altExp.type=Y.enhanced))
        altExp(sce.enhanced, altExp.type) <- Y.enhanced
    } else {
        diagnostic <- attr(Y.enhanced, "diagnostic")
        
        ## If we only enhanced a subset of features, need to add NA vectors for
        ## the remaining features so the number of rows within the SCE remains
        ## consistent
        if (length(feature_names) != nrow(Y.ref)) {
            Y.full <- matrix(data=NA, nrow=nrow(Y.ref), ncol=ncol(sce.enhanced))
            rownames(Y.full) <- rownames(Y.ref)
            colnames(Y.full) <- colnames(Y.enhanced)
            Y.full[feature_names, ] <- Y.enhanced
            assay(sce.enhanced, assay.type) <- Y.full
            
            ## Fill in the diagnostic/error values as above
            for (name in names(diagnostic)) {
                diagnostic.full <- rep(NA, nrow(sce.ref))
                names(diagnostic.full) <- rownames(sce.ref)
                diagnostic.full[feature_names] <- diagnostic[[name]]
                col.name <- sprintf("enhanceFeatures.%s", name)
                rowData(sce.enhanced)[[col.name]] <- diagnostic.full
            }
        } else {
            assay(sce.enhanced, assay.type) <- Y.enhanced
            for (name in names(diagnostic)) {
                col.name <- sprintf("enhanceFeatures.%s", name)
                rowData(sce.enhanced)[[col.name]] <- diagnostic[[name]]
            }
        }
    }
    
    return(sce.enhanced)
}
