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
#' @importFrom stats lm predict
#' @export
predictExpression <- function(sce, newdata, dimred = "PCA",
    genes = rownames(sce), components = ncol(newdata)) {
    
    ## TODO: swap args so the enhanced sce we're predicting on is `sce`
    actual_data <- data.frame(reducedDim(sce, dimred))[, seq_len(components)]
    newdata <- as.data.frame(newdata)
    if (ncol(actual_data) != ncol(newdata)) {
        stop("number of components do not match")
    }
    if (!all(colnames(newdata) == colnames(actual_data))) {
        warning("colnames of reducedDim and newdata do not match.")
        warning("Setting newdata colnames to match reducedDim.")
        colnames(newdata) <- colnames(actual_data)
    }
    rsquared <- numeric(length(genes))
    names(rsquared) <- genes
    deconv_expression <- matrix(nrow=length(genes), ncol=nrow(newdata))
    rownames(deconv_expression) <- genes
    colnames(deconv_expression) <- rownames(newdata)
    for (gene in genes) {
        train <- lm(logcounts(sce)[gene, ] ~ ., data=actual_data)
        rsquared[gene] <- summary(train)$r.squared
        deconv_expression[gene, ] <- predict(train, newdata=newdata)
    }
    
    ## TODO: store deconvolved data in sce assays or altExps or reducedDims
    list(expression=deconv_expression, r2=rsquared)
}
