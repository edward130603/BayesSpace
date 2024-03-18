#' Preprocess a spatial dataset for BayesSpace
#' 
#' Adds metadata required for downstream analyses, and (optionally) performs PCA
#' on log-normalized expression of top HVGs.
#' 
#' @param sce SingleCellExperiment to preprocess
#' @param platform Spatial sequencing platform. Used to determine spot layout
#'   and neighborhood structure (Visium = hex, VisiumHD = square, ST = square).
#' @param n.PCs Number of principal components to compute. We suggest using the
#'   top 15 PCs in most cases.
#' @param n.HVGs Number of highly variable genes to run PCA upon.
#' @param skip.PCA Skip PCA (if dimensionality reduction was previously
#'   computed.)
#' @param log.normalize Whether to log-normalize the input data with scater. May
#'   be omitted if log-normalization previously computed.
#' @param assay.type Name of assay in \code{sce} containing normalized counts.
#'   Leave as "logcounts" unless you explicitly pre-computed a different
#'   normalization and added it to \code{sce} under another assay. Note that we
#'   do not recommend running BayesSpace on PCs computed from raw counts.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying which
#'   algorithm should be used to perform the PCA. By default, an exact PCA is
#'   performed, as current spatial datasets are generally small (<10,000 spots).
#'   To perform a faster approximate PCA, please specify
#'   \code{FastAutoParam()} and set a random seed to ensure
#'   reproducibility.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#'   to compute PCA in parallel or not (default to \code{SerialParam()}).
#'   To perform faster PCA decomposition, please specify \code{SnowParam()} or
#'   \code{MulticoreParam()}.
#' 
#' @return SingleCellExperiment with PCA and BayesSpace metadata
#' 
#' @examples
#' sce <- exampleSCE()
#' sce <- spatialPreprocess(sce)
#' 
#' @export
#' @importFrom scater logNormCounts runPCA
#' @importFrom scran modelGeneVar getTopHVGs
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom BiocSingular ExactParam
#' @importFrom BiocParallel SerialParam
spatialPreprocess <- function(sce, platform=c("Visium", "VisiumHD", "ST"),
                              n.PCs=15, n.HVGs=2000, skip.PCA=FALSE,
                              log.normalize=TRUE, assay.type="logcounts",
                              BSPARAM=ExactParam(), BPPARAM = SerialParam()) {
    
    ## Set BayesSpace metadata
    metadata(sce)$BayesSpace.data <- .bsData(sce, name = NULL, default = list())
    metadata(sce)$BayesSpace.data$platform <- match.arg(platform)
    metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
    # metadata(sce)$BayesSpace.data$use_dimred <- use.dimred
    # metadata(sce)$BayesSpace.data$d <- n.PCs

    ## Run PCA on HVGs, log-normalizing if necessary
    if (!skip.PCA) {
        if (log.normalize)
            sce <- logNormCounts(sce)
   
        dec <- modelGeneVar(sce, assay.type=assay.type)
        top <- getTopHVGs(dec, n=n.HVGs)
        sce <- runPCA(sce, subset_row=top, ncomponents=n.PCs, 
                      exprs_values=assay.type, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
        rowData(sce)[["is.HVG"]] <- (rownames(sce) %in% top)
    }

    sce
}
