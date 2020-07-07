#' Run SPOTlight cell type deconvolution on a spatial SCE
#' 
#' desc
#' 
#' @param spatial SingleCellExperiment or expression matrix from the spatial
#'   transcriptomic experiment.
#' @param scRNA SingleCellExperiment or expression matrix from a reference
#'   single-cell RNA-seq experiment.
#' @param cell_types Name of column in \code{colData(scRNA)} containing cell
#'   type labels, or a list of labels.
#' @param sctransform.spatial Run sctransform to normalize spatial expression.
#'   Recommended by SPOTlight developers, but set to false if providing
#'   log-normalized counts. (Will reduce runtime.)
#' @param sctransform.scRNA Run sctransform to normalize scRNA expression.
#'   Recommended by SPOTlight developers, but set to false if providing
#'   normalized counts. (Published reference datasets may only provide
#'   normalized counts.)
#' 
#' @export
#' @name runSPOTlight
NULL

.run_spotlight <- function(spatial, scRNA,
    sctransform.spatial = TRUE, sctransform.scRNA = FALSE, ...) {
    
    .check_spotlight_environment()
    
    # TODO: add cell types, filter out excluded types
    
    sc <- .spotlight_preprocess(scRNA, "RNA", sctransform.scRNA)
    ST <- .spotlight_preprocess(spatial, "Spatial", sctransform.spatial)
    
    markers <- .seurat_find_markers(sobj, cell_type_col)
    
    deconv <- spotlight_deconvolution(se_sc=sc, counts_spatial=ST@assays$Spatial@counts,
        clust_vr="cell_type", cluster_markers=cluster_markers_all, cl_n=100,
        hvg=3000, ntop=NULL, transf="uv", method="nsNMF", min_cont=0.09)
    
}

## Make Seurat object and run sctransform
## The spotlight demo calls for DR and clustering, but doesn't seem necessary
.spotlight_preprocess <- function(expr, assay = "RNA", sctransform = FALSE) {
    sobj <- Seurat::CreateSeuratObject(as.matrix(expr), assay=assay)
    
    if (sctransform)
        sobj <- Seurat::SCTransform(sobj, assay=assay, verbose=FALSE)
    
    # sobj <- Seurat::RunPCA(sobj, verbose = FALSE)
    # sobj <- Seurat::RunUMAP(sobj, dims = 1:30, verbose = FALSE)
    # sobj <- Seurat::FindNeighbors(sobj, dims = 1:30, verbose = FALSE)
    # sobj <- Seurat::FindClusters(sobj, verbose = FALSE)
}

.seurat_find_markers <- function(sobj) {
    Seurat::Idents(object = sobj) <- sobj@meta.data$cell_type
    markers <- Seurat::FindAllMarkers(object=sobj, assay="RNA", slot="data",
        verbose = FALSE, only.pos = TRUE, logfc.threshold = 1, min.pct = 0.9)
}

#' @export
#' @rdname runSPOTlight
setGeneric("runSPOTlight", signature=c("spatial", "scRNA"), function(spatial, scRNA, ...) standardGeneric("runSPOTlight"))

#' @export
#' @rdname runSPOTlight
setMethod("runSPOTlight", signature("ANY", "ANY"), .run_spotlight)

#' @export
#' @rdname runSPOTlight
#' @importFrom SummarizedExperiment assay
setMethod("runSPOTlight", signature(spatial="SummarizedExperiment", scRNA="ANY"), function(spatial, scRNA, ..., spatial.assay="counts") {
    .run_spotlight(assay(spatial, spatial.assay), scRNA, ...)
})

#' @export
#' @rdname runSPOTlight
#' @importFrom SummarizedExperiment assay
setMethod("runSPOTlight", signature(spatial="SummarizedExperiment", scRNA="SummarizedExperiment"), function(spatial, scRNA, ..., spatial.assay="counts", scRNA.assay="counts") {
    .run_spotlight(assay(spatial, spatial.assay), assay(scRNA, scRNA.assay), ...)
})

.extract_cell_types <- function(sce, cell_types) {
    if (length(cell_types) == 1) {
        if (cell_types %in% colnames(colData(sce))) {
            colData(sce)[[cell_types]]
        } else {
            stop("The provided cell type column does not exist.")
        }
    } else if (length(cell_types) == ncol(sce)) {
        
    } else {
        stop("The provided list of cell types does not match the number of cells.")
    }
}

## Check we have Seurat and SPOTlight installed
.check_spotlight_environment <- function() {
    msg <- list("Cannot run SPOTlight without missing dependencies. Please install before continuing:")
    
    if (!("Seurat" %in% rownames(installed.packages()))) {
        msg <- append(msg, "  > install.packages(\"Seurat\")")
    }
    
    if (!("SPOTlight" %in% rownames(installed.packages()))) {
        msg <- append(msg, "  > devtools::install_github(\"https://github.com/MarcElosua/SPOTlight\")")
    }
    
    if (length(msg) > 1) {
        stop(paste(msg, collapse='\n'))
    }
}
