
#' Load a Visium spatial dataset as a SingleCellExperiment.
#' 
#' @param dirname Path to spaceranger output directory (e.g. "sampleID/outs/").
#'   This directory must contain the counts matrix and feature/barcode TSVs in
#'   \code{filtered_feature_bc_matrix/}, and the spot positions at
#'   \code{spatial/tissue_positions_list.csv}. (These are default locations for
#'   spaceranger outputs.)
#'   
#' @return SingleCellExperiment containing the counts matrix in \code{counts}
#'   and spatial data in \code{colData}. Array coordinates for each spot are
#'   stored in columns \code{row} and \code{col}, while image coordinates are
#'   stored in columns \code{imagerow} and \code{imagecol}.
#'   
#' @details We store two variables associated with downstream BayesSpace
#'   functions in a list called \code{BayesSpace.data} in the
#'   SingleCellExperiment's \code{metadata}.
#'   \itemize{
#'     \item \code{platform} is set to "Visium", and is used to determine spot
#'       layout and neighborhood structure.
#'     \item \code{is.enhanced} is set to \code{FALSE} to denote the object
#'       contains spot-level data.
#'   }
#'   
#' @examples
#' \dontrun{
#' sce <- readVisium("path/to/outs/")
#' }
#' 
#' @export
#' @importFrom utils read.csv read.table
#' @importFrom scater uniquifyFeatureNames
#' @importFrom Matrix readMM
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors metadata metadata<-
readVisium <- function(dirname) {
    spatial_dir <- file.path(dirname, "spatial")
    matrix_dir <- file.path(dirname, "filtered_feature_bc_matrix")
    
    if (!dir.exists(matrix_dir)) stop(paste0("Matrix directory does not exist: ", matrix_dir))
    if (!dir.exists(spatial_dir)) stop(paste0("Spatial directory does not exist: ", spatial_dir))
    
    colData <- read.csv(file.path(spatial_dir, "tissue_positions_list.csv"), header=FALSE)
    
    ## TODO: using spatialLIBD conventions here (legacy), but should eventually
    ##   update to canonical spaceRanger names:
    ##   c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    ##   https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/images
    colnames(colData) <- c("spot", "in_tissue", "row", "col", "imagerow", "imagecol")
    rownames(colData) <- colData$spot
    colData <- colData[colData$in_tissue > 0, ]
    
    rowData <- read.table(file.path(matrix_dir, "features.tsv.gz"), header=FALSE)
    colnames(rowData) <- c("gene_id", "gene_name", "feature_type")
    rowData <- rowData[, c("gene_id", "gene_name")]
    rownames(rowData) <- scater::uniquifyFeatureNames(rowData$gene_id, rowData$gene_name)
    
    counts <- Matrix::readMM(file.path(matrix_dir, "matrix.mtx.gz"))
    barcodes <- read.table(file.path(matrix_dir, "barcodes.tsv.gz"), header=FALSE)
    colnames(counts) <- barcodes$V1
    rownames(counts) <- rownames(rowData)
    counts <- counts[, rownames(colData)]
    
    sce <- SingleCellExperiment(assays=list(counts=counts),
                                rowData=rowData,
                                colData=colData)
    
    metadata(sce)$BayesSpace.data <- list()
    metadata(sce)$BayesSpace.data$platform <- "Visium"
    metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
    
    sce
}
