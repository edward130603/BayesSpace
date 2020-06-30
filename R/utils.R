#' Compute pairwise distances between all spots and return list of neighbors
#' for each spot.
#' 
#' @param positions (n x 2) matrix of spot coordinates.
#' @param radius The maximum distance for two spots to be considered neighbors.
#' @param method Distance metric to use.
#' 
#' @return List df_j, where \code{df_j[[i]]} is a vector of zero-indexed
#'   neighbors of i.
#'         
#' @importFrom stats dist
find_neighbors <- function(positions, radius,
    method = c("manhattan", "euclidean")) {
    
    method <- match.arg(method)
    
    message("Calculating neighbors...")
    pdist <- as.matrix(stats::dist(positions, method=method))
    neighbors <- (pdist <= radius & pdist > 0)
    df_j <- sapply(seq_len(nrow(positions)), 
        function(x) as.vector(which(neighbors[x, ])) - 1)
    
    msg <- "Neighbors were identified for %d out of %d spots." 
    msg <- sprintf(msg, sum(rowSums(neighbors) > 0), nrow(positions))
    message(msg)
    
    df_j
}

#' Estimate the distance between two neighboring spots
#' 
#' Fit linear models between each image pixel coordinate and its corresponding
#' array coordinate to estimate the pixel distance between two spots along
#' each axis. Add these distances to estimate the L1 distance between two
#' spots, then add a small buffer.
#' 
#' @param sce SingleCellExperiment (must include row, col, imagerow, imagecol 
#'   in colData)
#' @param scale.factor Scale estimated L1 difference up by this amount.
#' 
#' @return doubles xdist, ydist, radius
#' 
#' @importFrom stats lm coef
.compute_interspot_distances <- function(sce, scale.factor = 1.02) {
    ## TODO: remove hardcoding of columns
    dists <- list()
    
    dists$xdist <- coef(lm(sce$imagecol ~ sce$col))[2]
    dists$ydist <- coef(lm(sce$imagerow ~ sce$row))[2]
    dists$radius <- (dists$xdist + dists$ydist) * scale.factor
    
    dists
}

#' Find the mode
#' 
#' Used for finding the most frequent cluster for each z
#' 
#' @param x Numeric vector
#'
#' @return mode Numeric scalar, most frequent element in x
#'
Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

#' Run 'blocked' PCA on HVGs and gene sets of interest
#' 
#' PCA is computed for each specified gene set. The PCs from each gene set are 
#' concatenated into a single matrix, and their rotation matrices are combined
#' into a block diagonal rotation matrix.
#' 
#' If a set of highly variable genes (HVGs) is not included in the list of 
#' gene sets (using the name \code{hvgs}), the top \code{n_hvgs} will be
#' computed using scran and added to the list of gene sets.
#' 
#' @param expr Expression/counts (matrix or SingleCellExperiment)
#' @param n_pcs Number of PCs per gene set (TODO: optionally specify list, n 
#' per gene set)
#' @param n_hvgs Number of HVGs to use in HVG set
#' @param genesets List of other gene sets to include
#' 
#' @return Returns list with names:
#'   * \code{x} Concatenated principal components from each gene set
#'   * \code{rotation} Rotation matrix (block diagonal for each set of PCs)
#'   * \code{genesets} List of gene sets in each block
#'   
#' @importFrom Matrix bdiag
runGenesetPCA <- function(expr, n_pcs = 20, n_hvgs = 2000, genesets = list()) {
    ## Compute top HVGs, excluding any that appear in user-specified gene sets
    if (n_hvgs >= 1 && !("hvgs" %in% names(genesets))) {
        dec <- scran::modelGeneVar(expr)
        hvgs <- scran::getTopHVGs(dec, n=n_hvgs)
        others <- purrr::reduce(genesets, union)
        genesets[["hvgs"]] <- setdiff(hvgs, others)
    }
    
    ## Run PCA on each set of genes
    .runSubsetPCA <- function(geneset, name) {
        subset <- expr[geneset, ]
        pca <- BiocSingular::runPCA(t(subset), rank=n_pcs)
        
        colnames(pca$x) <- paste0(name, "_", colnames(pca$x))
        colnames(pca$rotation) <- paste0(name, "_", colnames(pca$rotation))
        
        pca
    }
    pcs <- purrr::imap(genesets, .runSubsetPCA)
    
    ## Concatenate PCs and assemble rotation matrices
    out <- list()
    out$x <- purrr::reduce(purrr::map(pcs, "x"), cbind)
    out$rotation <- purrr::reduce(purrr::map(pcs, "rotation"), bdiag)
    out$genesets <- genesets
    
    ## TODO: add PCs to reducedDims of SCE
    ## (for compatibility with other methods in addPCA)
    
    return(out)
}

#' Add PCA output to a SingleCellExperiment
#' 
#' Supports 'basic' PCA via scater, denoised PCA via scran, and blocked geneset
#' PCA. Results will be stored in \code{reducedDim(sce, 'PCA')}.
#' 
#' @param sce SingleCellExperiment
#' @param assay.type Assay in \code{sce} to compute PCA on
#' @param pca.method PCA method to apply
#' @param d Number of principal components to keep
#' 
#' @return Returns \code{sce} with PCs added to its \code{reducedDims}
addPCA <- function(sce, assay.type, pca.method, d = 15) {
    if (pca.method == "PCA") {
        sce <- scater::runPCA(sce, exprs_values=assay.type, ncomponents=d)
        
    } else if (pca.method == "denoised") {
        dec <- scran::modelGeneVar(sce, assay.type=assay.type)
        top <- scran::getTopHVGs(dec, prop=0.1)  # TODO: parameterize
        sce <- scran::denoisePCA(sce, technical=dec, subset.row=top, 
            assay.type=assay.type)
        
    } else if (pca.method == "geneset") {
        stop("geneset PCA not yet supported")
        
    } else {
        stop("Unsupported PCA method: %s", pca.method)
    }
    
    sce
}

## Prepare cluster/deconvolve inputs from SingleCellExperiment object
#' @importFrom SingleCellExperiment reducedDimNames
#' @importFrom purrr imap
.prepare_inputs <- function(sce, use.dimred = "PCA", d = 15, 
    positions = NULL, position.cols = c("imagecol", "imagerow"), 
    radius = NULL, xdist = NULL, ydist = NULL) {
    
    inputs <- list()
    
    if (!(use.dimred %in% reducedDimNames(sce))) 
        stop(sprintf("reducedDim %s not found in input SCE", use.dimred))
    
    PCs <- reducedDim(sce, use.dimred)
    d <- min(d, ncol(PCs))
    inputs$PCs <- PCs[, seq_len(d)]
    
    if (is.null(positions)) 
        positions <- as.matrix(colData(sce)[position.cols])
    
    colnames(positions) <- c("x", "y")
    inputs$positions <- positions
    
    dists <- .compute_interspot_distances(sce)
    dists <- imap(dists, function(d, n) ifelse(is.null(get(n)), d, get(n)))
    inputs <- c(inputs, dists)
    
    inputs
}
