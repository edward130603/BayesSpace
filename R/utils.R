#' Compute pairwise distances between all spots and return list of neighbors
#' for each spot.
#' 
#' @param positions (n x 2) matrix of spot coordinates.
#' @param radius The maximum distance for two spots to be considered neighbors.
#' @param method Distance metric to use.
#' 
#' @return List df_j, where df_j[[i]] is a vector of zero-indexed neighbors of i.
#'         
#' @importFrom stats dist
find_neighbors <- function(positions, radius,
                           method=c("manhattan", "euclidean")) {
  method <- match.arg(method)
  
  message("Calculating neighbors...")
  pdist <- as.matrix(stats::dist(positions, method=method))
  neighbors <- (pdist <= radius & pdist > 0)
  df_j <- sapply(1:nrow(positions),
                 function(x) as.vector(which(neighbors[x, ])) - 1)
  
  msg <- sprintf("Neighbors were identified for %d out of %d spots.",
                 sum(rowSums(neighbors) > 0), 
                 nrow(positions))
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
#' @return double radius
#' 
#' @importFrom stats lm coef
compute_neighborhood_radius <- function(sce, scale.factor=1.02) {
  # TODO: remove hardcoding of columns
  xdist <- coef(lm(sce$imagecol~sce$col))[2]  # x distance between neighbors
  ydist <- coef(lm(sce$imagerow~sce$row))[2]  # y distance between neighbors
  radius <- xdist + ydist
  radius <- radius * scale.factor
  
  radius
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

#' Run "blocked" PCA on HVGs and gene sets of interest
#' 
#' PCA is computed for each specified gene set. The PCs from each gene set are 
#' concatenated into a single matrix, and their rotation matrices are combined
#' into a block diagonal rotation matrix.
#' 
#' If a set of highly variable genes (HVGs) is not included in the list of 
#' gene sets (using the name `hvgs`), the top `n_hvgs` will be computed using 
#' scran and added to the list of gene sets.
#' 
#' @param expr Expression/counts (matrix or SingleCellExperiment)
#' @param n_pcs Number of PCs per gene set (TODO: optionally specify list, n 
#' per gene set)
#' @param n_hvgs Number of HVGs to use in HVG set
#' @param genesets List of other gene sets to include
#' 
#' @return Returns list with names:
#'   * `x` Concatenated principal components from each gene set
#'   * `rotation` Rotation matrix (block diagonal for each set of PCs)
#'   * `genesets` List of gene sets in each block
#'   
#' @importFrom Matrix bdiag
runGenesetPCA <- function(expr, n_pcs=20, n_hvgs=2000, genesets=list()) {
  # Compute top n HVGs and exclude any that appear in user-specified gene sets
  if (n_hvgs >= 1 && !("hvgs" %in% names(genesets))) {
    dec <- scran::modelGeneVar(expr)
    hvgs <- scran::getTopHVGs(dec, n=n_hvgs)
    others <- purrr::reduce(genesets, union)
    genesets[["hvgs"]] <- setdiff(hvgs, others)
  }
  
  # Run PCA on each set of genes
  .runSubsetPCA <- function(geneset, name) {
    subset <- expr[geneset, ]
    pca <- BiocSingular::runPCA(t(subset), rank=n_pcs)
    
    colnames(pca$x) <- paste0(name, "_", colnames(pca$x))
    colnames(pca$rotation) <- paste0(name, "_", colnames(pca$rotation))
    
    pca
  }
  pcs <- purrr::imap(genesets, .runSubsetPCA)
  
  # Concatenate PCs and assemble rotation matrices
  out <- list()
  out$x <- purrr::reduce(purrr::map(pcs, "x"), cbind)
  out$rotation <- purrr::reduce(purrr::map(pcs, "rotation"), bdiag)
  out$genesets <- genesets
  
  # TODO: add PCs to reducedDims of SCE 
  #       (for compatibility with other methods in addPCA)
  
  return(out)
}

#' Add PCA output to a SingleCellExperiment
#' 
#' Supports "basic" PCA via scater, denoised PCA via scran, and blocked geneset
#' PCA. Results will be stored in `reducedDim(sce, "PCA")`
#' 
#' @param sce SingleCellExperiment
#' @param assay.type Assay in `sce` to compute PCA on
#' @param pca.method PCA method to apply
#' @param d Number of principal components to keep
#' 
#' @return Returns `sce` with PCs added to its `reducedDims`
addPCA <- function(sce, assay.type, pca.method, d=15) {
  if (pca.method == "PCA") {
    sce <- scater::runPCA(sce, exprs_values=assay.type, ncomponents=d)
    
  } else if (pca.method == "denoised") {
    dec <- scran::modelGeneVar(sce, assay.type=assay.type)
    top <- scran::getTopHVGs(dec, prop=0.1)  # TODO: parameterize
    sce <- scran::denoisePCA(sce, technical=dec, subset.row=top, assay.type=assay.type)
    
  } else if (pca.method =="geneset") {
    stop("geneset PCA not yet supported")
    
  } else {
    stop("Unsupported PCA method: %s", pca.method)
  }
  
  sce
}
