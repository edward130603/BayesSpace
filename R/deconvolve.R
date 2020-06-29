#' Deconvolution of spots into cells
#' 
#' Backend calls iterate_deconv(), written in Rcpp.
#' Inputs are the same as `cluster()` except you have to specify xdist and ydist 
#' instead of total dist...(maybe would be better to change `cluster()` to match this)
#' 
#' @param Y A matrix or dataframe with 1 row per spot and 1 column per outcome (e.g. principal component)
#' @param positions A matrix or dataframe with two columns (x, y) that gives the spatial coordinates of the spot
#' @param nrep The maximum number of mcmc iterations
#' @param every TODO: define
#' @param gamma Smoothing parameter. Values in range of 1-3 seem to work well generally
#' @param xdist The distance along x-axis between neighboring spots
#' @param ydist The distance along y-axis between neighboring spots
#' @param q The number of clusters
#' @param init Initial cluster assignments (z's). Must be a vector of length equal to the number of rows of Y and positions
#' @param seed Random state seed
#' @param mu0 Prior mean hyperparameter for mu
#' @param lambda0 Prior precision hyperparam for mu
#' @param alpha Hyperparameter for Wishart distributed precision lambda
#' @param beta Hyperparameter for Wishart distributed precision lambda
#' 
#' @return out TODO: specify
#' 
deconvolve = function(Y, positions, nrep = 1000, every = 1, gamma = 2, 
                      xdist, ydist, q, init, 
                      model = "normal", seed = 100, platform = "visium", 
                      verbose = TRUE, jitter_scale = 5, c = 0.01, 
                      mu0 = colMeans(Y), lambda0 = diag(0.01, nrow = ncol(Y)), 
                      alpha = 1, beta = 0.01) {
  
  d = ncol(Y)
  n0 = nrow(Y)
  positions = as.matrix(positions)
  Y = as.matrix(Y)
  colnames(positions) = c("x", "y")
  
  subspots = ifelse(platform == "visium", 7, 
                    ifelse(platform == "ST", 9, stop("platform should be either 'visium' or 'ST'")))
  
  init1 = rep(init, subspots)
  Y2 = Y[rep(seq_len(n0), subspots), ] #rbind 7 or 9 times
  positions2 = positions[rep(seq_len(n0), subspots), ] #rbind 7 times
  if (platform == "visium"){
    shift = rbind(expand.grid(c(1/3, -1/3), c(1/3,-1/3)), expand.grid(c(2/3, -2/3,0), 0))
  } else {
    shift = rbind(expand.grid(c(1/3, -1/3, 0), c(1/3,-1/3,0)))
  }
  
  shift = t(t(shift)*c(xdist, ydist))
  dist = max(rowSums(abs(shift)))*1.05
  if (platform == "ST"){
    dist = dist/2
  }
  shift_long = shift[rep(seq_len(subspots), each = n0), ]
  positions2[,"x"] = positions2[,"x"] + shift_long[,"Var1"]
  positions2[,"y"] = positions2[,"y"] + shift_long[,"Var2"]
  n = nrow(Y2)
  if (verbose){message("Calculating neighbors...")}
  df_j = sapply(1:n, function(x){which((abs(positions2[,1] -positions2[x,1]) + abs(positions2[,2] - positions2[x,2])) <= dist &  #L1 distance
                                         (abs(positions2[,1] -positions2[x,1]) + abs(positions2[,2] - positions2[x,2])) > 0)-1})
  if (verbose){message("Fitting model...")}
  tdist = ifelse(model == "t", TRUE, FALSE)
  out = iterate_deconv(Y = Y2, df_j = df_j, tdist = tdist, nrep = nrep, n = n, n0 = n0, d = d, gamma = gamma, q = q, init = init1, subspots = subspots, verbose = verbose, jitter_scale = jitter_scale, c = c, mu0 = mu0, lambda0 = lambda0, alpha = alpha, beta = beta)
  out$positions = positions2
  out
}

spatialEnhance <- function(sce, q,
                           assay.type="logcounts", d=15, use.dimred=NULL,
                           pca.method=c("PCA", "denoised", "geneset"),
                           init=NULL, init.method=c("spatialCluster"),
                           positions=NULL, position.cols=c("imagecol", "imagerow")) {
  
  if (is.null(use.dimred)) {
    use.dimred <- "PCA"
    pca.method <- match.arg(pca.method)
    sce <- addPCA(sce, assay.type, pca.method, d)
  }
  
  PCs <- as.matrix(reducedDim(sce, use.dimred))
  d <- min(d, ncol(PCs))
  PCs <- PCs[, 1:d]
  
  if (is.null(positions)) {
    positions <- as.matrix(colData(sce)[position.cols])
    colnames(positions) <- c("x", "y")
  }
  
  # TODO: add checks that all four cols are in SCE
  # TODO: sort out discrepancy between cluster and enhance here
  xdist <- coef(lm(sce$imagecol~sce$col))[2]  # x distance between neighbors
  ydist <- coef(lm(sce$imagerow~sce$row))[2]  # y distance between neighbors
  radius <- (xdist + ydist) * 1.02
  
  # Initialize cluster assignments (use k-means for now)
  if (is.null(init)) {
    init.method <- match.arg(init.method)
    if (init.method == "kmeans") {
      init <- kmeans(PCs, centers = q)$cluster
    } else if (init.method == "spatialCluster") {
      sce <- spatialCluster(sce, q, assay.type, d, use.dimred, pca.method,
                            init, "kmeans", positions, position.cols, radius)
      init <- sce$spatial.cluster
    }
  }
  
  # TODO: add other deconvolve arguments to this function
  deconv <- deconvolve(Y = PCs, positions = positions, nrep = 1000, gamma = 2, 
                       xdist = xdist, ydist = ydist, init = init, q = q)
  
  # Create enhanced SCE
  
  # TODO: auto-run predictExpression and include as assay
  # deconv_PCs <- deconv$Y[[length(deconv$Y)]]
  # expr <- predictExpression(sce, deconv_PCs)
  
  cdata <- as.data.frame(deconv$positions)
  colnames(cdata) <- c("imagecol", "imagerow")
  enhanced <- SingleCellExperiment(assays=list(), 
                                   rowData=rowData(sce), colData=cdata)
  
  reducedDim(enhanced, "PCA") <- deconv$Y[[length(deconv$Y)]]
  enhanced$spatial.cluster <- apply(deconv$z[900:1000, ], 2, Mode)
  
  enhanced
}
