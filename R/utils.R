#' Compute pairwise distances between all spots and return list of neighbors
#' for each spot
#' 
#' @param positions 
#'        (n x 2) matrix of spot coordinates
#' @param radius
#'        The maximum distance for two spots to be considered neighbors 
#' @param method
#'        Distance metric to use
#' 
#' @return List df_j, 
#'         where df_j[[i]] is a vector of zero-indexed neighbors of i 
#'         
#' @importFrom stats dist
find_neighbors <- function(positions, neighborhood.radius,
                           method=c("manhattan", "euclidean")) {
  method <- match.arg(method)
  
  message("Calculating neighbors...")
  pdist <- as.matrix(stats::dist(positions, method=method))
  neighbors <- (pdist <= neighborhood.radius & pdist > 0)
  df_j <- sapply(1:nrow(positions),
                 function(x) as.vector(which(neighbors[x, ])) - 1)
  
  msg <- sprintf("Neighbors were identified for %d out of %d spots.",
                 sum(rowSums(neighbors) > 0), 
                 nrow(positions))
  message(msg)
  
  df_j
}

#' Compute the distance between two neighboring spots
#' 
#' @return double radius
#' 
#' @importFrom stats lm
#' @importFrom stats coef
compute_neighborhood_radius <- function(sce) {
  # TODO: remove hardcoding of columns
  xdist <- coef(lm(sce$imagecol~sce$col))[2]  # x distance between neighbors
  ydist <- coef(lm(sce$imagerow~sce$row))[2]  # y distance between neighbors
  radius <- xdist + ydist + 0.2
  
  radius
}