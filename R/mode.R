#' Find the mode
#' 
#' Used for finding the most frequent cluster for each z
#' 
#' @param x Numeric vector
#' @return mode Numeric scalar, most frequent element in x
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}