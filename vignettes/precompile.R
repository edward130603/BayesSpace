#' Pre-compile vignette to avoid recomputing cluster/enhancement
#' 
#' Adapted from examples given in
#' https://ropensci.org/technotes/2019/12/08/precompute-vignettes/
#' https://github.com/ropensci/eia/blob/master/vignettes/precompile.R

knitr::knit("vignettes/BayesSpace.Rmd.orig", output="vignettes/BayesSpace.Rmd")
system("mv BayesSpace-*.png vignettes/")
