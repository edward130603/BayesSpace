#' Pre-compile vignette to avoid recomputing cluster/enhancement
#' 
#' Adapted from examples given in
#' https://ropensci.org/technotes/2019/12/08/precompute-vignettes/
#' https://github.com/ropensci/eia/blob/master/vignettes/precompile.R

# knitr::knit("vignettes/BayesSpace.Rmd.orig", output="vignettes/BayesSpace.Rmd")
# system2("mv", c("figures/BayesSpace-*.png", "vignettes/figures/"))

# knitr::knit("vignettes/thrane_melanoma.Rmd.orig", output="vignettes/thrane_melanoma.Rmd")
# system2("mv", c("figures/melanoma-*", "vignettes/figures/"))

# knitr::knit("vignettes/maynard_DLPFC.Rmd.orig", output="vignettes/maynard_DLPFC.Rmd")
# system2("mv", c("figures/maynard_DLPFC-*.png", "vignettes/figures/"))

knitr::knit("vignettes/ji_SCC.Rmd.orig", output="vignettes/ji_SCC.Rmd")
system2("mv", c("figures/ji_SCC-*.png", "vignettes/figures/"))