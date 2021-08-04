#' Pre-compile vignettes to avoid recomputing cluster/enhancement on install
#' 
#' Adapted from examples given in
#' https://ropensci.org/technotes/2019/12/08/precompute-vignettes/
#' https://github.com/ropensci/eia/blob/master/vignettes/precompile.R
#'
#' Vignettes must be pre-compiled for the package build time to meet
#' Bioconductor requirements. Use the following to add or update a vignette:
#' 
#' 1) Save the vignette (in Rmarkdown) to vignettes/${vignette}.Rmd.orig.
#'    Make sure to include `fig.path=figures/${vignette}-"` in the knitr
#'    options to make migrating the figures easier.
#'
#' 2) Following the below examples, knit the vignette to Rmarkdown. This will
#'    convert the code blocks to plaintext and save the figures to disk. Make
#'    sure to copy the figures over afterwards.
#'
#' 3) Add the new vignette to the list of "Examples" in _pkgdown.yml, so it
#'    will show up in the website's menu bar.
#'
#' 4) Rebuild the website by running `pkgdown::build_site()`. This will render
#'    the knitted vignettes to HTML.
#'
#' 5) Optionally, preview the website with `pkgdown::preview_site()`.
#'
#' 6) Finally, push the updated site with `pkgdown::deploy_to_branch()`.


# knitr::knit("vignettes/BayesSpace.Rmd.orig", output="vignettes/BayesSpace.Rmd")
# system2("mv", c("figures/BayesSpace-*.png", "vignettes/figures/"))

# knitr::knit("vignettes/thrane_melanoma.Rmd.orig", output="vignettes/thrane_melanoma.Rmd")
# system2("mv", c("figures/melanoma-*", "vignettes/figures/"))

# knitr::knit("vignettes/maynard_DLPFC.Rmd.orig", output="vignettes/maynard_DLPFC.Rmd")
# system2("mv", c("figures/maynard_DLPFC-*.png", "vignettes/figures/"))

# knitr::knit("vignettes/ji_SCC.Rmd.orig", output="vignettes/ji_SCC.Rmd")
# system2("mv", c("figures/ji_SCC-*.png", "vignettes/figures/"))

knitr::knit("vignettes/joint_clustering.Rmd.orig", output="vignettes/joint_clustering.Rmd")
system2("mv", c("figures/joint_clustering-*.png", "vignettes/figures/"))

system2("rmdir", "figures/")
