# BayesSpace (development version)

## New features

- `getRDS()` gains a `cache` parameter. When `TRUE`, the RDS is cached locally
  using `BiocFileCache`.

## Minor improvements and fixes

- Addressed reviewer concerns (https://github.com/Bioconductor/Contributions/issues/1624)
    * Updated stop/warning/message statements to remove redundancies and
      unnecessary use of `paste()`.
    * Removed inline conditional statements.
    * Cache downloaded RDS in `getRDS()` (see above).

# BayesSpace 0.99.5

## Minor improvements and fixes

- In `spatialCluster()` and `spatialEnhance()`, setting `burn.in` equal to
  `nrep` now raises an error.

# BayesSpace 0.99.4

## New features

- `enhanceFeatures()` now takes an `nrounds` parameter that corresponds to the
  same parameter in xgboost. If `nrounds` is set to 0, we automatically tune
  the parameter using a train/test split for improved feature prediction.
- `spatialCluster()` and `spatialEnhance()` both gain a `burn.in` parameter
  specifying the number of MCMC iterations to exclude when aggregating cluster
  labels and enhanced PCs.
- In `clusterPlot()`, `label` now accepts factors and vectors of strings, in
  addition to numeric vectors or a column name in `colData`.
- Additional vignettes provided for reproducing the analyses of the melanoma,
  dorsolateral prefrontal cortex, and squamous cell carcinoma datasets presented
  in the bioRxiv manuscript.

## Minor improvements and fixes

- The internal layout of subspots is now correctly oriented (accounting for
  vertical flip of spot coordinates) when using spatial plot functions on
  enhanced Visium data.
- In `spatialEnhance()`, PCs are now averaged over the MCMC iterations
  (excluding the burn-in period).
- In `enhanceFeatures()`, negative expression is now clipped to 0.
- `spatialPreprocess()` now adds a boolean `is.HVG` column to `rowData`.
- In `featurePlot()`, additional arguments to `geom_polygon()` are correctly
  passed through.

# BayesSpace 0.99.3

## Minor improvements and fixes

- Updated `README.md` to include system requirements, additional installation
  details, and link to vignette with demonstration of package functions, per
  journal guidelines.

# BayesSpace 0.99.2

## Minor improvements and fixes

- `spatialEnhance()` incorrectly added row offset to spot column coordinate
  when generating subspot colData, and vice versa. This resulted in subspots
  being reflected over y=x in spatial plots, and has been fixed.
- Figures in the demonstration vignette have been updated with this fix.

# BayesSpace 0.99.1

## Minor improvements and fixes

- Removed Maintainer field from DESCRIPTION to adhere to Bioconductor
  guidelines.

# BayesSpace 0.99.0

## New features

- Initial Bioconductor submission
