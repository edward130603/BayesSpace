# CHANGELOG - BayesSpace

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- The `nrounds` parameter in xgboost can now be tuned automatically within
  `enhanceFeatures()` for improved feature prediction.
- Additional vignettes provided for reproducing the analyses of the melanoma,
  dorsolateral prefrontal cortex, and squamous cell carcinoma datasets presented
  in the bioRxiv manuscript.

### Fixed

- When using the spatial plot functions on enhanced Visium data, the internal
  layout of the subspots was incorrectly flipped vertically. This has been
  fixed.

## [0.99.3] - 2020-09-21

### Changed

- Updated README to include system requirements, additional installation
  details, and link to vignette with demonstration of package functions, per
  journal guidelines.

## [0.99.2] - 2020-09-09

### Fixed

- `spatialEnhance()` incorrectly added row offset to spot column coordinate
  when generating subspot colData, and vice versa. This resulted in subspots
  being reflected over y=x in spatial plots, and has been fixed.
- Figures in the demonstration vignette have been updated with this fix.

## [0.99.1] - 2020-09-08

### Changed

- Removed Maintainer field from DESCRIPTION to adhere to Bioconductor
  guidelines.

## [0.99.0] - 2020-09-06

### Added

- Initial Bioconductor submission
