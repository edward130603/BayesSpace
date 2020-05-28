# Spatial transcriptome tools

- `cluster()`: cluster neighboring spots 
- `deconvolve()`: impute cell-level expression vectors from spots


See vignettes/maynard.Rmd for an example that should be fully compatible with the current code implementing the method

## TODOs to pass checks

### R CMD check
#### Errors (0)
#### Warnings (1)
- [ ] Dependency check.  
    * Remove mvnfast code if we don't need it
    * SingleCellExperiment declared but not used - should fix after refactoring
      `cluster()` to operate on SCE
#### Notes (3)
- [ ] Installed package size  
    * Currently at 10.8 Mb due to 10.1 Mb testdata
- [ ] DESCRIPTION meta-information.  
    * Write an actual Description and pick a license.
- [ ] Unstated vignette dependencies.  
    * Current vignette is just a working example, will clean up later

### R CMD BiocCheck
#### Errors (3)
- [ ] Package source tarball exceeds size requirement  
    * Currently at 10.6 MB (testdata is 10.1 MB); needs to be <5 MB
- [ ] Need runnable examples for exported objects  
    * Wait to do this until after refactoring is complete
- [ ] Maintainer must register at bioconductor support site  
    * Wait to do this until near publication
#### Warnings (5)
- [ ] Update R version dependency from 3.5.0 to 4.0  
    * TODO: determine minimum R version and run checks under it
- [ ] Test data over 5 MB in size
- [ ] Use TRUE/FALSE instead of T/F
- [ ] Remove set.seed usage from inside functions
- [ ] Add non-empty value sections to deconvolve docs
#### Notes (7)
- [ ] Avoid sapply(); use vapply()
- [ ] Avoid 1:...; use seq_len() or seq_along()
- [ ] Keep functions below 50 lines
    * Offenders are old implementations w/o Rcpp, will probably remove
- [ ] Add NEWS file
- [ ] Keep lines <= 80 characters
- [ ] Use multiples of 4 spaces for indents
- [ ] Subscribe maintainer to bioc-devel mailing list
