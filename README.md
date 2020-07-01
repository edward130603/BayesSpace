# Spatial transcriptome tools

- `cluster()`: cluster neighboring spots 
- `deconvolve()`: impute cell-level expression vectors from spots


See vignettes/maynard.Rmd for an example that should be fully compatible with the current code implementing the method

## TODOs to pass checks

Last updated 2020-06-29.

### R CMD check
#### Errors (0)
#### Warnings (0)
- [X] ~~Dependency check.~~  
    * ~~Remove mvnfast code if we don't need it~~ DONE
    * ~~SingleCellExperiment declared but not used - should be fixed after
      refactoring `cluster()` to operate on SCE~~ all imports cleaned up
- [X] ~~Undocumented arguments for deconvolve~~
    * ~~Get definitions for new args: model, platform, verbose, jitter_scale, c~~
#### Notes (0)
- [X] ~~Installed package size~~  
    * ~~Currently at 10.8 Mb due to 10.1 Mb testdata~~ Removed full RDS and
      replaced with CSVs of PCs and relevant colData
- [X] ~~DESCRIPTION meta-information.~~  
    * ~~Write an actual Description and pick a license.~~ MIT license
- [X] ~~Unstated vignette dependencies.~~  
    * Removed Maynard vignette to speed up build times during testing
    * Need to pick which example R scripts to include as vignettes, then fix

### R CMD BiocCheck
#### Errors (3)
- [ ] No 'vignettes' directory  
    * TODO: draft vignettes
- [X] ~~Package source tarball exceeds size requirement~~  
    * ~~Currently at 10.6 MB (testdata is 10.1 MB); needs to be <5 MB~~ see above
- [ ] Need runnable examples for exported objects  
    * Wait to do this until after refactoring is complete
- [ ] Maintainer must register at bioconductor support site  
    * Wait to do this until near publication
#### Warnings (0)
- [X] ~~Update R version dependency from 3.5.0 to 4.0~~  
    * Need minimum of 4.0.0 for latest Bioconductor release (3.11)
- [X] ~~Test data over 5 MB in size~~
- [X] ~~Use TRUE/FALSE instead of T/F~~
- [X] ~~Remove set.seed usage from inside functions~~
- [X] ~~Add non-empty value sections to deconvolve docs~~
- [X] ~~Add non-empty value sections to readChain docs~~
#### Notes (5)
- [ ] Avoid sapply(); use vapply()
- [X] ~~Avoid 1:...; use seq_len() or seq_along()~~
- [ ] Keep functions below 50 lines
    * ~~TODO: fix by refactoring out shared input parsing from spatialCluster/spatialEnhance~~
    * TODO: shorten spatialEnhance
- [X] ~~Add NEWS file~~
- [ ] Keep lines <= 80 characters
    * Remaining lines are auto-generated in RcppExports
- [ ] Use multiples of 4 spaces for indents
    * ~~TODO: Run FormatR~~
    * Remaining lines are in `man/`
- [ ] Subscribe maintainer to bioc-devel mailing list

## Notes

Under MacOS Catalina, may be necessary to install gfortran
[directly](https://github.com/fxcoudert/gfortran-for-macOS/releases) rather
than through homebrew in order for libraries to be linkable. I used 8.2 for
Mojave and it worked.
