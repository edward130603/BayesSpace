# BayesSpace

BayesSpace provides tools for clustering and enhancing the resolution of spatial gene expression experiments. 

BayesSpace clusters a low-dimensional representation of the gene expression
matrix, incorporating a spatial prior to encourage neighboring spots to cluster
together. The method can enhance the resolution of the low-dimensional
representation into "sub-spots", for which features such as gene expression or
cell type composition can be imputed.

## Installation

BayesSpace will be available on BioConductor soon. Until then, it can be installed with `devtools`:

```
# Install devtools if necessary
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("edward130603/BayesSpace")
```

### Notes

Under MacOS Catalina, it may be necessary to install gfortran
[directly](https://github.com/fxcoudert/gfortran-for-macOS/releases) rather
than through homebrew in order for libraries to be linkable. I used 8.2 for
Mojave and it worked.
