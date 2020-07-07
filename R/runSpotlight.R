#' Run SPOTlight cell type deconvolution on a spatial SCE
#' 
#' desc
#' 
#' @param 
#' 
#' @export
#' @importFrom 

runSPOTlight <- function() {
    .check_spotlight_environment()
}

.check_spotlight_environment <- function() {
    msg <- list("Cannot run SPOTlight without missing dependencies. Please install before continuing:")
    
    if (!("Seurat" %in% rownames(installed.packages()))) {
        msg <- append(msg, "  > install.packages(\"Seurat\")")
    }
    
    if (!("SPOTlight" %in% rownames(installed.packages()))) {
        msg <- append(msg, "  > devtools::install_github(\"https://github.com/MarcElosua/SPOTlight\")")
    }
    
    if (length(msg) > 1) {
        stop(paste(msg, collapse='\n'))
    }
}
