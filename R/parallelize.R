#' Parallelization
#'
#' A convenient wrapper function of \code{BiocParallel} providing easy parallelization.
#'
#' @param BPPARAM 
#' @param cores
#' @param FUN
#' @param type
#' @param verbose
#' @param X
#'
#' @name parallelize
NULL

#' 
#' @return 
#' 
#' @importFrom assertthat assert_that
#'
#' @keywords internal
list2vec <- function(X, sep = ":", collapse = ",", use_names = TRUE) {
  assert_that(!missing(X) && !is.null(X) && is.list(X))
  
  if (is.null(names(X)) && use_names)
    names(X) <- X
  
  if (use_names)
    ret <- sapply(names(X), function(x) paste(x, X[x], sep = sep), USE.NAMES = FALSE)
  else
    ret <- sapply(X, function(x) X[x], USE.NAMES = FALSE)
  
  if (!is.null(collapse) && length(collapse) > 0)
    ret <- paste(ret, collapse = collapse)
  
  ret
}

#' 
#' @return a list containing results of calling "bplapply"
#'
#' @importFrom assertthat assert_that
#' @importFrom BiocParallel SerialParam MulticoreParam SnowParam bplapply
#' @importFrom purrr compact
#'
#' @export
#' @rdname paraLapply
paraLapply <- function(X, FUN, BPPARAM = NULL, cores = 1L, type = c("serial", "fork", "sock", "mpi"), verbose = FALSE, ...) {
  # Sanity check.
  assert_that(!missing(X) && !is.null(X) && length(X) > 0)
  assert_that(!missing(FUN) && !is.null(FUN))
  assert_that(is.null(BPPARAM) || inherits(BPPARAM, "BiocParallelParam"))
  assert_that(length(cores) == 1L && is.integer(cores) && cores > 0L)
  
  # Common arguments for "BiocParallelParam".
  .common.args <- c("stop.on.error", "progressbar", "RNGseed", "timeout", "threshold", "log", "logdir", "resultdir", "jobname")
  
  .args <- list(...)
  
  # Initialize back end for parallization is "BPPARAM" is not provided.
  if (is.null(BPPARAM)) {
    if (verbose)
      message("[DEBUG] Back end for parallelization is not provided. Initializing...")
    
    type <- match.arg(type)
    
    # Get common arguments for "BiocParallelParam".
    .init.args <- compact(.args[.common.args])
    
    # Set the number of cores to be used.
    if (!is.null(cores))
      .init.args[["workers"]] = cores
    
    if (verbose)
      message(sprintf("[DEBUG] Provided effective arguments for creating a parallelization back end: %s", list2vec(.init.args)))
    
    # Check if the operating system is Windows.
    .is.windows <- switch(Sys.info()[["sysname"]], Windows = TRUE, Linux = FALSE, Darwin = FALSE)
    
    # Create an object for parallelization reasonably.
    if (type == "serial") {
      BPPARAM <- do.call(
        SerialParam,
        discard(.init.args, function(x) x %in% c("workers"))
      )
      
      if (verbose)
        message("[DEBUG] A serial back end is created.")
    } else if (type == "mpi") {
      # if (is.null(.init.args[["workers"]]))
      #   .init.args[["workers"]] = mpi.universe.size - 1
      # 
      # BPPARAM <- do.call(
      #   SnowParam,
      #   c(.init.args, list("type" = "MPI"))
      # )
      
      # if (verbose)
      #   message("[DEBUG] An MPI back end is created.")
      
      stop("Error! The MPI mode is not yet available.")
    } else if (type == "sock" || (type == "fork" && .is.windows)) {
      BPPARAM <- do.call(
        SnowParam,
        c(.init.args, list("type" = "SOCK"))
      )
      
      if (verbose)
        message("[DEBUG] A socket back end is created.")
    } else {
      BPPARAM <- do.call(
        MulticoreParam,
        .init.args
      )
      
      if (verbose)
        message("[DEBUG] A fork back end is created.")
    }
  }

  # Extra arguments for FUN.
  .extra.args <- .args[discard(names(.args), function(x) x %in% .common.args)]
  if (verbose)
    message(sprintf("[DEBUG] Provided effective arguments for the customized function: %s", list2vec(.extra.args)))
  
  do.call(
    bplapply,
    c(
      list(
        X = X,
        FUN = FUN,
        BPPARAM = BPPARAM
      ),
      .extra.args
    )
  )
}










