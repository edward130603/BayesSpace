#' Parallelization
#'
#' A convenient wrapper function of \code{BiocParallel} providing easy
#'   parallelization.
#'
#' @param X Any object for which methods \code{length}, \code{[}, and \code{[[}
#'   are implemented (passed to \code{bplapply}).
#' @param FUN The \code{function} to be applied to each element of X
#'   (passed to \code{bplapply}).
#' @param BPPARAM An optional \code{BiocParallelParam} instance determining
#'   the parallel back-end to be used during evaluation, or a list of
#'   \code{BiocParallelParam} instances, to be applied in sequence for
#'   nested calls to \code{BiocParallel} functions.
#' @param cores The number of threads to use. The results are invariate to the
#'   value of \code{cores}.
#' @param type One of "serial", "fork", or "sock". When \code{cores} is one,
#'   \code{type} is always "serial". Both "fork" and "sock" are for
#'   multi-threading. "fork" is faster, but only supports linux and macos.
#'   "sock" supports linux, macos, and windows.
#' @param verbose Whether to print debug information or not.
#' @param ... Additional parameters passed to \code{bplapply}.
#'
#' @name parallelize
NULL

#' @return See \code{lapply}.
#'
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom BiocParallel SerialParam MulticoreParam SnowParam bplapply
#' @importFrom purrr compact
#'
#' @rdname parallelize
paraLapply <- function(X, FUN, BPPARAM = NULL, cores = 1L, type = c("serial", "fork", "sock"), verbose = FALSE, ...) {
  # Sanity check.
  assert_that(!missing(X) && !is.null(X) && length(X) > 0)
  assert_that(!missing(FUN) && !is.null(FUN))
  assert_that(is.null(BPPARAM) || inherits(BPPARAM, "BiocParallelParam"))
  assert_that(length(cores) == 1L && is.numeric(cores) && cores > 0L)

  # Common arguments for "BiocParallelParam".
  .common.args <- c("stop.on.error", "progressbar", "RNGseed", "timeout", "threshold", "log", "logdir", "resultdir", "jobname")

  .args <- list(...)

  # Initialize back end for parallization is "BPPARAM" is not provided.
  if (is.null(BPPARAM)) {
    if (verbose) {
      message("[DEBUG] Back end for parallelization is not provided. Initializing...")
    }

    type <- match.arg(type)

    # Get common arguments for "BiocParallelParam".
    .init.args <- compact(.args[.common.args])

    # Set the number of cores to be used.
    .init.args[["workers"]] <- as.integer(cores)

    if (verbose) {
      message(sprintf("[DEBUG] Provided effective arguments for creating a parallelization back end: %s", .list2vec(.init.args)))
    }

    # Check if the operating system is Windows.
    .is.windows <- switch(Sys.info()[["sysname"]],
      Windows = TRUE,
      Linux = FALSE,
      Darwin = FALSE
    )

    # Create an object for parallelization reasonably.
    if (type == "serial" || cores == 1L) {
      BPPARAM <- do.call(
        SerialParam,
        .init.args[discard(names(.init.args), function(x) x %in% c("workers"))]
      )

      if (verbose) {
        message("[DEBUG] A serial back end is created.")
      }
    } else if (type == "mpi") {
      stop("Error! The MPI mode is not yet available.")
    } else if (type == "sock" || (type == "fork" && .is.windows)) {
      BPPARAM <- do.call(
        SnowParam,
        c(.init.args, list("type" = "SOCK"))
      )

      if (verbose) {
        message("[DEBUG] A socket back end is created.")
      }
    } else {
      BPPARAM <- do.call(
        MulticoreParam,
        .init.args
      )

      if (verbose) {
        message("[DEBUG] A fork back end is created.")
      }
    }
  }

  # Extra arguments for FUN.
  .extra.args <- .args[discard(names(.args), function(x) x %in% .common.args)]
  if (verbose) {
    message(sprintf("[DEBUG] Provided effective arguments for the customized function: %s", .list2vec(.extra.args)))
  }

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
