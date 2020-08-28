#' Tuning the choice of q (number of clusters) before running spatialCluster
#'
#' Before running \code{spatialCluster()}, we recommend tuning the choice of
#' \code{q} by choosing the \code{q} that maximizes the model's negative log
#' likelihood over early iterations. \code{qTune()} computes the average
#' negative log likelihood for a range of q values over iterations 100:1000, and
#' \code{qPlot()} displays the results.
#' 
#' @param sce A SingleCellExperiment object containing the spatial data.
#' @param qs The values of q to evaluate.
#' @param min_rep,max_rep Integers specifying the range of repetitions to
#'   compute 
#' @param force.retune If specified, existing tuning values in \code{sce} will
#'   be overwritten.
#' @param ... Other parameters are passed to \code{spatialCluster()}.
#' 
#' @return \code{qTune()} returns a modified \code{sce} with tuning log 
#'   likelihoods stored as an attribute named \code{"q.logliks"}.
#'   
#'   \code{qPlot()} returns a ggplot object.
#'         
#' @details 
#' \code{qTune()} takes the same parameters as \code{spatialCluster()} and will
#'   run the MCMC clustering algorithm up to \code{max_rep} iterations for each
#'   value of \code{q}. The first \code{min_rep} iterations are discarded as
#'   burn-in and the log likelihood is averaged over the remaining iterations.
#'   
#' \code{qPlot()} plots the computed negative log likelihoods as a function of
#'   q. If \code{qTune()} was run previously, i.e. there exists an attribute of
#'   \code{sce} named \code{"q.logliks"}, the pre-computed results are
#'   displayed. Otherwise, or if \code{force.retune} is specified,
#'   \code{qplot()} will automatically run \code{qTune()} before plotting (and
#'   can take the same parameters as \code{spatialCluster()}.
#' 
#' @examples
#' set.seed(149)
#' sce <- exampleSCE()
#' sce <- qTune(sce, seq(3, 7))
#' qPlot(sce)
#'
#' @name qTune
NULL

#' @importFrom ggplot2 ggplot geom_line geom_point xlab ylab labs theme_bw
#' 
#' @export
#' @rdname qTune
qPlot <- function(sce, qs=seq(3, 7), force.retune=FALSE, ...) {
    if (!("q.logliks" %in% names(attributes(sce))) || force.retune) {
        sce <- qTune(sce, qs, ...)        
    }
    
    logliks <- attr(sce, "q.logliks")
    qplot <- ggplot(data=logliks, aes(x=q, y=-loglik)) +
        geom_line() +
        geom_point() +
        xlab("Number of clusters (q)") +
        ylab("Negative log likelihood") +
        labs(title="spatialCluster likelihood as a function of q") +
        theme_bw()
    
    qplot
}


#' @importFrom purrr compact discard
#' 
#' @export
#' @rdname qTune
qTune <- function(sce, qs=seq(3, 7), min_rep=100, max_rep=1000, ...) {
    ## TODO: refactor args into a ClusterConfig object and store as sce attribute
    args <- list(...)
    print(args)
    
    input.args <- c("use.dimred", "d", "positions", "position.cols", "radius")
    input.args <- compact(args[input.args])
    inputs <- do.call(.prepare_inputs, c(list(sce=sce), input.args))
    
    init.args <- c("init", "init.method")
    init.args <- compact(args[init.args])
    
    cluster.args <- discard(names(args), function(x) {x %in% c(names(input.args), names(init.args))})
    cluster.args <- compact(args[cluster.args])
    cluster.args$nrep <- max_rep
    
    logliks <- list()
    for (q in qs) {
        sce <- do.call(.init_cluster, c(list(sce=sce, q=q, inputs=inputs), init.args))
        
        input.args <- list(Y=inputs$PCs, positions=inputs$positions, 
                           radius=inputs$radius, q=q, init=sce$cluster.init)
        
        results <- do.call(cluster, c(input.args, cluster.args))
        logliks[[q]] <- data.frame(q=q, loglik=mean(results$plogLik[min_rep:max_rep]))
    }
    
    logliks <- do.call(rbind, logliks)
    attr(sce, "q.logliks") <- logliks
    sce
}
