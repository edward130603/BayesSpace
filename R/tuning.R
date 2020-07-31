
#' @importFrom ggplot2 ggplot geom_line geom_point xlab ylab labs
#' @export
qPlot <- function(sce, qs=seq(3, 7), force.retune=FALSE, ...) {
    if (!("q.logliks" %in% names(attributes(sce))) || force.retune) {
        sce <- qTune(sce, qs, ...)        
    }

    logliks <- attribute(sce, "q.logliks")
    qplot <- ggplot(data=logliks, aes(x=q, y=-loglik)) +
        geom_line() +
        geom_point() +
        xlab("Number of clusters (q)") +
        ylab("Negative log likelihood") +
        labs(title="spatialCluster likelihood as a function of q")
    
    qplot
}


#' @importFrom purrr compact discard
#' @export
qTune <- function(sce, qs=seq(3, 7), min_rep=100, max_rep=1000, ...) {
    ## TODO: refactor args into a ClusterConfig object and store as sce attribute
    args <- list(...)
    
    input.args <- c("use.dimred", "d", "positions", "position.cols", "radius")
    input.args <- compact(args[input.args])
    inputs <- do.call(.prepare_inputs, c(list(sce=sce), input.args))
    
    init.args <- c("init", "init.method")
    init.args <- compact(args[init.args])
    
    cluster.args <- discard(names(args), function(x) {x %in% c(input.args, init.args)})
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