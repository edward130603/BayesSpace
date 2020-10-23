suppressPackageStartupMessages({
    library(BayesSpace)
    library(mclust)
})

## =============================================================================
## Prepare inputs
## =============================================================================
prepare_inputs <- function(sce) {
    ## Cluster data
    q <- 5
    inputs <- BayesSpace:::.prepare_inputs(sce, use.dimred="PCA")
    init <- mclust::Mclust(inputs$PCs, q, "EEE", verbose=FALSE)$classification
    Y <- as.matrix(inputs$PCs)
    positions <- as.matrix(inputs$positions)
    df_j <- BayesSpace:::find_neighbors(positions, inputs$radius, "manhattan")
    
    cluster.args <- list(Y=Y, q=q, df_j=df_j, d=ncol(Y), n=nrow(Y),
                         init=init, mu0=colMeans(Y), lambda0=diag(0.01, nrow = ncol(Y)), 
                         gamma=2, alpha=1, beta=0.01)
    
    ## Deconvolution data
    n0 <- nrow(Y)
    d <- ncol(Y)
    c <- 0.3 * 1 / (2 * mean(diag(cov(Y))))
    subspots <- 9
    
    init1 <- rep(init, subspots)
    Y2 <- Y[rep(seq_len(n0), subspots), ]
    positions2 <- positions[rep(seq_len(n0), subspots), ]
    
    shift <- BayesSpace:::.make_subspot_offsets(subspots)
    shift <- t(t(shift) * c(inputs$xdist, inputs$ydist))
    dist <- (max(rowSums(abs(shift))) * 1.05) / 2
    
    shift_long <- shift[rep(seq_len(subspots), each=n0), ]
    positions2[, "x"] <- positions2[, "x"] + shift_long[, "Var1"]
    positions2[, "y"] <- positions2[, "y"] + shift_long[, "Var2"]
    n <- nrow(Y2)
    
    df_j2 <- BayesSpace:::find_neighbors(positions2, dist, "manhattan")
    
    deconv.args <- list(Y=Y2, df_j=df_j2, tdist=TRUE, n=nrow(Y2), n0=nrow(Y),
                 d=ncol(Y), gamma=2, q=q, init=init1, subspots=subspots, verbose=FALSE,
                 jitter_scale=5, c=c, mu0=colMeans(Y), lambda0=diag(0.01, nrow=ncol(Y)),
                 alpha=1, beta=0.01)
    
    return(list(cluster.args=cluster.args, deconv.args=deconv.args))
}

## =============================================================================
## Timing
## =============================================================================
time_list <- list()
for (nspots in c(100, 200, 1000)) {
    set.seed(149)
    if (nspots == 100) {
        sce <- exampleSCE(nrow=10, ncol=10)
    } else if (nspots == 200) {
        sce <- exampleSCE(nrow=20, ncol=10)
    } else {
        sce <- exampleSCE(nrow=20, ncol=50)
    }
    
    inputs <- prepare_inputs(sce)
    cluster.args <- inputs$cluster.args
    deconv.args <- inputs$deconv.args
    
    for (nrep in c(1000, 2000, 10000)) {
        set.seed(151)
        tres <- system.time(c.refactor <- do.call(BayesSpace:::iterate_t_refactor, c(cluster.args, list(nrep=nrep))))
        cluster.refactor <- data.frame(method="refactor", algorithm="cluster", nspots=nspots, nrep=nrep, user=tres[[1]], system=tres[[2]], elapsed=tres[[3]])
        
        set.seed(151)
        tres <- system.time(c.original <- do.call(BayesSpace:::iterate_t, c(cluster.args, list(nrep=nrep))))
        cluster.original <- data.frame(method="original", algorithm="cluster", nspots=nspots, nrep=nrep, user=tres[[1]], system=tres[[2]], elapsed=tres[[3]])
        
        set.seed(151)
        tres <- system.time(c.refactor <- do.call(BayesSpace:::iterate_deconv_refactor, c(deconv.args, list(nrep=nrep))))
        deconv.refactor <- data.frame(method="refactor", algorithm="deconvolve", nspots=nspots, nrep=nrep, user=tres[[1]], system=tres[[2]], elapsed=tres[[3]])
        
        set.seed(151)
        tres <- system.time(c.original <- do.call(BayesSpace:::iterate_deconv, c(deconv.args, list(nrep=nrep))))
        deconv.original <- data.frame(method="original", algorithm="deconvolve", nspots=nspots, nrep=nrep, user=tres[[1]], system=tres[[2]], elapsed=tres[[3]])
        
        time_list <- c(time_list, list(cluster.refactor), list(cluster.original), list(deconv.refactor), list(deconv.original))
    }
}
times <- do.call(rbind, time_list)
write.csv(times, "inst/profiling/refactor_timing.csv", quote=FALSE, row.names=FALSE)
