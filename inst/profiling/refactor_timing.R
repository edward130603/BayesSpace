library(BayesSpace)
library(mclust)

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
    
    q <- 5
    inputs <- BayesSpace:::.prepare_inputs(sce, use.dimred="PCA")
    init <- mclust::Mclust(inputs$PCs, q, "EEE", verbose=FALSE)$classification
    Y <- as.matrix(inputs$PCs)
    positions <- as.matrix(inputs$positions)
    df_j <- BayesSpace:::find_neighbors(positions, inputs$radius, "manhattan")
    
    for (nrep in c(1000, 2000, 10000)) {
        set.seed(151)
        tres <- system.time(c.refactor <- BayesSpace:::iterate_t_refactor(Y=Y, q=q, df_j=df_j, d=ncol(Y), n=nrow(Y),
                                                 init=init, mu0=colMeans(Y), lambda0=diag(0.01, nrow = ncol(Y)), 
                                                 gamma=2, alpha=1, beta=0.01, nrep=nrep))
        stats.refactor <- data.frame(method="refactor", nspots=nspots, nrep=nrep, user=tres[[1]], system=tres[[2]], elapsed=tres[[3]])
        
        set.seed(151)
        tres <- system.time(c.original <- BayesSpace:::iterate_t_orig(Y=Y, q=q, df_j=df_j, d=ncol(Y), n=nrow(Y),
                                                 init=init, mu0=colMeans(Y), lambda0=diag(0.01, nrow = ncol(Y)), 
                                                 gamma=2, alpha=1, beta=0.01, nrep=nrep))
        stats.original <- data.frame(method="original", nspots=nspots, nrep=nrep, user=tres[[1]], system=tres[[2]], elapsed=tres[[3]])
        
        time_list <- c(time_list, list(stats.refactor), list(stats.original))
    }
}
times <- do.call(rbind, time_list)
write.csv(times, "inst/profiling/refactor_timing.csv", quote=FALSE, row.names=FALSE)
