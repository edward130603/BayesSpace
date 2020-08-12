library(BayesSpace)

sce <- exampleSCE()
BayesSpace:::start_profiler("inst/profiling/spatialCluster_profile.out")
sce <- spatialCluster(sce, 5)
BayesSpace:::stop_profiler()

## =============================================================================
## Profiling
## =============================================================================
set.seed(149)
sce <- exampleSCE(nrow=20, ncol=20)

q <- 5
inputs <- BayesSpace:::.prepare_inputs(sce, use.dimred="PCA")
init <- mclust::Mclust(inputs$PCs, q, "EEE", verbose=FALSE)$classification
Y <- as.matrix(inputs$PCs)
positions <- as.matrix(inputs$positions)
df_j <- BayesSpace:::find_neighbors(positions, inputs$radius, "manhattan")

for (nrep in c(50000)) {
    BayesSpace:::start_profiler(sprintf("inst/profiling/curr/iterate_t.nrep_%d.out", nrep))
    BayesSpace:::iterate_t(Y=Y, q=q, df_j=df_j, d=ncol(Y), n=nrow(Y),
                              init=init, mu0=colMeans(Y), lambda0=diag(0.01, nrow = ncol(Y)), 
                              gamma=2, alpha=1, beta=0.01, nrep=nrep)
    BayesSpace:::stop_profiler()
}


z.refactor <- apply(c.refactor$z[1000:10000, ], 2, BayesSpace:::Mode)
z.original <- apply(c.original$z[1000:10000, ], 2, BayesSpace:::Mode)

## =============================================================================
## R profiling
## =============================================================================
library(lineprof)
l <- lineprof(spatialCluster(sce, 5, nrep=200))
