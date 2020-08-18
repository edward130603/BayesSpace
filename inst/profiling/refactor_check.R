library(mclust)
library(BayesSpace)

set.seed(149)
sce <- exampleSCE(nrow=20, ncol=10)

q <- 5
inputs <- BayesSpace:::.prepare_inputs(sce, use.dimred="PCA")
init <- mclust::Mclust(inputs$PCs, q, "EEE", verbose=FALSE)$classification
Y <- as.matrix(inputs$PCs)
positions <- as.matrix(inputs$positions)
df_j <- BayesSpace:::find_neighbors(positions, inputs$radius, "manhattan")
    
nrep <- 1000
set.seed(151)
tres.refactor <- system.time(c.refactor <- BayesSpace:::iterate_t_refactor(Y=Y, q=q, df_j=df_j, d=ncol(Y), n=nrow(Y),
                                                         init=init, mu0=colMeans(Y), lambda0=diag(0.01, nrow = ncol(Y)), 
                                                         gamma=2, alpha=1, beta=0.01, nrep=nrep, model="t"))
z.refactor <- apply(c.refactor$z[2:1000, ], 2, BayesSpace:::Mode)

set.seed(151)
tres.original <- system.time(c.original <- BayesSpace:::iterate_t_orig(Y=Y, q=q, df_j=df_j, d=ncol(Y), n=nrow(Y),
                                                              init=init, mu0=colMeans(Y), lambda0=diag(0.01, nrow = ncol(Y)), 
                                                              gamma=2, alpha=1, beta=0.01, nrep=nrep))
z.original <- apply(c.original$z[2:1000, ], 2, BayesSpace:::Mode)

if (all.equal(z.refactor, z.original)) {
    print("t: MATCH") 
} else {
    print("t: MISMATCH")
}

print(rbind(tres.refactor, tres.original))

set.seed(151)
tres.refactor <- system.time(c.refactor <- BayesSpace:::iterate_t_refactor(Y=Y, q=q, df_j=df_j, d=ncol(Y), n=nrow(Y),
                                                         init=init, mu0=colMeans(Y), lambda0=diag(0.01, nrow = ncol(Y)), 
                                                         gamma=2, alpha=1, beta=0.01, nrep=nrep, model="normal"))
z.refactor <- apply(c.refactor$z[2:1000, ], 2, BayesSpace:::Mode)

set.seed(151)
tres.original <- system.time(c.original <- BayesSpace:::iterate(Y=Y, q=q, df_j=df_j, d=ncol(Y), n=nrow(Y),
                                                              init=init, mu0=colMeans(Y), lambda0=diag(0.01, nrow = ncol(Y)), 
                                                              gamma=2, alpha=1, beta=0.01, nrep=nrep))
z.original <- apply(c.original$z[2:1000, ], 2, BayesSpace:::Mode)

if (all.equal(z.refactor, z.original)) {
    print("normal: MATCH") 
} else {
    print("normal: MISMATCH")
}

print(rbind(tres.refactor, tres.original))
