suppressPackageStartupMessages({
    library(mclust)
    library(BayesSpace)
})

## Helper to compare refactored/original results
check_refactor <- function(fn_refactor, fn_orig, args, name) {
    set.seed(151)
    c.refactor <- do.call(fn_refactor, args)
    z.refactor <- apply(c.refactor$z, 2, BayesSpace:::Mode)
    
    set.seed(151)
    c.original <- do.call(fn_orig, args)
    z.original <- apply(c.original$z, 2, BayesSpace:::Mode)
    
    if (all.equal(z.refactor, z.original)) {
        print(paste0(name, ": MATCH"))
    } else {
        print(paste0(name, "t: MISMATCH"))
    }
}

## =============================================================================
## Clustering
## =============================================================================

## Example data
set.seed(149)
sce <- exampleSCE(nrow=20, ncol=10)
q <- 5
inputs <- BayesSpace:::.prepare_inputs(sce, use.dimred="PCA")
init <- mclust::Mclust(inputs$PCs, q, "EEE", verbose=FALSE)$classification
Y <- as.matrix(inputs$PCs)
positions <- as.matrix(inputs$positions)
df_j <- BayesSpace:::find_neighbors(positions, inputs$radius, "manhattan")
nrep <- 1000

## Equal precision
args <- list(Y=Y, q=q, df_j=df_j, d=ncol(Y), n=nrow(Y),
             init=init, mu0=colMeans(Y), lambda0=diag(0.01, nrow = ncol(Y)), 
             gamma=2, alpha=1, beta=0.01, nrep=nrep)

check_refactor(BayesSpace:::iterate_t_refactor, BayesSpace:::iterate_t, args, "t equal")
check_refactor(BayesSpace:::iterate_refactor, BayesSpace:::iterate, args, "normal equal")

## Variable precision (need to bump alpha/beta for numerical instability)
args <- list(Y=Y, q=q, df_j=df_j, d=ncol(Y), n=nrow(Y),
             init=init, mu0=colMeans(Y), lambda0=diag(0.01, nrow = ncol(Y)), 
             gamma=2, alpha=10, beta=0.11, nrep=nrep)

check_refactor(BayesSpace:::iterate_t_vvv_refactor, BayesSpace:::iterate_t_vvv, args, "t variable")
check_refactor(BayesSpace:::iterate_vvv_refactor, BayesSpace:::iterate_vvv, args, "normal variable")

## =============================================================================
## Deconvolution
## =============================================================================
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

args <- list(Y=Y2, df_j=df_j2, tdist=TRUE, nrep=nrep, n=nrow(Y2), n0=nrow(Y),
             d=ncol(Y), gamma=2, q=q, init=init1, subspots=subspots, verbose=FALSE,
             jitter_scale=5, c=c, mu0=colMeans(Y), lambda0=diag(0.01, nrow=ncol(Y)),
             alpha=1, beta=0.01)

check_refactor(BayesSpace:::iterate_deconv_refactor, BayesSpace:::iterate_deconv, args, "deconvolution")
